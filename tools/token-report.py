#!/usr/bin/env python3
"""token-report.py — where the tokens actually go.

Reads per-request usage straight out of Claude Code's session transcripts
(`~/.claude/projects/**/*.jsonl`) and reports the numbers that drive cost.

Why not `~/.claude/stats-cache.json`: it is aggregate-only, recomputed
periodically (it was 12 days stale when this analysis started), and cannot show
per-request context size — which is the number that matters. The transcripts
carry `usage` on every assistant message, so they are the ground truth.

The headline metric is **cache read per request**: the size of the prompt prefix
the API re-read on each request. Cache reads run ~99% of token volume and ~88%
of cost, so context size x request count is the whole game. Output tokens are a
rounding error by comparison.

Also reports the behavioural metrics that the 2026-08-10 policy change targets:
how often skills are invoked, and how often work is delegated to a subagent
instead of being done in the main loop.

Usage:
    ./tools/token-report.py                    # full report vs. recorded baseline
    ./tools/token-report.py --top 15           # more per-session detail
    ./tools/token-report.py --since 2026-08-11 # only sessions modified since
    ./tools/token-report.py --json             # machine-readable
    ./tools/token-report.py --projects-dir DIR # alternate transcript root

Stdlib only, matching the repo's dependency floor (PyYAML is for sync-config.py).
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Iterable

# Measured 2026-08-10 over 311 sessions / 34,889 requests, before the
# user-level CLAUDE.md policy and the two context-hygiene hooks landed.
# Deltas are shown against this so a later run says whether it worked.
BASELINE = {
    "label": "2026-08-10 (pre-policy)",
    "sessions": 311,
    "requests": 34889,
    "cache_read_per_request": 163573,
    "output_per_request": 852,
    "pct_sessions_with_skill": 15.1,
    "skill_calls": 75,
    "agent_dispatches": 262,
}

USAGE_KEYS = (
    "input_tokens",
    "output_tokens",
    "cache_read_input_tokens",
    "cache_creation_input_tokens",
)
DISPATCH_TOOLS = ("Task", "Agent")


class SessionStats:
    """Aggregated usage and behaviour for one transcript file."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.requests = 0
        self.totals: Counter[str] = Counter()
        self.skills: Counter[str] = Counter()
        self.dispatches: Counter[str] = Counter()

    @property
    def cache_read(self) -> int:
        return self.totals["cache_read_input_tokens"]

    @property
    def per_request(self) -> float:
        return self.cache_read / self.requests if self.requests else 0.0


def iter_events(path: Path) -> Iterable[dict]:
    """Yield parsed JSON events, skipping blank and malformed lines."""
    try:
        with path.open(errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue
    except OSError:
        return


def scan(path: Path) -> SessionStats:
    stats = SessionStats(path)
    for event in iter_events(path):
        message = event.get("message")
        if not isinstance(message, dict):
            continue

        usage = message.get("usage")
        if isinstance(usage, dict):
            stats.requests += 1
            for key in USAGE_KEYS:
                value = usage.get(key)
                if isinstance(value, int):
                    stats.totals[key] += value

        content = message.get("content")
        if not isinstance(content, list):
            continue
        for item in content:
            if not isinstance(item, dict) or item.get("type") != "tool_use":
                continue
            name = item.get("name", "")
            tool_input = item.get("input") or {}
            if name == "Skill":
                stats.skills[tool_input.get("skill", "?")] += 1
            elif name in DISPATCH_TOOLS:
                stats.dispatches[tool_input.get("subagent_type", "(default)")] += 1
    return stats


def collect(projects_dir: Path, since: datetime | None) -> list[SessionStats]:
    out = []
    for path in sorted(projects_dir.rglob("*.jsonl")):
        if since is not None:
            try:
                if datetime.fromtimestamp(path.stat().st_mtime) < since:
                    continue
            except OSError:
                continue
        stats = scan(path)
        if stats.requests:
            out.append(stats)
    return out


def _delta(current: float, base: float, *, lower_is_better: bool = True) -> str:
    """Human-readable delta against the baseline."""
    if not base:
        return ""
    pct = 100.0 * (current - base) / base
    if abs(pct) < 0.5:
        return "  (flat vs baseline)"
    better = (pct < 0) if lower_is_better else (pct > 0)
    return f"  ({pct:+.0f}% vs baseline — {'better' if better else 'worse'})"


def build_summary(sessions: list[SessionStats]) -> dict:
    totals: Counter[str] = Counter()
    skills: Counter[str] = Counter()
    dispatches: Counter[str] = Counter()
    requests = 0
    with_skill = 0
    for s in sessions:
        requests += s.requests
        totals.update(s.totals)
        skills.update(s.skills)
        dispatches.update(s.dispatches)
        if s.skills:
            with_skill += 1
    n = len(sessions) or 1
    r = requests or 1
    return {
        "sessions": len(sessions),
        "requests": requests,
        "totals": dict(totals),
        "cache_read_per_request": totals["cache_read_input_tokens"] / r,
        "output_per_request": totals["output_tokens"] / r,
        "pct_sessions_with_skill": 100.0 * with_skill / n,
        "skill_calls": sum(skills.values()),
        "agent_dispatches": sum(dispatches.values()),
        "skills": dict(skills),
        "dispatches": dict(dispatches),
    }


def render(sessions: list[SessionStats], summary: dict, top: int) -> None:
    t = summary["totals"]
    volume = sum(t.get(k, 0) for k in USAGE_KEYS) or 1

    print(f"sessions {summary['sessions']}   requests {summary['requests']:,}")
    print(f"baseline: {BASELINE['label']}\n")

    print("token volume")
    for key in ("cache_read_input_tokens", "cache_creation_input_tokens",
                "output_tokens", "input_tokens"):
        value = t.get(key, 0)
        print(f"  {key:32} {value:>16,}  {100.0 * value / volume:5.1f}%")

    cr = summary["cache_read_per_request"]
    op = summary["output_per_request"]
    print("\nper request  <- the number that matters")
    print(f"  cache read (context size)        {cr:>12,.0f}"
          f"{_delta(cr, BASELINE['cache_read_per_request'])}")
    print(f"  output                           {op:>12,.0f}"
          f"{_delta(op, BASELINE['output_per_request'])}")
    if op:
        print(f"  ratio cache-read : output        {cr / op:>12,.0f}x")

    print("\ndelegation and skills")
    pct = summary["pct_sessions_with_skill"]
    print(f"  sessions invoking >=1 skill      {pct:>11.1f}%"
          f"{_delta(pct, BASELINE['pct_sessions_with_skill'], lower_is_better=False)}")
    print(f"  skill invocations                {summary['skill_calls']:>12,}")
    print(f"  agent dispatches                {summary['agent_dispatches']:>12,}")
    gp = summary["dispatches"].get("general-purpose", 0)
    if summary["agent_dispatches"]:
        share = 100.0 * gp / summary["agent_dispatches"]
        print(f"  ... general-purpose             {gp:>12,}  ({share:.1f}%"
              f" — inherits session model unless model= is passed)")

    if summary["skills"]:
        print("\nskills invoked")
        for name, count in Counter(summary["skills"]).most_common(12):
            print(f"  {name:44} {count:5}")

    print(f"\ntop {top} sessions by context re-read")
    ranked = sorted(sessions, key=lambda s: -s.cache_read)[:top]
    for s in ranked:
        print(f"  {s.cache_read / 1e6:8.1f}M  {s.requests:5} reqs  "
              f"{s.per_request / 1000:6.0f}K/req  {s.path.parent.name[:46]}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Report where Claude Code tokens go, from session transcripts."
    )
    parser.add_argument("--projects-dir", type=Path,
                        default=Path.home() / ".claude" / "projects")
    parser.add_argument("--top", type=int, default=8,
                        help="how many sessions to list (default 8)")
    parser.add_argument("--since", metavar="YYYY-MM-DD",
                        help="only transcripts modified on/after this date")
    parser.add_argument("--json", action="store_true", dest="as_json",
                        help="emit the summary as JSON")
    args = parser.parse_args(argv)

    if not args.projects_dir.is_dir():
        print(f"no transcript directory at {args.projects_dir}", file=sys.stderr)
        return 1

    since = None
    if args.since:
        try:
            since = datetime.strptime(args.since, "%Y-%m-%d")
        except ValueError:
            print(f"--since expects YYYY-MM-DD, got {args.since!r}", file=sys.stderr)
            return 1

    sessions = collect(args.projects_dir, since)
    if not sessions:
        print("no transcripts with usage data found", file=sys.stderr)
        return 1

    summary = build_summary(sessions)
    if args.as_json:
        print(json.dumps({"baseline": BASELINE, "current": summary}, indent=2))
    else:
        render(sessions, summary, args.top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
