#!/usr/bin/env python3
"""session-start-context.py — UserPromptSubmit hook, fires once per session.

On the first user prompt of a session, emits a compact context block to
stdout: current working directory, git branch + ahead/behind, last 2
commits, dirty-file count, and seconds since the last commit. Claude
Code prepends UserPromptSubmit stdout to the conversation context so
the model sees this without explicit explanation.

When run inside the Claude config repo (detected by a `sync-config.py`
at cwd), it also appends a one-line drift warning if `~/.claude`
diverges from `claude-config/` — surfacing untracked live edits before
they get clobbered by the next push. Best-effort and fail-silent.

Subsequent prompts in the same session are no-ops (a marker file in
the per-session state dir ensures one-shot behavior).

Exit code 0 always — informational only, never blocks.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path


def _git(args: list[str], cwd: str | None = None) -> str:
    try:
        return subprocess.check_output(
            ["git"] + args,
            text=True,
            stderr=subprocess.DEVNULL,
            cwd=cwd,
            timeout=2,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return ""


MARKER_TTL_DAYS = 30


def _marker_path(session_id: str) -> Path:
    base = Path(os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd()))
    marker_dir = base / ".claude" / "session-state"
    marker_dir.mkdir(parents=True, exist_ok=True)
    safe_id = re.sub(r"[^A-Za-z0-9_\-]", "_", session_id)[:64] or "default"
    return marker_dir / f"{safe_id}-start.marker"


def _cleanup_old_markers(marker_dir: Path) -> None:
    """Reap session-state marker files older than MARKER_TTL_DAYS.

    Best-effort: silently skip files that can't be stat'd or unlinked.
    Runs on every session start so the dir stays bounded without any
    external cron job.
    """
    import time
    if not marker_dir.exists():
        return
    cutoff = time.time() - MARKER_TTL_DAYS * 86400
    try:
        for entry in marker_dir.iterdir():
            if not entry.name.endswith("-start.marker"):
                continue
            try:
                if entry.stat().st_mtime < cutoff:
                    entry.unlink()
            except OSError:
                continue
    except OSError:
        return


def _build_context(cwd: str) -> str | None:
    if not _git(["rev-parse", "--is-inside-work-tree"], cwd=cwd):
        return None
    branch = _git(["branch", "--show-current"], cwd=cwd) or "DETACHED"
    upstream = _git(["rev-parse", "--abbrev-ref", "@{upstream}"], cwd=cwd) or "(no upstream)"
    counts = _git(["rev-list", "--left-right", "--count", "@{upstream}...HEAD"], cwd=cwd)
    ahead_behind = ""
    if counts:
        parts = counts.split()
        if len(parts) == 2:
            behind, ahead = parts
            if int(ahead) or int(behind):
                ahead_behind = f"  ahead {ahead}, behind {behind}"
    dirty = _git(["status", "--porcelain"], cwd=cwd).splitlines()
    last_commits = _git(["log", "--oneline", "-2"], cwd=cwd)
    last_commit_ts = _git(["log", "-1", "--format=%ct"], cwd=cwd)
    seconds_since = ""
    if last_commit_ts.isdigit():
        import time
        delta = int(time.time()) - int(last_commit_ts)
        if delta < 60:
            seconds_since = f"{delta}s ago"
        elif delta < 3600:
            seconds_since = f"{delta // 60}m ago"
        elif delta < 86400:
            seconds_since = f"{delta // 3600}h ago"
        else:
            seconds_since = f"{delta // 86400}d ago"

    lines = [
        "<session-context>",
        f"cwd: {cwd}",
        f"branch: {branch} (upstream: {upstream}){ahead_behind}",
        f"dirty files: {len(dirty)}" + (f" ({dirty[0].strip()[:60]}...)" if dirty else ""),
        f"last commit: {seconds_since}" if seconds_since else "no commits yet",
    ]
    if last_commits:
        lines.append("recent commits:")
        for c in last_commits.splitlines():
            lines.append(f"  {c}")
    drift = _drift_warning(cwd)
    if drift:
        lines.append(drift)
    lines.append("</session-context>")
    return "\n".join(lines)


def _drift_warning(cwd: str) -> str | None:
    """If cwd is the Claude config repo and ~/.claude diverges from
    claude-config/, return a one-line warning; else None.

    Self-scoping: only runs when a `sync-config.py` sits at cwd (true
    only for this repo). Best-effort — any error/timeout returns None so
    the session-context block is never blocked or delayed past 15s.
    """
    sync = Path(cwd) / "sync-config.py"
    if not sync.is_file():
        return None
    try:
        out = subprocess.run(
            ["python3", str(sync), "status"],
            capture_output=True,
            text=True,
            timeout=15,
            cwd=cwd,
        ).stdout
    except (subprocess.SubprocessError, OSError):
        return None
    # Only inspect the user-wide section (before project configs); strip ANSI.
    section = re.sub(r"\x1b\[[0-9;]*m", "", out.split("Project configurations")[0])
    changes = [ln.strip() for ln in section.splitlines() if re.match(r"\s*\*\s", ln)]
    m = re.search(r"Orphaned files \((\d+)\)", section)
    orphans = int(m.group(1)) if m else 0
    if not changes and not orphans:
        return None
    bits = list(changes)
    if orphans:
        bits.append(f"{orphans} orphan(s)")
    return (
        "⚠ ~/.claude diverges from this repo ["
        + "; ".join(bits)
        + "] — run ./sync-config.py status; pull live edits or push repo "
        "changes before editing config (repo is source of truth)."
    )


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0

    session_id = payload.get("session_id", "")
    marker = _marker_path(session_id)
    # Reap stale markers before deciding whether to fire. Cheap (one
    # iterdir + stat per file), bounded by MARKER_TTL_DAYS.
    _cleanup_old_markers(marker.parent)
    if marker.exists():
        return 0  # already fired this session

    cwd = payload.get("cwd") or os.getcwd()
    block = _build_context(cwd)
    if block:
        print(block)  # goes into the model's view

    try:
        marker.touch()
    except OSError:
        pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
