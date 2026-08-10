#!/usr/bin/env python3
"""read-burst-nudge.py — PostToolUse hook (Read, Grep, Glob).

Non-blocking. Fires once per session when the main thread has done a lot of
file reading without ever dispatching a subagent.

Everything read into the main context is re-read on every later request in the
session. A twenty-file sweep done inline is charged on every turn that follows;
the same sweep in a subagent is charged once and comes back as a summary. This
hook exists because that cost is invisible at the moment it is incurred.

Cost discipline: the read counter is an O(1) increment against a small JSON
state file, so the common path does no transcript I/O. The transcript is walked
at most once per session — when the counter first crosses the threshold — to
check whether any subagent has been dispatched. After that the session is
marked done and the hook returns immediately forever.

Subagents are exempt: a subagent reading many files IS the recommended pattern.

Suppress with `noqa-read-burst` in the tool call's description.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _hook_lib import iter_tool_uses  # noqa: E402  # pyright: ignore[reportMissingImports]

WATCHED_TOOLS = ("Read", "Grep", "Glob")
DISPATCH_TOOLS = ("Task", "Agent")

# Reads allowed in the main thread before the nudge. Set from observed
# behaviour: sessions that stay under this rarely reach problem context sizes.
THRESHOLD = 15

STATE_TTL_DAYS = 30

MESSAGE = (
    "\n[read-burst-nudge] {n} file reads in the main context this session, with\n"
    "  no subagent dispatched. Each one is re-read on every later request, so a\n"
    "  broad sweep done inline compounds for the rest of the session.\n"
    "  For the next sweep, dispatch a subagent and have it return findings —\n"
    "  pass model=\"sonnet\" (or \"haiku\" for a mechanical scan) so it does not\n"
    "  inherit Opus. Keep the decision and the edits here.\n"
    "  Also: prefer Read with offset/limit and Grep -A/-B over whole files.\n"
    "  Suppress: add `noqa-read-burst` to the tool call's description.\n"
)


def _state_path(session_id: str) -> Path:
    base = Path(os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd()))
    state_dir = base / ".claude" / "session-state"
    state_dir.mkdir(parents=True, exist_ok=True)
    safe_id = re.sub(r"[^A-Za-z0-9_\-]", "_", session_id)[:64] or "default"
    return state_dir / f"{safe_id}-readburst.json"


def _cleanup_old_state(state_dir: Path) -> None:
    """Reap readburst state files older than STATE_TTL_DAYS.

    Best-effort, same contract as session-start-context.py's marker reaping:
    silently skip anything that cannot be stat'd or unlinked.
    """
    if not state_dir.exists():
        return
    cutoff = time.time() - STATE_TTL_DAYS * 86400
    try:
        for entry in state_dir.iterdir():
            if not entry.name.endswith("-readburst.json"):
                continue
            try:
                if entry.stat().st_mtime < cutoff:
                    entry.unlink()
            except OSError:
                continue
    except OSError:
        return


def _load(path: Path) -> dict:
    try:
        data = json.loads(path.read_text())
        return data if isinstance(data, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _save(path: Path, state: dict) -> None:
    try:
        path.write_text(json.dumps(state))
    except OSError:
        pass


def _is_subagent(payload: dict) -> bool:
    """True when this hook fires inside a subagent rather than the main thread.

    Mirrors paper-context-guard.py. Falls back to False, since a spurious
    nudge in a subagent is cheap and a missed one in the main thread is not.
    """
    for key in ("is_subagent", "isSubagent", "subagent"):
        if payload.get(key):
            return True
    if payload.get("agent_type") not in (None, "", "main"):
        return True
    return bool(os.environ.get("CLAUDE_AGENT_TYPE"))


def _has_dispatched(transcript_path: str) -> bool:
    for name, _ in iter_tool_uses(transcript_path):
        if name in DISPATCH_TOOLS:
            return True
    return False


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0

    if payload.get("tool_name") not in WATCHED_TOOLS:
        return 0

    tool_input = payload.get("tool_input", {}) or {}
    if "noqa-read-burst" in (tool_input.get("description", "") or ""):
        return 0
    if _is_subagent(payload):
        return 0

    path = _state_path(payload.get("session_id", "") or "")
    _cleanup_old_state(path.parent)

    state = _load(path)
    if state.get("done"):
        return 0

    count = int(state.get("reads", 0)) + 1
    if count < THRESHOLD:
        _save(path, {"reads": count, "done": False})
        return 0

    # Threshold reached: this is the one and only transcript walk. Whatever we
    # find, the session is done being watched.
    _save(path, {"reads": count, "done": True})
    if _has_dispatched(payload.get("transcript_path", "") or ""):
        return 0

    print(MESSAGE.format(n=count), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
