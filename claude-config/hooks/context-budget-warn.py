#!/usr/bin/env python3
"""context-budget-warn.py — Stop hook.

Non-blocking. Reports when the session's context has grown past a band worth
acting on, and says so once per band.

Context size is the dominant cost term in a long session: every request re-reads
the whole prefix, so a session running at 400K is paying 400K on each of
hundreds of turns. The number is not visible anywhere during normal work — this
hook surfaces it at a natural pause.

The measurement is the last `cache_read_input_tokens` in the transcript, which is
the size of the prefix the API actually re-read on the most recent request. That
is the real number, not an estimate.

Cost discipline: the transcript is read from the tail (last TAIL_BYTES only) and
scanned backwards for the first usage record, so this stays cheap on the very
large transcripts it is most relevant to.

Suppress for a session by creating the state file with {"muted": true}.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

# Bands in tokens. Warn once when context first exceeds each.
BANDS = (150_000, 250_000, 350_000)

# Tail window. Comfortably larger than one transcript line, small enough that
# reading it is trivial next to the context being reported on.
TAIL_BYTES = 262_144

STATE_TTL_DAYS = 30

ADVICE = {
    150_000: (
        "Finish the current step, then consider /clear before starting anything\n"
        "  new. Delegate the next broad sweep instead of reading files inline."
    ),
    250_000: (
        "This is expensive per turn. Wrap up and /clear; carry forward only what\n"
        "  the next task needs. Long sessions are the single largest cost driver."
    ),
    350_000: (
        "Near the practical limit. Every remaining turn re-reads all of this.\n"
        "  Stop here, write down state, and /clear."
    ),
}


def _state_path(session_id: str) -> Path:
    base = Path(os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd()))
    state_dir = base / ".claude" / "session-state"
    state_dir.mkdir(parents=True, exist_ok=True)
    safe_id = re.sub(r"[^A-Za-z0-9_\-]", "_", session_id)[:64] or "default"
    return state_dir / f"{safe_id}-ctxbudget.json"


def _cleanup_old_state(state_dir: Path) -> None:
    """Reap ctxbudget state files older than STATE_TTL_DAYS. Best-effort."""
    if not state_dir.exists():
        return
    cutoff = time.time() - STATE_TTL_DAYS * 86400
    try:
        for entry in state_dir.iterdir():
            if not entry.name.endswith("-ctxbudget.json"):
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


def _last_context_size(transcript_path: str) -> int:
    """Largest context size on the last request that reported one.

    Reads only the tail of the transcript and scans backwards. Returns 0 when
    the transcript is missing, unreadable, or carries no usage records.
    """
    if not transcript_path:
        return 0
    p = Path(transcript_path)
    try:
        size = p.stat().st_size
        with p.open("rb") as f:
            if size > TAIL_BYTES:
                f.seek(size - TAIL_BYTES)
                f.readline()  # discard the partial line
            chunk = f.read()
    except OSError:
        return 0

    for line in reversed(chunk.decode("utf-8", errors="ignore").splitlines()):
        line = line.strip()
        if not line or "cache_read_input_tokens" not in line:
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        msg = event.get("message")
        usage = msg.get("usage") if isinstance(msg, dict) else None
        if not isinstance(usage, dict):
            continue
        total = 0
        for key in ("cache_read_input_tokens", "cache_creation_input_tokens", "input_tokens"):
            value = usage.get(key)
            if isinstance(value, int):
                total += value
        if total:
            return total
    return 0


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0

    path = _state_path(payload.get("session_id", "") or "")
    _cleanup_old_state(path.parent)

    state = _load(path)
    if state.get("muted"):
        return 0

    size = _last_context_size(payload.get("transcript_path", "") or "")
    if not size:
        return 0

    band = max((b for b in BANDS if size >= b), default=0)
    if not band or band <= int(state.get("warned", 0)):
        return 0

    _save(path, {"warned": band, "muted": False})
    print(
        f"\n[context-budget-warn] Context is ~{size // 1000}K tokens, re-read on\n"
        f"  every request from here.\n"
        f"  {ADVICE[band]}\n",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
