"""Shared helpers for hooks under ~/.claude/hooks/.

Currently provides transcript-walking utilities. Each hook is otherwise
self-contained — this is the one piece of common infrastructure. If
more shared behavior accumulates, it lives here.

Import pattern from a hook script:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _hook_lib import iter_tool_uses, iter_user_messages
"""
from __future__ import annotations

import json
import sys as _sys
from pathlib import Path
from typing import Iterator


# ----- Diagnose mode (used by blocking hooks for debugging false positives) -----


def is_diagnose_mode() -> bool:
    """True if --diagnose appears in this script's argv.

    Blocking hooks check this to convert a "would block" decision into
    an "informational print + exit 0" without actually interfering with
    the tool call. Use by piping a synthetic JSON payload into the
    hook with --diagnose appended.
    """
    return "--diagnose" in _sys.argv[1:]


def report(message: str, *, block: bool = True) -> int:
    """Emit a hook decision and return the appropriate exit code.

    In normal mode:
      - block=True:  prints to stderr, returns 2 (Claude Code blocks)
      - block=False: prints to stderr, returns 0 (informational only)

    In diagnose mode (--diagnose in argv):
      - prints to stdout with a clear DIAGNOSE prefix, returns 0
        regardless of block flag. This lets the user inspect what the
        hook would do without it actually intervening.
    """
    if is_diagnose_mode():
        prefix = "[DIAGNOSE] would block:\n" if block else "[DIAGNOSE] would warn:\n"
        print(prefix + message)
        return 0
    print(message, file=_sys.stderr)
    return 2 if block else 0


def iter_transcript_events(transcript_path: str) -> Iterator[dict]:
    """Yield one parsed JSON event per line of the transcript JSONL.

    Silently skips empty lines and malformed JSON entries. Yields the
    raw event dict — callers navigate to message.content as needed.
    """
    if not transcript_path:
        return
    p = Path(transcript_path)
    if not p.exists():
        return
    try:
        with p.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue
    except OSError:
        return


def iter_tool_uses(transcript_path: str) -> Iterator[tuple[str, dict]]:
    """Yield (tool_name, tool_input) for every tool_use entry across
    every assistant message in the transcript.

    Handles both legacy (top-level fields) and v2 (`message.content[]`)
    transcript shapes.
    """
    for event in iter_transcript_events(transcript_path):
        msg = event.get("message", event)
        if not isinstance(msg, dict):
            continue
        content = msg.get("content")
        if not isinstance(content, list):
            continue
        for item in content:
            if not isinstance(item, dict) or item.get("type") != "tool_use":
                continue
            name = item.get("name", "")
            inp = item.get("input", {}) or {}
            yield name, inp


def iter_user_messages(transcript_path: str) -> Iterator[str]:
    """Yield the text body of every user-role message.

    For messages whose content is a list of blocks (text + tool_use +
    tool_result), only the `text` blocks are joined into the returned
    string. Empty messages are skipped.
    """
    for event in iter_transcript_events(transcript_path):
        msg = event.get("message", event)
        if not isinstance(msg, dict):
            continue
        role = msg.get("role") or event.get("role")
        if role != "user":
            continue
        content = msg.get("content", "")
        if isinstance(content, list):
            text = " ".join(
                c.get("text", "") for c in content
                if isinstance(c, dict) and c.get("type") == "text"
            )
        else:
            text = str(content)
        if text:
            yield text
