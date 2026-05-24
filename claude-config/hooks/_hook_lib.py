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
from pathlib import Path
from typing import Iterator


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
