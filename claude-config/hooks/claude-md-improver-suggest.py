#!/usr/bin/env python3
"""claude-md-improver-suggest.py — Stop hook.

Scans the session transcript for user messages containing phrases that
suggest a CLAUDE.md update would be useful (durable preferences,
"from now on", "next time", "remember to", etc.). If found, emits a
gentle reminder to consider invoking the `claude-md-management:
revise-claude-md` skill before ending the session.

Non-blocking — exit 0 always.
"""
from __future__ import annotations

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _hook_lib import iter_user_messages  # noqa: E402  # pyright: ignore[reportMissingImports]

TRIGGER_PATTERNS = [
    re.compile(r"\bfrom now on\b", re.IGNORECASE),
    re.compile(r"\bnext time\b", re.IGNORECASE),
    re.compile(r"\bremember to\b", re.IGNORECASE),
    re.compile(r"\balways (?:use|prefer|do|skip|avoid)\b", re.IGNORECASE),
    re.compile(r"\bnever (?:use|do|commit|push)\b", re.IGNORECASE),
    re.compile(r"\bdon't (?:use|do|forget)\b", re.IGNORECASE),
    re.compile(r"\bplease (?:always|never|don't)\b", re.IGNORECASE),
    re.compile(r"\bgoing forward\b", re.IGNORECASE),
    re.compile(r"\bmake (?:a|this a) (?:rule|habit|convention)\b", re.IGNORECASE),
]


def _scan_user_messages(transcript_path: str) -> list[tuple[str, str]]:
    """Return list of (matched_phrase, surrounding_context) for user
    messages containing a trigger phrase."""
    hits: list[tuple[str, str]] = []
    for text in iter_user_messages(transcript_path):
        for pat in TRIGGER_PATTERNS:
            m = pat.search(text)
            if not m:
                continue
            start = max(0, m.start() - 30)
            end = min(len(text), m.end() + 60)
            excerpt = text[start:end].replace("\n", " ").strip()
            hits.append((m.group(0), excerpt))
            break  # one hit per message is enough
    return hits


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    transcript_path = payload.get("transcript_path", "")
    hits = _scan_user_messages(transcript_path)
    if not hits:
        return 0

    # Dedupe by excerpt; keep at most 3
    seen: set[str] = set()
    unique: list[tuple[str, str]] = []
    for phrase, excerpt in hits:
        if excerpt in seen:
            continue
        seen.add(excerpt)
        unique.append((phrase, excerpt))
        if len(unique) >= 3:
            break

    lines = [
        "[claude-md-improver-suggest] This session contained durable-preference phrases:",
    ]
    for phrase, excerpt in unique:
        lines.append(f"  - matched '{phrase}': \"...{excerpt}...\"")
    lines.append("")
    lines.append("  Consider invoking `claude-md-management:revise-claude-md` to capture")
    lines.append("  any new conventions in CLAUDE.md so they persist across sessions.")
    print("\n".join(lines), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
