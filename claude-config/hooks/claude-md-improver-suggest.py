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
import re
import sys
from pathlib import Path

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


def _walk_user_messages(transcript_path: str) -> list[tuple[str, str]]:
    """Return list of (matched_phrase, surrounding_context) for user
    messages containing a trigger phrase."""
    hits: list[tuple[str, str]] = []
    if not transcript_path or not Path(transcript_path).exists():
        return hits
    try:
        with open(transcript_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                # Only user-role messages
                msg = event.get("message", event)
                if not isinstance(msg, dict):
                    continue
                role = msg.get("role") or event.get("role")
                if role != "user":
                    continue
                content = msg.get("content", "")
                # content may be string or list of dicts
                if isinstance(content, list):
                    text = " ".join(
                        c.get("text", "") for c in content
                        if isinstance(c, dict) and c.get("type") == "text"
                    )
                else:
                    text = str(content)
                if not text:
                    continue
                for pat in TRIGGER_PATTERNS:
                    m = pat.search(text)
                    if not m:
                        continue
                    start = max(0, m.start() - 30)
                    end = min(len(text), m.end() + 60)
                    excerpt = text[start:end].replace("\n", " ").strip()
                    hits.append((m.group(0), excerpt))
                    break  # one hit per message is enough
    except OSError:
        pass
    return hits


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    transcript_path = payload.get("transcript_path", "")
    hits = _walk_user_messages(transcript_path)
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
