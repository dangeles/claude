#!/usr/bin/env python3
"""stale-todo-tracker.py — Stop hook.

At session end, scan every Write/Edit performed in the conversation
and emit a summary of TODO/FIXME/XXX/HACK markers added or removed.
Non-blocking — never refuses to stop.

Detects additions in Write content and Edit new_string; detects
removals in Edit old_string. The report shows net change so I can
see at a glance whether the session left more debt than it cleared.

Exit code 0 always.
"""
from __future__ import annotations

import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _hook_lib import iter_tool_uses  # noqa: E402  # pyright: ignore[reportMissingImports]

MARKER = re.compile(r"\b(TODO|FIXME|XXX|HACK)\b[^\n]{0,120}")


def _walk_transcript(
    transcript_path: str,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Return (added, removed) — each a list of (file_path, marker_line)."""
    added: list[tuple[str, str]] = []
    removed: list[tuple[str, str]] = []
    for name, inp in iter_tool_uses(transcript_path):
        path = inp.get("file_path") or inp.get("notebook_path") or ""
        if not path:
            continue
        if name == "Write":
            for m in MARKER.finditer(inp.get("content", "") or ""):
                added.append((path, m.group(0).strip()))
        elif name == "Edit":
            for m in MARKER.finditer(inp.get("new_string", "") or ""):
                added.append((path, m.group(0).strip()))
            for m in MARKER.finditer(inp.get("old_string", "") or ""):
                removed.append((path, m.group(0).strip()))
        elif name == "NotebookEdit":
            for m in MARKER.finditer(inp.get("new_source", "") or ""):
                added.append((path, m.group(0).strip()))
    return added, removed


def _diff_lists(
    added: list[tuple[str, str]], removed: list[tuple[str, str]]
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Subtract overlapping (file, line) tuples — they represent marker
    moves within an Edit, not net additions/removals."""
    rem_counter: dict[tuple[str, str], int] = {}
    for k in removed:
        rem_counter[k] = rem_counter.get(k, 0) + 1
    net_added: list[tuple[str, str]] = []
    for k in added:
        if rem_counter.get(k):
            rem_counter[k] -= 1
        else:
            net_added.append(k)
    net_removed = []
    for k, c in rem_counter.items():
        for _ in range(c):
            net_removed.append(k)
    return net_added, net_removed


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    transcript_path = payload.get("transcript_path", "")
    added, removed = _walk_transcript(transcript_path)
    net_added, net_removed = _diff_lists(added, removed)
    if not net_added and not net_removed:
        return 0

    lines = ["[stale-todo-tracker] session summary:"]
    if net_added:
        lines.append(f"  Added {len(net_added)} TODO/FIXME/XXX/HACK marker(s):")
        for path, marker in net_added[:10]:
            rel = os.path.relpath(path) if os.path.isabs(path) else path
            lines.append(f"    + {rel}: {marker[:100]}")
        if len(net_added) > 10:
            lines.append(f"    (... and {len(net_added) - 10} more)")
    if net_removed:
        lines.append(f"  Removed {len(net_removed)} marker(s):")
        for path, marker in net_removed[:5]:
            rel = os.path.relpath(path) if os.path.isabs(path) else path
            lines.append(f"    - {rel}: {marker[:100]}")
    print("\n".join(lines), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
