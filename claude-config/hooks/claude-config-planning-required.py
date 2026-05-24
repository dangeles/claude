#!/usr/bin/env python3
"""claude-config-planning-required.py — Stop hook.

Scoped to the claude-config repo. When Claude is about to stop a
session that modified anything under claude-config/skills/,
claude-config/agents/, claude-config/hooks/, or
claude-config/settings.json, check whether a new planning entry was
created in this session under planning/mac/. If not, emit a reminder
to stderr and force continuation (exit 2 in a Stop hook means "don't
stop yet").

The check is best-effort:
  - Identifies the claude-config repo via git rev-parse
  - Walks the conversation transcript for Write/Edit/NotebookEdit tool
    calls that target paths under the gated directories
  - Considers a planning entry "created" if any Write call this session
    targeted `planning/.+\\.md$`

Bypass: SKIP_PLANNING_REQUIRED=1 env var, or include
SKIP_PLANNING_REQUIRED in the user's most recent prompt.

Exit codes:
  0  — let the session stop normally
  2  — block stop; Claude must continue (typically: create the entry)
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _hook_lib import iter_tool_uses, report  # noqa: E402  # pyright: ignore[reportMissingImports]

GATED_PATH_PATTERN = re.compile(
    r"claude-config/(?:skills/|agents/|hooks/|settings\.json|plugins/)"
)
PLANNING_ENTRY_PATTERN = re.compile(r"planning/[^/]+/.+\.md$")
EXPECTED_REPO_BASENAME = "claude"


def _git_repo_basename(cwd: str) -> str:
    try:
        root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=cwd,
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=2,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return ""
    return os.path.basename(root)


def _walk_transcript(transcript_path: str) -> tuple[bool, bool, list[str]]:
    """Return (touched_gated, created_planning_entry, sample_paths)."""
    touched_gated = False
    created_planning = False
    samples: list[str] = []
    for name, inp in iter_tool_uses(transcript_path):
        if name not in ("Write", "Edit", "NotebookEdit"):
            continue
        path = inp.get("file_path") or inp.get("notebook_path") or ""
        if not path:
            continue
        if GATED_PATH_PATTERN.search(path):
            touched_gated = True
            if len(samples) < 5:
                samples.append(path)
        if PLANNING_ENTRY_PATTERN.search(path) and name == "Write":
            created_planning = True
    return touched_gated, created_planning, samples


def main() -> int:
    if os.environ.get("SKIP_PLANNING_REQUIRED") == "1":
        return 0
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    cwd = payload.get("cwd") or os.getcwd()
    if _git_repo_basename(cwd) != EXPECTED_REPO_BASENAME:
        return 0  # only enforce in the claude-config repo

    transcript_path = payload.get("transcript_path", "")
    touched, created, samples = _walk_transcript(transcript_path)
    if not touched:
        return 0
    if created:
        return 0

    msg = [
        "[claude-config-planning-required] This session modified claude-config/ but no",
        "planning entry was created.",
        "",
        "  Per CLAUDE.md (step 2 of the workflow), every change to the synced",
        "  configuration should have a corresponding entry under planning/mac/.",
        "  Sample touched paths:",
    ]
    for s in samples:
        msg.append(f"    - {s}")
    msg.append("")
    msg.append("  Fix: `./sync-config.py plan --title \"...\"`, then populate the entry")
    msg.append("  and add it to the next commit.")
    msg.append("")
    msg.append("  Bypass: include SKIP_PLANNING_REQUIRED in your message, or export")
    msg.append("  SKIP_PLANNING_REQUIRED=1.")
    return report("\n".join(msg), block=True)


if __name__ == "__main__":
    sys.exit(main())
