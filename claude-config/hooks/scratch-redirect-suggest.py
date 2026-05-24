#!/usr/bin/env python3
"""scratch-redirect-suggest.py — PostToolUse hook (Write).

Non-blocking. When a file is written to a repo's root with a scratch-y
filename (`temp_*`, `scratch_*`, `wip_*`, `debug_*`, `Untitled*`),
emit a notice suggesting `~/repos/scratch/` instead.

This encodes the existing convention: the dedicated `scratch/` repo
exists for exactly this purpose; other repos shouldn't accumulate
ephemeral test files at root.

Exit code 0 always — informational only.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

SCRATCHY = re.compile(
    r"^(?:temp|tmp|scratch|wip|debug|untitled|test_quick|throwaway|delete_me|todo)"
    r"[_\-. ]?\w*\.\w+$",
    re.IGNORECASE,
)


def _repo_root(path: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=os.path.dirname(path) or ".",
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=2,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return None


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    if payload.get("tool_name") != "Write":
        return 0
    tool_input = payload.get("tool_input", {}) or {}
    file_path = tool_input.get("file_path", "") or ""
    if not file_path:
        return 0

    name = os.path.basename(file_path)
    if not SCRATCHY.match(name):
        return 0

    # Only suggest if writing to a repo's root (or one directory deep)
    abs_path = os.path.abspath(file_path)
    repo_root = _repo_root(abs_path)
    if not repo_root:
        return 0
    try:
        rel = Path(abs_path).relative_to(repo_root).as_posix()
    except ValueError:
        return 0
    # Top-level only — if it's already buried deep, don't nag
    if "/" in rel.rstrip("/"):
        return 0
    # If we're already in ~/repos/scratch, no nag
    if repo_root.rstrip("/") == os.path.expanduser("~/repos/scratch"):
        return 0

    print(
        f"\n[scratch-redirect-suggest] '{name}' looks scratch-y but was written to "
        f"the root of {os.path.basename(repo_root)}/.\n"
        f"  Consider ~/repos/scratch/ for ephemeral files (matches your existing convention).\n",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
