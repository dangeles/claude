#!/usr/bin/env python3
"""rm-rf-sanity-guard.py — PreToolUse hook (Bash).

Refuses destructive `rm`/`find -delete` invocations targeting:
  - $HOME or paths matching $HOME/repos/* exactly (not within)
  - filesystem root or `/*`
  - empty / suspicious tokens (`""`, `*`, `.`, `./`)
  - paths outside the current git repo's root (when in a repo)

This is the "typo guard" — most legitimate destruction is targeted
inside the current working tree; bulk-destruction outside is almost
always a typo.

Exit codes:
  0  — allow
  2  — block

Bypass: SKIP_RM_GUARD in tool_input.description, or
        SKIP_RM_GUARD=1 in env.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys

# Match destructive operations. Conservative: only `rm -r*`, `rm -f*r*`,
# and `find ... -delete` qualify. Plain `rm file` is fine.
RM_RECURSIVE = re.compile(r"\brm\s+(?:-[a-zA-Z]*[rRf][a-zA-Z]*\s+)+")
FIND_DELETE = re.compile(r"\bfind\b.*-delete\b")
COMMAND_SEPARATORS = re.compile(r"\n|&&|\|\||;")


def _statements(command: str) -> list[str]:
    cleaned = re.sub(r"\\\s*\n", " ", command)
    return [s.strip() for s in COMMAND_SEPARATORS.split(cleaned) if s.strip()]


def _extract_rm_paths(stmt: str) -> list[str]:
    """Pull non-flag tokens after `rm`."""
    m = re.match(r"\s*rm\s+(.+)$", stmt)
    if not m:
        return []
    try:
        tokens = shlex.split(m.group(1))
    except ValueError:
        tokens = m.group(1).split()
    return [t for t in tokens if t and not t.startswith("-")]


def _extract_find_paths(stmt: str) -> list[str]:
    """Find usually starts with one or more search roots before predicates."""
    m = re.match(r"\s*find\s+(.+?)\s+-", stmt)
    if not m:
        return []
    try:
        tokens = shlex.split(m.group(1))
    except ValueError:
        tokens = m.group(1).split()
    return tokens


def _is_dangerous(path: str) -> tuple[bool, str]:
    """Return (is_dangerous, reason)."""
    home = os.path.expanduser("~")
    repos_root = os.path.join(home, "repos")

    if not path:
        return True, "empty path"
    raw = path
    expanded = os.path.expanduser(path).rstrip("/")

    # Suspicious literal tokens
    if raw.strip() in ("", "*", ".", "./", "/", "$HOME"):
        return True, f"suspicious bare token: {raw!r}"
    # Absolute roots
    if expanded in ("", "/"):
        return True, "filesystem root"
    # Home dir exactly
    if expanded == home.rstrip("/"):
        return True, "$HOME exactly"
    # repos root exactly (e.g. `rm -rf ~/repos/`)
    if expanded == repos_root:
        return True, "~/repos exactly"
    # `~/repos/<single-dir>` — wiping a whole project repo
    if (
        expanded.startswith(repos_root + "/")
        and expanded[len(repos_root) + 1 :].rstrip("/").count("/") == 0
    ):
        return True, f"a top-level repo directory ({expanded})"

    return False, ""


def _outside_git_repo(path: str) -> tuple[bool, str]:
    """If we're inside a git repo, refuse rm on paths outside that repo."""
    try:
        repo_root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False, ""  # not in a repo; skip this check
    expanded = os.path.abspath(os.path.expanduser(path))
    if not (expanded == repo_root or expanded.startswith(repo_root + "/")):
        return True, f"outside repo root {repo_root}"
    return False, ""


def main() -> int:
    if os.environ.get("SKIP_RM_GUARD") == "1":
        return 0
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    if payload.get("tool_name") != "Bash":
        return 0
    tool_input = payload.get("tool_input", {}) or {}
    command = tool_input.get("command", "") or ""
    description = tool_input.get("description", "") or ""
    if "SKIP_RM_GUARD" in description:
        return 0

    flagged: list[tuple[str, str, str]] = []  # (path, command_snippet, reason)
    for stmt in _statements(command):
        if RM_RECURSIVE.search(stmt):
            paths = _extract_rm_paths(stmt)
        elif FIND_DELETE.search(stmt):
            paths = _extract_find_paths(stmt)
        else:
            continue
        for p in paths:
            danger, reason = _is_dangerous(p)
            if danger:
                flagged.append((p, stmt[:160], reason))
                continue
            outside, reason = _outside_git_repo(p)
            if outside:
                flagged.append((p, stmt[:160], reason))

    if not flagged:
        return 0

    lines = ["[rm-rf-sanity-guard] Refusing destructive command — dangerous target(s):"]
    for path, snippet, reason in flagged:
        lines.append(f"  - path: {path!r}  reason: {reason}")
        lines.append(f"    in statement: {snippet}")
    lines.append("")
    lines.append("  If this is intentional: include SKIP_RM_GUARD in the description,")
    lines.append("  or export SKIP_RM_GUARD=1.")
    print("\n".join(lines), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
