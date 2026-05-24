#!/usr/bin/env python3
"""session-scope-guard.py — PreToolUse hook that prevents Claude from
committing files it didn't modify this session.

Hooks into Bash tool calls. For `git add` / `git commit` commands:
  1. Refuses bulk-stage patterns outright (-A, --all, ., -u, --update).
  2. For specific-path `git add`, refuses if any path isn't in the
     session set derived from this conversation's transcript.
  3. For `git commit`, refuses if any currently-staged path isn't in
     the session set.

Session set is derived from transcript JSONL events: Write/Edit/
NotebookEdit tool inputs, plus best-effort parse of Bash mutations
(mv/cp/rm dst, > redirects, tee, sed -i).

Exit codes (Claude Code PreToolUse contract):
  0  — allow the tool call to proceed
  2  — block the tool call; stderr is shown back to Claude

Bypass (use sparingly, escalate to user first):
  - export SESSION_SCOPE_BYPASS=1
  - include the literal string "SKIP_SESSION_SCOPE_GUARD" in the
    Bash tool call's `description` field

Conservative posture: if the session set cannot be determined (no
transcript, no parseable events), fail OPEN — never block when
uncertain. Bulk-stage patterns are the one exception: those are
refused regardless of session-set status.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

BULK_STAGE = re.compile(r"\bgit\s+add\s+(?:-A\b|--all\b|\.(?:\s|$|;|&|\|)|-u\b|--update\b)")
IS_GIT_ADD = re.compile(r"\bgit\s+add\b")
IS_GIT_COMMIT = re.compile(r"\bgit\s+commit\b")

# Best-effort Bash mutation patterns (path-producing operations).
# Each pattern's group(1) is treated as a modified path.
BASH_MUTATIONS = [
    re.compile(r"\bmv\s+(?:-[a-zA-Z]+\s+)*\S+\s+(\S+)"),
    re.compile(r"\bcp\s+(?:-[a-zA-Z]+\s+)*\S+\s+(\S+)"),
    re.compile(r"\brm\s+(?:-[a-zA-Z]+\s+)*(\S+)"),
    re.compile(r"\btee\s+(?:-[a-zA-Z]+\s+)*(\S+)"),
    re.compile(r"\bsed\s+-i(?:\s*'\S*')?\s+\S+\s+(\S+)"),
    re.compile(r"(?:^|[\s;&|])>\s*(\S+)"),
    re.compile(r"(?:^|[\s;&|])>>\s*(\S+)"),
]


def block(message: str) -> int:
    print(f"\n[session-scope-guard] {message}\n", file=sys.stderr)
    return 2


def canonicalize(path: str, repo_root: str | None) -> str:
    """Normalize to repo-relative form when possible."""
    p = path.strip().strip("'\"")
    if not p:
        return p
    p = os.path.expanduser(p)
    if p.startswith("/") and repo_root:
        try:
            p = str(Path(p).relative_to(repo_root))
        except ValueError:
            pass
    return p


def get_repo_root() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def derive_session_files(transcript_path: str, repo_root: str | None) -> set[str]:
    """Walk transcript JSONL and collect paths Claude has touched."""
    files: set[str] = set()
    if not transcript_path:
        return files
    p = Path(transcript_path)
    if not p.exists():
        return files
    try:
        with p.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                _collect_from_event(event, files, repo_root)
    except OSError:
        pass
    return files


def _collect_from_event(event: dict, files: set[str], repo_root: str | None) -> None:
    """Extract tool_use entries; handles both legacy and v2 transcript shapes."""
    msg = event.get("message", event)
    if not isinstance(msg, dict):
        return
    content = msg.get("content")
    if not isinstance(content, list):
        return
    for item in content:
        if not isinstance(item, dict) or item.get("type") != "tool_use":
            continue
        name = item.get("name", "")
        inp = item.get("input", {}) or {}
        if name in ("Write", "Edit"):
            p = inp.get("file_path") or ""
            if p:
                files.add(canonicalize(p, repo_root))
        elif name == "NotebookEdit":
            p = inp.get("notebook_path") or ""
            if p:
                files.add(canonicalize(p, repo_root))
        elif name == "Bash":
            cmd = inp.get("command", "") or ""
            for pat in BASH_MUTATIONS:
                for m in pat.finditer(cmd):
                    files.add(canonicalize(m.group(1), repo_root))


def get_staged_paths(repo_root: str | None) -> list[str]:
    try:
        out = subprocess.check_output(
            ["git", "diff", "--cached", "--name-only"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    return [canonicalize(p, repo_root) for p in out.splitlines() if p]


def extract_git_add_paths(command: str) -> list[str]:
    """Parse `git add <paths>` segment; returns the listed paths."""
    m = re.search(r"\bgit\s+add\s+(.+)", command)
    if not m:
        return []
    # Stop at command-chaining boundaries
    args_str = re.split(r"&&|;|\|\||\|", m.group(1), maxsplit=1)[0]
    try:
        tokens = shlex.split(args_str)
    except ValueError:
        tokens = args_str.split()
    # Drop flags
    return [t for t in tokens if t and not t.startswith("-")]


def render_paths_list(paths: list[str] | set[str], max_show: int = 12) -> str:
    lst = sorted(paths)
    body = "".join(f"  - {p}\n" for p in lst[:max_show])
    if len(lst) > max_show:
        body += f"  (... and {len(lst) - max_show} more)\n"
    return body


def main() -> int:
    # Bypass via env (set by user, not Claude)
    if os.environ.get("SESSION_SCOPE_BYPASS") == "1":
        return 0

    # Read PreToolUse payload from stdin
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0  # malformed payload — fail open

    if payload.get("tool_name") != "Bash":
        return 0

    tool_input = payload.get("tool_input", {}) or {}
    command = tool_input.get("command", "") or ""
    description = tool_input.get("description", "") or ""

    # Bypass via description annotation
    if "SKIP_SESSION_SCOPE_GUARD" in description:
        return 0

    # Only act on git add / git commit
    if not (IS_GIT_ADD.search(command) or IS_GIT_COMMIT.search(command)):
        return 0

    # 1. Refuse bulk-stage patterns regardless of session-set status
    if BULK_STAGE.search(command):
        return block(
            "Refusing bulk-stage pattern. Stage specific paths only.\n"
            f"  Offending command: {command[:240]}\n"
            "  Fix: list paths explicitly (e.g. `git add path/a path/b`).\n"
            "  Bypass: include SKIP_SESSION_SCOPE_GUARD in the description,\n"
            "  or export SESSION_SCOPE_BYPASS=1 if intentional."
        )

    transcript_path = payload.get("transcript_path", "")
    repo_root = get_repo_root()
    session_files = derive_session_files(transcript_path, repo_root)

    # If we can't determine the session set, fail open (uncertain → allow).
    if not session_files:
        return 0

    # 2. For specific-path `git add`, check each path
    if IS_GIT_ADD.search(command):
        paths = [canonicalize(p, repo_root) for p in extract_git_add_paths(command)]
        out_of_session = [p for p in paths if p and p not in session_files]
        if out_of_session:
            return block(
                "Refusing `git add` — paths not modified this session:\n"
                + render_paths_list(out_of_session)
                + f"\n  Session-modified files ({len(session_files)}):\n"
                + render_paths_list(session_files)
                + "\n  Bypass: SKIP_SESSION_SCOPE_GUARD in description,\n"
                "  or export SESSION_SCOPE_BYPASS=1 if intentional."
            )

    # 3. For `git commit`, check currently-staged paths
    if IS_GIT_COMMIT.search(command):
        staged = get_staged_paths(repo_root)
        out_of_session = [p for p in staged if p not in session_files]
        if out_of_session:
            unstage_cmd = "git reset HEAD " + " ".join(out_of_session[:8])
            if len(out_of_session) > 8:
                unstage_cmd += "  # ...and more"
            return block(
                "Refusing `git commit` — staged paths not modified this session:\n"
                + render_paths_list(out_of_session)
                + f"\n  Unstage them first:\n    {unstage_cmd}\n"
                "\n  Bypass: SKIP_SESSION_SCOPE_GUARD in description,\n"
                "  or export SESSION_SCOPE_BYPASS=1 if intentional."
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
