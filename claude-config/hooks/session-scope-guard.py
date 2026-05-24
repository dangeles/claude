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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _hook_lib import iter_tool_uses  # noqa: E402  # pyright: ignore[reportMissingImports]

BULK_STAGE = re.compile(r"\bgit\s+add\s+(?:-A\b|--all\b|\.(?:\s|$|;|&|\|)|-u\b|--update\b)")
IS_GIT_ADD = re.compile(r"\bgit\s+add\b")
IS_GIT_COMMIT = re.compile(r"\bgit\s+commit\b")

# Command separators for splitting a multi-statement bash command into
# individual statements. Pipe (`|`) is intentionally excluded: a pipeline
# is still one logical statement for our purposes.
COMMAND_SEPARATORS = re.compile(r"\n|&&|\|\||;")

# Best-effort patterns for path-producing Bash sub-invocations.
# Each pattern's group(1) is the path or paths-segment; multi-arg
# commands (rm, mv, cp targets) are post-processed via shlex split.
# Order matters: rm/mv/cp first, then sed -i, then redirects/tee.
BASH_MUTATIONS_SINGLE = [
    # rm: capture everything after flags up to a command boundary;
    # caller splits via shlex to handle `rm -rf a b c`.
    (re.compile(r"\brm\s+(?:-[a-zA-Z]+\s+)*([^\s].*?)(?:&&|;|\|\||\||$)", re.DOTALL), "all"),
    # mv/cp <flags> <src...> <dst>: the destination is the last token.
    # We capture the args segment and pull the last non-flag token.
    (re.compile(r"\bmv\s+(?:-[a-zA-Z]+\s+)*([^\s].*?)(?:&&|;|\|\||\||$)", re.DOTALL), "last"),
    (re.compile(r"\bcp\s+(?:-[a-zA-Z]+\s+)*([^\s].*?)(?:&&|;|\|\||\||$)", re.DOTALL), "last"),
    (re.compile(r"\btee\s+(?:-[a-zA-Z]+\s+)*([^\s].*?)(?:&&|;|\|\||\||$)", re.DOTALL), "non-flag"),
    # sed -i ... <file>: the last token after the sed expression.
    (re.compile(r"\bsed\s+-i(?:\s*'[^']*')?\s+\S+\s+(\S+)", re.DOTALL), "first"),
    # Shell redirects target a single path.
    (re.compile(r"(?:^|[\s;&|])>\s*(\S+)"), "first"),
    (re.compile(r"(?:^|[\s;&|])>>\s*(\S+)"), "first"),
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
    for name, inp in iter_tool_uses(transcript_path):
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
            for p in _extract_bash_mutation_paths(cmd):
                files.add(canonicalize(p, repo_root))
    return files


def _extract_bash_mutation_paths(command: str) -> list[str]:
    """Return every path the command appears to mutate (best-effort)."""
    command = _strip_line_continuations(command)
    out: list[str] = []
    for pat, mode in BASH_MUTATIONS_SINGLE:
        for m in pat.finditer(command):
            segment = m.group(1).strip()
            tokens = _shlex_safe_split(segment)
            args = [t for t in tokens if t and not t.startswith("-")]
            if not args:
                continue
            if mode == "all":
                out.extend(args)
            elif mode == "last":
                out.append(args[-1])
            elif mode == "non-flag":
                out.extend(args)
            elif mode == "first":
                out.append(args[0])
    return out


def _strip_line_continuations(command: str) -> str:
    """Collapse backslash-newline continuations into single spaces."""
    return re.sub(r"\\\s*\n", " ", command)


def _shlex_safe_split(s: str) -> list[str]:
    try:
        return shlex.split(s)
    except ValueError:
        return s.split()


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


def split_commands(command: str) -> list[str]:
    """Split a multi-statement bash command into individual statements.

    Line continuations are collapsed first; then split on \\n, &&, ||, ;.
    Pipe (`|`) is preserved because a pipeline is one logical statement.
    """
    cleaned = _strip_line_continuations(command)
    return [s.strip() for s in COMMAND_SEPARATORS.split(cleaned) if s.strip()]


def extract_git_add_paths(command: str, repo_root: str | None = None) -> list[str]:
    """Parse `git add <paths>` from each statement in the command;
    returns every listed path across all statements.

    Multi-statement input like `git status; git add foo; git status`
    correctly picks only `foo` (not `status`) by splitting on command
    separators before regex matching.

    Filters out tokens that don't look like real paths — this guards
    against multi-line commands where heredoc'd commit-message text
    happens to start a "line" with `git add capture regex ...` and
    yields English words as fake paths. Tokens are considered real
    paths only if they contain a path indicator (`/` or `.`) OR exist
    on disk OR are tracked by git.
    """
    out: list[str] = []
    for stmt in split_commands(command):
        m = re.match(r"\s*git\s+add\s+(.+)$", stmt)
        if not m:
            continue
        args_str = m.group(1).split("|", 1)[0]
        tokens = _shlex_safe_split(args_str)
        for t in tokens:
            if t and not t.startswith("-") and t != "\\" and looks_like_path(t, repo_root):
                out.append(t)
    return out


def looks_like_path(token: str, repo_root: str | None) -> bool:
    """True if token plausibly identifies something git could stage.

    Accepts tokens that contain a path separator or extension dot, or
    that resolve to an existing filesystem entry, or that git already
    tracks (covers tracked-but-deleted files). Rejects bare English
    words like 'capture' or 'regex' that crept in from heredoc parsing.
    """
    if not token:
        return False
    if "/" in token or "\\" in token:
        return True
    # An interior dot (foo.txt) is a strong path indicator. Trailing
    # punctuation (e.g. 'capture,') is not.
    if "." in token and not token.endswith((".", ",", ";", ":", "!")):
        return True
    # Filesystem existence (cheap-ish check)
    if os.path.exists(token):
        return True
    if repo_root and os.path.exists(os.path.join(repo_root, token)):
        return True
    # Tracked-but-deleted: git knows about it even though it's gone.
    try:
        subprocess.check_call(
            ["git", "ls-files", "--error-unmatch", "--", token],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def has_git_add(command: str) -> bool:
    return any(re.match(r"\s*git\s+add\b", s) for s in split_commands(command))


def has_git_commit(command: str) -> bool:
    return any(re.match(r"\s*git\s+commit\b", s) for s in split_commands(command))


def has_bulk_stage(command: str) -> bool:
    """True if any statement IS a `git add` invocation that uses a
    bulk-stage flag. Statements that merely mention `git add -A` in
    their text (e.g. inside a `git commit -m "..."` body or a heredoc)
    do not count — only statements where the bulk-stage applies to a
    real `git add` call.
    """
    for s in split_commands(command):
        if not re.match(r"\s*git\s+add\b", s):
            continue
        if BULK_STAGE.search(s):
            return True
    return False


def expand_path_for_session_check(
    path: str, repo_root: str | None
) -> list[str]:
    """If `path` is a directory or directory-like, expand to the actual
    files git would stage under it. Otherwise return [path].

    Uses `git status --porcelain -- <path>` to get the real list of
    dirty files git would pick up. This is what `git add <dir>` stages.
    """
    if not path:
        return []
    # Quick check: treat trailing slash as directory; also probe the
    # filesystem for explicit directories.
    is_dir_arg = path.endswith("/") or _looks_like_dir(path, repo_root)
    if not is_dir_arg:
        return [path]
    try:
        # Use rstrip("\n") (not strip()) so we don't drop the leading
        # space that porcelain emits for unstaged-modified files (" M path").
        out = subprocess.check_output(
            ["git", "status", "--porcelain", "--", path.rstrip("/")],
            text=True,
            stderr=subprocess.DEVNULL,
        ).rstrip("\n")
    except (subprocess.CalledProcessError, FileNotFoundError):
        return [path]
    expanded: list[str] = []
    for line in out.splitlines():
        # Porcelain format: 2 status chars + 1 space + path; rename:
        # "XY old -> new". Status chars may include a literal space
        # (e.g. " M" for unstaged-modified), so the prefix length is
        # fixed at 3 regardless of which char is the space.
        if len(line) < 4:
            continue
        rest = line[3:]
        if " -> " in rest:
            rest = rest.split(" -> ", 1)[1]
        expanded.append(rest.strip().strip('"'))
    return expanded if expanded else [path]


def _looks_like_dir(path: str, repo_root: str | None) -> bool:
    candidates = [path]
    if repo_root and not path.startswith("/"):
        candidates.append(os.path.join(repo_root, path))
    return any(os.path.isdir(c) for c in candidates)


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

    # Only act on git add / git commit. Use statement-aware detection
    # so a `git status; git add foo` command isn't mis-parsed.
    if not (has_git_add(command) or has_git_commit(command)):
        return 0

    # 1. Refuse bulk-stage patterns regardless of session-set status
    if has_bulk_stage(command):
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

    # 2. For specific-path `git add`, check each path. Directory args
    #    are expanded to their constituent files via `git status -- <dir>`
    #    so `git add some/dir/` checks each file git would actually stage.
    if has_git_add(command):
        raw_paths = extract_git_add_paths(command, repo_root)
        expanded: list[str] = []
        for p in raw_paths:
            for f in expand_path_for_session_check(p, repo_root):
                expanded.append(canonicalize(f, repo_root))
        out_of_session = [p for p in expanded if p and p not in session_files]
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
    if has_git_commit(command):
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
