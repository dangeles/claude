#!/usr/bin/env python3
"""genesis-token-guard.py — PreToolUse hook (Write, Edit).

Non-blocking. Exit code 0 always, on every path.

Fires when a stylesheet belonging to a Genesis surface is about to be
written, and points at the design system that governs it. The surfaces
are ~/repos/white_papers, ~/repos/papers and ~/repos/kol; the system is
~/repos/aesthetic_and_web_guidelines.

WHY THIS IS A NUDGE AND NOT A CHECKER
-------------------------------------
An earlier design inspected the CSS being written and cited specific
rules — hex literals, font weights, focus rings. It was measured against
the real files and reimplemented those rules wrongly: the weight check
alone produced six false positives against six true findings in
papers/web/static/library.css, because SN Leif genuinely is OS/2 class
600 and the canonical checker exempts it.

The rules live in check_surface.py in the guidelines repository. A copy
here would be a second source of truth in a third repository, with no
shared fixtures pinning the two together — which is the exact failure
the design system exists to prevent. So this hook does not judge
content. It routes to the checker, which is canonical and takes about
125 ms on a whole repository.

The one content-free exception is tokens.css. That file is meant to be
copied verbatim from the guidelines repo (ADOPT-01), so its filename
alone justifies a more specific message. A filename cannot drift.

Subagents are exempt. Suppress with `noqa-genesis` in the tool call's
description.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time

# Guarded: if _hook_lib is ever missing mid-push, this hook must still
# exit 0 rather than traceback on every Write on the machine.
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _hook_lib import report  # pyright: ignore[reportMissingImports]
except Exception:  # pragma: no cover - defensive
    def report(message: str, *, block: bool = True) -> int:  # type: ignore[misc]
        print(message, file=sys.stderr)
        return 2 if block else 0


GUIDELINES = "~/repos/aesthetic_and_web_guidelines"

# Surface repositories governed by the design system.
SURFACE_ROOTS = (
    "~/repos/white_papers",
    "~/repos/papers",
    "~/repos/kol",
)

STYLESHEET_SUFFIXES = (".css", ".scss")

# Mirrors SKIP_DIRS in the guidelines repo's check_surface.py, so the two
# tools agree on what is ours. The cautionary case is Playwright's bundled
# trace viewer under kol, which an unfiltered checker reported thousands of
# findings against.
VENDOR_PARTS = frozenset({
    ".venv", "venv", "site-packages", "node_modules", "build", "dist",
    "vendor", "third_party", ".git", "__pycache__",
})

MINIFIED = re.compile(r"\.(?:min|bundle)\.(?:css|scss)$", re.IGNORECASE)

MARKER_TTL_DAYS = 30


def _real(path: str) -> str:
    """Expand ~ and resolve symlinks. Returns '' on anything unusable."""
    try:
        return os.path.realpath(os.path.expanduser(path))
    except Exception:
        return ""


def _is_subagent(payload: dict) -> bool:
    """True when this fires inside a subagent rather than the main thread.

    Claude Code exposes no stable field, so check the documented-ish
    signals and fall back to False. Mirrors paper-context-guard.py.
    """
    for key in ("is_subagent", "isSubagent", "subagent"):
        if payload.get(key):
            return True
    if payload.get("agent_type") not in (None, "", "main"):
        return True
    return bool(os.environ.get("CLAUDE_AGENT_TYPE"))


def _governed(path: str) -> "tuple[str, str] | None":
    """Is this a Genesis stylesheet? Returns (repo_root, repo_slug) or None.

    A stylesheet qualifies if it sits under a known surface repository, or
    if some component of its path is exactly 'genesis'. The path component
    test is deliberately not a substring test: 'genesis' as a substring
    also matches ~/repos/organogenesis_notion, which has nothing to do
    with the design system and 185 files that are not stylesheets.
    """
    if not path.lower().endswith(STYLESHEET_SUFFIXES):
        return None
    if MINIFIED.search(path):
        return None

    parts = path.split(os.sep)
    if VENDOR_PARTS.intersection(parts):
        return None

    for root in SURFACE_ROOTS:
        real_root = _real(root)
        if real_root and (path == real_root or path.startswith(real_root + os.sep)):
            return real_root, os.path.basename(real_root)

    for part in parts:
        if part.lower() == "genesis":
            return os.path.dirname(path), "genesis"

    # tokens.css anywhere is in scope — it is the file the adoption
    # contract is about, and a fork of it is the failure worth catching.
    if os.path.basename(path).lower() == "tokens.css":
        return os.path.dirname(path), "tokens"

    return None


def _marker_dir():
    from pathlib import Path
    return Path(os.path.expanduser("~/.claude/session-state/genesis-token-guard"))


def _reap(marker_dir) -> None:
    """Drop markers older than MARKER_TTL_DAYS. Best-effort, create-path only."""
    cutoff = time.time() - MARKER_TTL_DAYS * 86400
    try:
        for entry in marker_dir.iterdir():
            try:
                if entry.stat().st_mtime < cutoff:
                    entry.unlink()
            except OSError:
                continue
    except OSError:
        return


def _claim_pointer(session_id: str, repo_slug: str) -> bool:
    """True if the generic pointer has not yet fired for this session+repo.

    State lives under ~/.claude/session-state/, which is outside every
    sync_rules.always path in sync.config.yaml and therefore invisible to
    sync-config.py's pull, push, --delete and status.

    With no session_id there is nothing to key on, so emit and write no
    marker: a repeated pointer is better than one that silently fires
    once ever and is never seen again.
    """
    if not session_id:
        return True

    safe_id = re.sub(r"[^A-Za-z0-9_\-]", "_", session_id)[:64]
    safe_repo = re.sub(r"[^A-Za-z0-9_\-]", "_", repo_slug)[:64] or "repo"

    try:
        marker_dir = _marker_dir()
        marker_dir.mkdir(parents=True, exist_ok=True)
        _reap(marker_dir)
        marker = marker_dir / "{}__{}.marker".format(safe_id, safe_repo)
        # Atomic create-or-fail: concurrent writes cannot both claim it.
        with open(str(marker), "x"):
            pass
        return True
    except FileExistsError:
        return False
    except OSError:
        # Unwritable ~/.claude. Emit rather than go silent.
        return True


def _adopt_message(path: str) -> str:
    return (
        "ADOPT-01 — `tokens.css` is copied verbatim from the design system, not forked.\n"
        "  If this surface needs a token the canonical file lacks, add it there first,\n"
        "  run `python3 {g}/tools/check_tokens.py`, then copy again.\n"
        "  Compare:  diff {g}/tokens/tokens.css {p}\n"
        "  Invoke the `genesis-design` skill before continuing."
    ).format(g=GUIDELINES, p=path)


def _pointer_message(repo_root: str) -> str:
    return (
        "This stylesheet is a Genesis surface, governed by {g}.\n"
        "  `python3 {g}/tools/check_surface.py --json {r}` reports by rule ID in ~0.1s\n"
        "  (or run /genesis-check). It exits 1 when it finds something — that is a\n"
        "  result, not a failure.\n"
        "  Invoke the `genesis-design` skill before continuing."
    ).format(g=GUIDELINES, r=repo_root)


def main() -> int:
    try:
        try:
            payload = json.loads(sys.stdin.read() or "{}")
        except json.JSONDecodeError:
            return 0
        if not isinstance(payload, dict):
            return 0

        if payload.get("tool_name") not in ("Write", "Edit"):
            return 0

        tool_input = payload.get("tool_input", {}) or {}
        if not isinstance(tool_input, dict):
            return 0

        if "noqa-genesis" in (tool_input.get("description", "") or ""):
            return 0
        if _is_subagent(payload):
            return 0

        path = _real(tool_input.get("file_path", "") or "")
        if not path:
            return 0

        governed = _governed(path)
        if governed is None:
            return 0
        repo_root, repo_slug = governed

        if os.path.basename(path).lower() == "tokens.css":
            return report(_adopt_message(path), block=False)

        if _claim_pointer(payload.get("session_id", "") or "", repo_slug):
            return report(_pointer_message(repo_root), block=False)

        return 0
    except Exception:
        # Fail open, always. A traceback here would land on every Write
        # and Edit on the machine.
        return 0


if __name__ == "__main__":
    sys.exit(main())
