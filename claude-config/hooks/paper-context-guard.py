#!/usr/bin/env python3
"""paper-context-guard.py — PreToolUse hook (Read, WebFetch).

Non-blocking. Fires when the main thread is about to pull a full-text
scientific paper into its own context: a local PDF, or a URL on a
publisher / preprint / repository host.

A full text runs 8-15k tokens, most of it methods boilerplate and
references that never get referred to again. Three papers read inline
cost more context than the conversation that follows. The paper-reader
skill's rule is one subagent per paper, notes back to the main thread.

This only nudges — reading a paper inline is sometimes right (short
comment, working through a passage with the user, returning to a
section already known to matter). Exit code 0 always.

Subagents are exempt: this hook is about main-thread context, and the
subagent reading the paper IS the recommended pattern.

Suppress with `noqa-paper-context` in the tool call's description.
"""
from __future__ import annotations

import json
import os
import re
import sys

PDF_SUFFIX = re.compile(r"\.pdf(\?|#|$)", re.IGNORECASE)

# Publisher, preprint, and full-text repository hosts. Matched against the
# URL host only, so a path segment like /pmc/ elsewhere won't trigger.
PAPER_HOSTS = (
    "biorxiv.org",
    "medrxiv.org",
    "arxiv.org",
    "pmc.ncbi.nlm.nih.gov",
    "ncbi.nlm.nih.gov",
    "nature.com",
    "science.org",
    "sciencedirect.com",
    "cell.com",
    "pnas.org",
    "elifesciences.org",
    "journals.plos.org",
    "embopress.org",
    "onlinelibrary.wiley.com",
    "link.springer.com",
    "academic.oup.com",
    "tandfonline.com",
    "frontiersin.org",
    "mdpi.com",
    "jbc.org",
    "ahajournals.org",
    "nejm.org",
    "thelancet.com",
    "bmj.com",
    "doi.org",
)

MESSAGE = (
    "\n[paper-context-guard] About to read a full-text paper into the main\n"
    "  context. A full text is ~8-15k tokens, mostly methods and references.\n"
    "  The paper-reader skill's default is one subagent per paper, returning\n"
    "  notes — dispatch several concurrently when there's a stack of them.\n"
    "  Read inline anyway if it's short, if you're working through a passage\n"
    "  with the user, or if you know which section you need.\n"
    "  Suppress: add `noqa-paper-context` to the tool call's description.\n"
)


def _host_of(url: str) -> str:
    """Lowercased host of a URL, empty string if unparseable."""
    match = re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*://([^/?#]+)", url.strip())
    if not match:
        return ""
    host = match.group(1).lower()
    # Strip userinfo and port.
    host = host.rsplit("@", 1)[-1].rsplit(":", 1)[0]
    return host


def _is_paper_host(host: str) -> bool:
    return any(host == h or host.endswith("." + h) for h in PAPER_HOSTS)


def _is_subagent(payload: dict) -> bool:
    """True when this hook fires inside a subagent rather than the main thread.

    Claude Code does not expose a stable field for this, so check the
    documented-ish signals and fall back to False (warn) — a spurious
    warning in a subagent is cheap, a missed one in the main thread is
    the thing we care about.
    """
    for key in ("is_subagent", "isSubagent", "subagent"):
        if payload.get(key):
            return True
    if payload.get("agent_type") not in (None, "", "main"):
        return True
    return bool(os.environ.get("CLAUDE_AGENT_TYPE"))


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0

    tool = payload.get("tool_name")
    if tool not in ("Read", "WebFetch"):
        return 0

    tool_input = payload.get("tool_input", {}) or {}
    if "noqa-paper-context" in (tool_input.get("description", "") or ""):
        return 0
    if _is_subagent(payload):
        return 0

    if tool == "Read":
        target = tool_input.get("file_path", "") or ""
        hit = bool(PDF_SUFFIX.search(target))
    else:
        target = tool_input.get("url", "") or ""
        hit = bool(PDF_SUFFIX.search(target)) or _is_paper_host(_host_of(target))

    if not hit:
        return 0

    print(MESSAGE, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
