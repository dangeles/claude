#!/usr/bin/env python3
"""secret-leak-guard.py — PreToolUse hook (Write/Edit).

Refuses to write content that contains known secret-shaped patterns.
Catches the common cases (cloud provider keys, GitHub tokens, JWT,
PEM private keys, hard-coded passwords with high-entropy values).

Exit codes:
  0  — allow the tool call
  2  — block; stderr shown to Claude

Bypass: SKIP_SECRET_LEAK_GUARD anywhere in the tool_input.description.
"""
from __future__ import annotations

import json
import os
import re
import sys

# Patterns ordered roughly by specificity. Each entry: (name, regex).
SECRET_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("Anthropic API key", re.compile(r"\bsk-ant-(?:api03|admin01)-[A-Za-z0-9_\-]{60,}\b")),
    ("OpenAI API key (sk-)", re.compile(r"\bsk-(?:proj-)?[A-Za-z0-9_\-]{32,}\b")),
    ("GitHub personal access token", re.compile(r"\bghp_[A-Za-z0-9]{36,}\b")),
    ("GitHub fine-grained PAT", re.compile(r"\bgithub_pat_[A-Za-z0-9_]{82,}\b")),
    ("GitHub OAuth token", re.compile(r"\bgho_[A-Za-z0-9]{36,}\b")),
    ("AWS access key id", re.compile(r"\bAKIA[0-9A-Z]{16}\b")),
    ("AWS secret access key", re.compile(r"(?i)aws[_\- ]?secret[_\- ]?(?:access[_\- ]?)?key['\"\s:=]+[A-Za-z0-9/+=]{40}\b")),
    ("Slack bot token", re.compile(r"\bxox[bpars]-[A-Za-z0-9-]{10,}\b")),
    ("Google API key", re.compile(r"\bAIza[0-9A-Za-z_\-]{35}\b")),
    ("Stripe live key", re.compile(r"\b(?:sk|rk)_live_[A-Za-z0-9]{20,}\b")),
    ("JWT", re.compile(r"\beyJ[A-Za-z0-9_\-]{10,}\.[A-Za-z0-9_\-]{10,}\.[A-Za-z0-9_\-]{5,}\b")),
    ("PEM private key", re.compile(r"-----BEGIN (?:RSA |EC |DSA |OPENSSH |ENCRYPTED |PGP )?PRIVATE KEY-----")),
    # Hard-coded password assignments with non-trivial values (heuristic;
    # excludes obvious placeholders like "password", "***", "<...>").
    (
        "hard-coded password/token assignment",
        re.compile(
            r"""(?ix)
            (?:password|passwd|api[_-]?key|secret|access[_-]?token)
            \s*[:=]\s*
            ['"]
            (?![\*<>{}\[\]xX]|null|None|TODO|REPLACE|YOUR_|EXAMPLE|placeholder|password|secret|changeme|test)
            [A-Za-z0-9_\-./+=]{8,}
            ['"]
            """
        ),
    ),
]


def _extract_content(tool_name: str, tool_input: dict) -> str:
    """Pull the content that's about to be written. Different tools
    expose content under different keys."""
    if tool_name == "Write":
        return tool_input.get("content", "") or ""
    if tool_name == "Edit":
        # Check the new value (what's being added). Old string check
        # would catch removing-secrets, which is fine.
        return tool_input.get("new_string", "") or ""
    if tool_name == "NotebookEdit":
        return tool_input.get("new_source", "") or ""
    return ""


def _scan(content: str) -> list[tuple[str, str]]:
    """Return list of (pattern_name, matched_excerpt) findings."""
    findings: list[tuple[str, str]] = []
    for name, pat in SECRET_PATTERNS:
        m = pat.search(content)
        if m:
            excerpt = m.group(0)
            if len(excerpt) > 80:
                excerpt = excerpt[:40] + "..." + excerpt[-20:]
            findings.append((name, excerpt))
    return findings


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0  # malformed → fail open

    tool_name = payload.get("tool_name", "")
    if tool_name not in ("Write", "Edit", "NotebookEdit"):
        return 0

    tool_input = payload.get("tool_input", {}) or {}
    description = tool_input.get("description", "") or ""
    if "SKIP_SECRET_LEAK_GUARD" in description:
        return 0
    if os.environ.get("SKIP_SECRET_LEAK_GUARD") == "1":
        return 0

    content = _extract_content(tool_name, tool_input)
    if not content:
        return 0

    findings = _scan(content)
    if not findings:
        return 0

    lines = ["[secret-leak-guard] Refusing write — secret-shaped content detected:"]
    for name, excerpt in findings:
        lines.append(f"  - {name}: {excerpt}")
    lines.append("")
    lines.append("  If this is a false positive, include SKIP_SECRET_LEAK_GUARD in the")
    lines.append("  tool call's description (or export SKIP_SECRET_LEAK_GUARD=1).")
    lines.append("  Otherwise: move the secret to an env var or a gitignored file,")
    lines.append("  rotate it immediately if it was ever committed elsewhere.")
    print("\n".join(lines), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
