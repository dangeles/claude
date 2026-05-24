#!/usr/bin/env python3
"""skill-frontmatter-validator.py — PreToolUse hook (Write/Edit).

For paths matching `claude-config/skills/*/SKILL.md`, validate that
the YAML frontmatter parses and contains the two required fields:
`name:` and `description:`. Catches the same problem that
`sync-config.py push` gates at push-time, but upstream — before the
bad file is written.

Exit codes:
  0  — allow
  2  — block (frontmatter missing or malformed)

Bypass: SKIP_FRONTMATTER_VALIDATION in description.
"""
from __future__ import annotations

import json
import os
import re
import sys

SKILL_PATH = re.compile(r"claude-config/skills/[^/]+/SKILL\.md$")
REQUIRED_FIELDS = ("name", "description")


def _extract_full_content(tool_name: str, tool_input: dict, file_path: str) -> str | None:
    """For Write, content is in tool_input. For Edit, we need the
    full resulting content — reconstruct from old + new strings if
    possible, falling back to reading the on-disk file and applying
    the substitution.
    """
    if tool_name == "Write":
        return tool_input.get("content")
    if tool_name == "Edit":
        old = tool_input.get("old_string", "")
        new = tool_input.get("new_string", "")
        if not os.path.exists(file_path):
            # Nothing to validate against; pass-through.
            return None
        try:
            with open(file_path) as f:
                current = f.read()
        except OSError:
            return None
        if old in current:
            return current.replace(old, new, 1)
        # Edit will fail at apply time; no point validating here.
        return None
    return None


def _validate_frontmatter(content: str) -> tuple[bool, str]:
    """Return (ok, message)."""
    if not content.startswith("---\n"):
        return False, "missing opening `---` line"
    try:
        end = content.index("\n---\n", 4)
    except ValueError:
        return False, "missing closing `---` line"
    fm_text = content[4:end]
    # Minimal manual scan — avoid importing yaml so the hook has zero deps.
    keys: set[str] = set()
    for line in fm_text.splitlines():
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_\-]*)\s*:", line)
        if m:
            keys.add(m.group(1))
    missing = [f for f in REQUIRED_FIELDS if f not in keys]
    if missing:
        return False, f"frontmatter missing required field(s): {', '.join(missing)}"
    # Try a real YAML parse if PyYAML is available; otherwise trust the regex.
    try:
        import yaml  # type: ignore
        meta = yaml.safe_load(fm_text)
        if not isinstance(meta, dict):
            return False, "frontmatter does not parse to a mapping"
        for f in REQUIRED_FIELDS:
            if not meta.get(f):
                return False, f"required field `{f}:` is empty or missing"
    except ImportError:
        pass
    except Exception as e:
        return False, f"YAML parse error: {e}"
    return True, ""


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    tool_name = payload.get("tool_name", "")
    if tool_name not in ("Write", "Edit"):
        return 0
    tool_input = payload.get("tool_input", {}) or {}
    description = tool_input.get("description", "") or ""
    if "SKIP_FRONTMATTER_VALIDATION" in description:
        return 0

    file_path = tool_input.get("file_path", "") or ""
    if not SKILL_PATH.search(file_path):
        return 0

    content = _extract_full_content(tool_name, tool_input, file_path)
    if content is None:
        return 0

    ok, msg = _validate_frontmatter(content)
    if ok:
        return 0
    print(
        f"\n[skill-frontmatter-validator] Refusing write to {file_path}\n"
        f"  Reason: {msg}\n"
        f"  Required fields: {', '.join(REQUIRED_FIELDS)}\n"
        f"  Bypass: include SKIP_FRONTMATTER_VALIDATION in the description.\n",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    sys.exit(main())
