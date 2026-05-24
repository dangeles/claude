#!/usr/bin/env python3
"""ipynb-output-warn.py — PreToolUse hook (Bash).

Non-blocking. When the command stages a `.ipynb` file via `git add`,
emit a stderr warning that notebook outputs are about to be committed.
The user's notebook-writer skill recommends jupytext markdown sources;
this nudges back toward that workflow without enforcing it (some
notebooks legitimately ship with outputs, e.g. tf_reprogramming_tools/
notebooks/recap_of_paper.ipynb at 17.8 MB).

Exit code 0 always — informational only.

Suppress with `# noqa-ipynb-warn` substring in the tool_input.description.
"""
from __future__ import annotations

import json
import re
import sys

IPYNB_ADD = re.compile(r"\bgit\s+add\b[^&;|\n]*\.ipynb\b")


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read() or "{}")
    except json.JSONDecodeError:
        return 0
    if payload.get("tool_name") != "Bash":
        return 0
    tool_input = payload.get("tool_input", {}) or {}
    command = tool_input.get("command", "") or ""
    description = tool_input.get("description", "") or ""
    if "noqa-ipynb-warn" in description:
        return 0
    if not IPYNB_ADD.search(command):
        return 0
    print(
        "\n[ipynb-output-warn] A `.ipynb` file is about to be staged. Notebooks\n"
        "  serialize their outputs; check before committing whether outputs are\n"
        "  intended or if a jupytext `.md` companion would be better.\n"
        "  Strip outputs: `jupyter nbconvert --clear-output --inplace <file>.ipynb`\n"
        "  Suppress: add `noqa-ipynb-warn` to the tool call's description.\n",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
