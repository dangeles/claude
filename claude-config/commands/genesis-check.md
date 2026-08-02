---
name: genesis-check
description: Audit a Genesis surface against the design system, reporting by rule ID
argument-hint: "[path — defaults to the current repository]"
allowed-tools: Bash(python3:*), Bash(git rev-parse:*), Bash(test:*), Read
---

Audit a Genesis surface against `~/repos/aesthetic_and_web_guidelines`, reporting findings
by rule ID.

Target: `$ARGUMENTS` (empty means the current repository).

Run this, then interpret the output:

```bash
GUIDE=~/repos/aesthetic_and_web_guidelines
TOOL="$GUIDE/tools/check_surface.py"

test -f "$TOOL" || { echo "genesis-check: $GUIDE not found — clone it, or skip this check."; exit 0; }

# check_surface.py globs <target>/**/*.css. Handed a single FILE it finds zero
# stylesheets and reports a clean bill of health, which is indistinguishable
# from a real pass. So always resolve to a repository root.
TARGET="${ARGUMENTS:-.}"
if [ -f "$TARGET" ]; then TARGET=$(dirname "$TARGET"); fi
TARGET=$(cd "$TARGET" 2>/dev/null && git rev-parse --show-toplevel 2>/dev/null || echo "$TARGET")

python3 "$TOOL" --json "$TARGET"
echo "---exit=$?"
```

## Reading the result

**Exit status is not success or failure.** `check_surface.py` exits 1 when it finds
something at or above the severity gate. That is a result. Judge by whether stdout parsed as
JSON: valid JSON means the tool worked, whatever the exit code. Only treat it as a tool error
if stdout is unparseable.

Pointed at `~/repos/aesthetic_and_web_guidelines` itself, the tool refuses — that repo checks
itself inward with `check_tokens.py`. Run that instead.

## Reporting back

Summarise by rule, worst first. For each finding give the rule ID, the file and line, and the
`guideline` anchor — every finding carries one, so the paragraph arguing for the rule is one
`Read` away.

Do not restate a rule from memory, and do not hand-compute a contrast ratio; get ratios from
`palette_audit.py`. If the findings need interpretation, invoke the `genesis-design` skill.

**Do not fix anything unless asked.** This command observes. Adoption is the surface owner's
decision.
