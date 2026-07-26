# Quality Gates for Skill Editing

The three gates used in the skill-editor workflow, with the commands that check them.

## Contents

- Gate 1: Specification Approval (Phase 1)
- Gate 2: Plan Approval (Phase 3)
- Gate 3: Execution Verification (Phase 4)
- Escalation

---

## Gate 1: Specification Approval (Phase 1)

**Purpose**: the request is specific and well-scoped before work begins.

**Owner**: whoever produced the spec — `skill-editor-request-refiner` if it was delegated,
otherwise the orchestrator.

A usable specification is specific rather than vague ("Add parallel execution to Phase 2
agents", not "Make the skill better"), states what done looks like in measurable terms,
draws the IN/OUT scope boundary, and does not conflict with existing skills or naming
conventions.

There is no automated check here. The gate is the user's approval.

**On failure**: ask the remaining clarifying questions, or present the options if the
request has several valid readings.

**Example**:

```markdown
Original Request: "Edit the researcher skill"

Refined Specification:
- Objective: Add parallel web search to researcher skill
- Success Criteria:
  - researcher can invoke 3 simultaneous web searches
  - Results aggregated before synthesis
  - No regression in existing functionality
- Scope:
  - IN: Modify researcher/SKILL.md workflow
  - OUT: No changes to external-researcher agent
- User Approval: confirmed
```

---

## Gate 2: Plan Approval (Phase 3)

**Purpose**: the implementation plan is concrete enough to execute and safe to execute.

**Owner**: `skill-editor-adversarial-reviewer` when review was run, plus the user.

The plan names exact file paths for every edit, creation, and deletion; specifies the git
workflow (which files to stage, the commit message format); identifies integration points
(other skills or agents, sync-config.py interaction, planning journal update); and carries
the reviewer's verdict if a review ran.

**Checks**:

```bash
# Plan has concrete file operations
grep -E "^(Edit|Create|Delete):" implementation-plan.md | wc -l
# Should be > 0

# All target paths are under claude-config/
grep -E "^(Edit|Create):" implementation-plan.md | grep -v "claude-config/"
# Should return empty

# Git workflow specified
grep -i "git add\|git commit" implementation-plan.md
```

When the target is an orchestrator, the orchestrator pattern coverage counts are advisory:
record them and note which absent patterns matter for this architecture. Missing patterns
do not fail this gate.

**On failure**: revise the plan and, if a review ran, re-review.

---

## Gate 3: Execution Verification (Phase 4)

**Purpose**: the change is valid, synced, and has not broken anything. This is the gate that
protects the live `~/.claude/` config, so it is entirely objective — validators and exit
codes, not a prose read-through.

**Owner**: `skill-editor-executor`, or the orchestrator for small edits.

### Pre-sync validation

```bash
# Validate YAML frontmatter
for skill in claude-config/skills/*/SKILL.md; do
  python3 << EOF
import yaml
with open('$skill', 'r') as f:
    content = f.read()
parts = content.split('---', 2)
yaml.safe_load(parts[1])
print("ok $skill")
EOF
done

# Validate JSON agent files
for agent in claude-config/agents/*.json; do
  python3 -m json.tool "$agent" > /dev/null && echo "ok $agent"
done

# Check for duplicate skill names
find claude-config/skills -name "SKILL.md" -exec grep "^name:" {} + | sort | uniq -d
# Should be empty

# Dry-run sync
./sync-config.py push --dry-run
```

Do not sync until these pass. A malformed skill breaks the live config, and the push-time
gates will refuse the push anyway.

### Post-sync verification

```bash
# Sync left no divergence
./sync-config.py status

# Skill file landed and parses
ls "$HOME/.claude/skills/$SKILL_NAME/SKILL.md"

# Planning journal entry exists
ls -l planning/$(hostname -s)/*.md | tail -1

# Working tree state
git status
```

Also confirm the skill invokes without error and that existing skills still work (smoke
test).

**On failure**: roll back with `git reset --hard HEAD`, re-sync via `./sync-config.py push`,
document the failure in the planning journal, and return to the phase that needs redoing.

---

## Escalation

From CONFIG_MANAGEMENT.md decision thresholds:

- **Ask the user (AskUserQuestion)** before adding a quality gate, changing a gate's
  pass/fail criteria, skipping a gate, re-running a phase after a gate failure, or modifying
  a gate's checks.
- **Decide yourself, then report**: fixing a validation syntax error, updating the planning
  entry, re-running a validation command.
