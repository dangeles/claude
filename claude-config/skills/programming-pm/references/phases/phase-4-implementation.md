# Phase 4: Implementation — Full Detail

This reference holds the complete implementation walkthrough for Phase 4 of the
programming-pm workflow. SKILL.md keeps the phase summary, gate criteria, and a
pointer here. Read this file when you are actively executing Phase 4.

## Contents

- Step 0: Architecture Implementability Check
- Step 1: Task Decomposition
- Step 1b: Junior-Developer Evaluation
- Step 2: Wave-Based Parallel Execution (and Wave-Based Specialist Launch)
- Step 3: Completion Tracking
- Step 4: Quality Gate 4a - Specialist Completion Check (and the Phase 4 → 5 handoff validation script)

**Objective**: Implement architecture with specialist agents in parallel.

**Mode-based execution**:
- **SIMPLE**: Sequential execution (one specialist at a time)
- **STANDARD/EXTENDED**: Wave-based parallel execution

**Steps**:

## Step 0: Architecture Implementability Check

**Skip condition**: SIMPLE mode projects with <= 1 component skip Step 0 and proceed directly to Step 1.

Before task decomposition, confirm the architecture handoff is implementable.

**Owner**: programming-pm. Read the handoff and check it yourself when it is small and the
components are straightforward. Dispatch senior-developer via Task tool (using the
"Dispatch to senior-developer (Phase 4 Step 0 Pre-flight)" template) when the architecture
is large, the dependency graph is non-obvious, or an initial read surfaces substantive gaps.

**Context size check**: When dispatching, if the architecture handoff YAML is <= 200 lines
paste it verbatim; if > 200 lines pass the file path and instruct senior-developer to read
it directly.

**What to validate**:
1. All component interfaces are fully specified (inputs have types, outputs have types)
2. Component dependencies are traversable without cycles (check by tracing dependency chains)
3. Technology choices are recognizable and compatible
4. Implementation order in the handoff is internally consistent

**Output** — a validation report at `{SESSION_DIR}/deliverables/phase4-preflight.yaml`:
```yaml
preflight_validation:
  status: "PASS" | "FAIL" | "PASS_WITH_WARNINGS"
  gaps:
    - component: string
      issue: string
      severity: "BLOCKING" | "WARNING"
  recommendations: []
```

**On PASS or PASS_WITH_WARNINGS**: Log any warnings, proceed to Step 1.

**On FAIL**:
1. Escalate back to Phase 3 -- re-invoke systems-architect via Task tool with the gap report
2. systems-architect addresses BLOCKING gaps only
3. Re-run Step 0 (max 2 escalation cycles)
4. If still FAIL after 2 cycles: present to user for override (proceed with warnings, or abort)

If a dispatched Step 0 check fails to return, log "Step 0 pre-flight unavailable" and proceed to Step 1 with a warning.

## Step 1: Task Decomposition

Parse architecture handoff to identify components and assign specialists:

```bash
ARCHITECTURE_FILE="${SESSION_DIR}/handoffs/phase3-architecture-handoff.yaml"

# Extract components
if command -v yq &> /dev/null; then
  COMPONENT_COUNT=$(yq eval '.handoff.components | length' "$ARCHITECTURE_FILE")

  # Initialize task list
  > "$SESSION_DIR/task-assignments.txt"

  # Iterate through components
  for i in $(seq 0 $((COMPONENT_COUNT - 1))); do
    COMPONENT_NAME=$(yq eval ".handoff.components[$i].name" "$ARCHITECTURE_FILE")
    COMPONENT_DESC=$(yq eval ".handoff.components[$i].responsibility" "$ARCHITECTURE_FILE")
    DEPENDENCIES=$(yq eval ".handoff.components[$i].dependencies[]" "$ARCHITECTURE_FILE" 2>/dev/null || echo "")

    SPECIALIST="senior-developer"  # default

    # PRIMARY SIGNAL: Explicit specialist flags from architecture handoff (v1.4+)
    REQUIRES_MATH=$(yq eval ".handoff.components[$i].specialist_flags.requires_mathematician // \"null\"" "$ARCHITECTURE_FILE" 2>/dev/null | tr '[:upper:]' '[:lower:]')
    REQUIRES_STATS=$(yq eval ".handoff.components[$i].specialist_flags.requires_statistician // \"null\"" "$ARCHITECTURE_FILE" 2>/dev/null | tr '[:upper:]' '[:lower:]')
    REQUIRES_NOTEBOOK=$(yq eval ".handoff.components[$i].specialist_flags.requires_notebook_writer // \"null\"" "$ARCHITECTURE_FILE" 2>/dev/null | tr '[:upper:]' '[:lower:]')

    # Normalize YAML boolean variants: true/yes/on -> "true"; false/no/off/null -> not-true
    [[ "$REQUIRES_MATH" =~ ^(true|yes|on)$ ]] && REQUIRES_MATH="true"
    [[ "$REQUIRES_STATS" =~ ^(true|yes|on)$ ]] && REQUIRES_STATS="true"
    [[ "$REQUIRES_NOTEBOOK" =~ ^(true|yes|on)$ ]] && REQUIRES_NOTEBOOK="true"

    if [ "$REQUIRES_MATH" = "true" ]; then
      SPECIALIST="mathematician"
    elif [ "$REQUIRES_STATS" = "true" ]; then
      SPECIALIST="statistician"
    elif [ "$REQUIRES_NOTEBOOK" = "true" ]; then
      SPECIALIST="notebook-writer"
    elif [ "$REQUIRES_MATH" = "null" ] && [ "$REQUIRES_STATS" = "null" ] && [ "$REQUIRES_NOTEBOOK" = "null" ]; then
      # FALLBACK: keyword-based assignment for pre-v1.4 architecture handoffs
      echo "  WARN: No specialist_flags found for '$COMPONENT_NAME'. Using keyword-based fallback (upgrade to v1.4 handoff schema recommended)."
      if echo "$COMPONENT_DESC" | grep -qiE "algorithm|optimization|complexity"; then
        SPECIALIST="mathematician"
      elif echo "$COMPONENT_DESC" | grep -qiE "statistic|hypothesis|regression|bayesian"; then
        SPECIALIST="statistician"
      elif echo "$COMPONENT_DESC" | grep -qiE "notebook|jupyter|ipynb|jupytext|interactive.analysis|parameter.sweep|analysis.report|visualization.notebook|reproducible.analysis|data.exploration"; then
        SPECIALIST="notebook-writer"
      fi
    fi

    # Junior-developer assignment (keyword-based, independent of specialist flags)
    if [ "$SPECIALIST" = "senior-developer" ]; then
      if echo "$COMPONENT_DESC" | grep -qiE "simple|utility|helper|wrapper|validation|config|parser|formatter|converter|serializer|loader|constants|enum|data[._-]class|type[._-]definition|test[._-]fixture|test[._-]helper|boilerplate|scaffold"; then
        SPECIALIST="junior-developer"
      fi
    fi

    # Record task assignment
    TASK_ID="TASK-$(printf "%03d" $((i + 1)))"
    echo "$TASK_ID|$COMPONENT_NAME|$SPECIALIST|$DEPENDENCIES" >> "$SESSION_DIR/task-assignments.txt"
  done
else
  echo "⚠️  yq not found. Manual task decomposition required."
fi
```

**Specialist assignment logic**:
- **Algorithm design** → `mathematician`
- **Statistical methods** → `statistician`
- **Notebook/Jupyter creation** → `notebook-writer`
- **Complex implementation** → `senior-developer`
- **Routine implementation** → `junior-developer` (supervised by senior)

**Task assignment format**:
```yaml
task:
  id: "TASK-001"
  description: string
  assigned_to: skill_name
  dependencies: []
  estimated_duration: "2h"
  acceptance_criteria: []
  handoff_format: "See handoff-schema.md"
  architecture_context:
    path: "/path/to/.architecture/context.md"  # Absolute path if document exists
    component: "module_name"  # Component/module being implemented
    tier: 0  # 0=foundation, 1=core, 2=application (extracted from context doc)
```

## Step 1b: Junior-Developer Evaluation (STANDARD/EXTENDED mode; SIMPLE mode only when task count >= 6)

The grep-based assignment in Step 1a is a first-pass heuristic. Step 1b is a second-pass evaluation that reads and may update `task-assignments.txt`. Even if Step 1a assigned no tasks to junior-developer, Step 1b may reassign some.

**Evaluation criteria** -- consider junior-developer for a task if any of these are true:
1. Total task count >= 4 (any mode) or >= 6 (SIMPLE mode only)
2. The task has ALL of the following:
   - Single function or single class scope
   - Clear input/output specification derivable from the architecture handoff
   - No cross-module dependencies
   - Does not require design judgment or architectural decisions
3. Architecture identifies utility modules, helper functions, data validation, configuration handling, or test fixtures
4. The task is test-writing separate from implementation logic

**Process**:
1. Count total tasks from `task-assignments.txt`
2. For each task currently assigned to senior-developer, evaluate against criteria above
3. If any tasks qualify AND junior-developer skill is available (`~/.claude/skills/junior-developer/SKILL.md` exists): reassign in `task-assignments.txt`
4. If any tasks qualify BUT junior-developer is unavailable: document "junior-developer skill not available, all tasks retained by senior-developer"
5. Document the evaluation in session state:

```json
{
  "junior_developer_evaluation": {
    "evaluated": true,
    "total_tasks": "<int>",
    "tasks_evaluated_for_junior": "<int>",
    "tasks_reassigned_to_junior": "<int>",
    "reassigned_tasks": [],
    "rationale": "<Why tasks were or were not reassigned>"
  }
}
```

Use jq to write this to session state:
```bash
jq '.junior_developer_evaluation = {"evaluated": true, "total_tasks": '$TOTAL', "tasks_evaluated_for_junior": '$EVALUATED', "tasks_reassigned_to_junior": '$REASSIGNED', "reassigned_tasks": [], "rationale": "'"$RATIONALE"'"}' \
  "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp" && \
  mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
```

**If NO tasks qualify**: Document "No tasks suitable for junior-developer: all require senior judgment" and proceed.
**SIMPLE mode**: Skip this evaluation unless task count >= 6.

## Step 2: Wave-Based Parallel Execution

**SIMPLE mode** (`$PROGRAMMING_PM_MODE` = SIMPLE): skip waves. Walk
`$SESSION_DIR/task-assignments.txt` and dispatch one specialist at a time, recording each
deliverable before moving to the next task.

**STANDARD/EXTENDED mode**: launch in waves as below.

## Wave-Based Specialist Launch

Launch specialists in three waves to respect dependency ordering. Track launched tasks in
`running-agents.txt` to prevent double-launches, and send the agents of a single wave in one
message so they run concurrently.

**Wave 1 (immediate)** -- Launch specialists whose output feeds other tasks:

If mathematician tasks exist in task-assignments.txt:
  For each mathematician task:
    Launch mathematician via Task tool.
    Description: "Mathematician: Design algorithm for {component_name}"
    Prompt: Include algorithm requirements from architecture phase, performance targets, and constraints.
    Write output to: `{session_dir}/deliverables/{task_id}-math-analysis.md`
  Record each launched task in running-agents tracking.

If statistician tasks exist in task-assignments.txt:
  For each statistician task:
    Launch statistician via Task tool.
    Description: "Statistician: Design statistical approach for {component_name}"
    Prompt: Include statistical requirements, data characteristics, and validation criteria.
    Write output to: `{session_dir}/deliverables/{task_id}-stats-analysis.md`
  Record each launched task in running-agents tracking.

These launch first because implementation specialists may depend on their output.

**Wave 2 (after Wave 1 launches)** -- Launch implementation specialists for independent tasks:

Identify tasks with no dependencies (or dependencies already satisfied from prior phases).
For each independent task not already launched in Wave 1:
  Launch senior-developer or junior-developer via Task tool (per task assignment).
  Description: "{Specialist}: Implement {component_name}"
  Prompt: Include architecture spec, coding standards, test requirements, and output path.
  Write output to: `{session_dir}/deliverables/{task_id}-implementation/`
Skip any task already tracked in running-agents (prevents double-launch).

**Wave 3 -- Bounded dependency pass**:

Run this pass after each batch of specialists returns, since a returning specialist may have
unblocked a pending task. It is bounded: a task still blocked after 3 passes escalates.

```bash
MAX_WAVE3_PASSES=3
WAVE3_PASS=$(( $(cat "$SESSION_DIR/wave3-pass-count.txt" 2>/dev/null || echo 0) + 1 ))
echo "$WAVE3_PASS" > "$SESSION_DIR/wave3-pass-count.txt"

for TASK_ID in $(cat "$SESSION_DIR/wave3-pending.txt" 2>/dev/null); do
  # Read dependencies for this task
  DEPS=$(grep "^$TASK_ID|" "$SESSION_DIR/task-assignments.txt" | cut -d'|' -f4 | tr ',' ' ')
  ALL_DEPS_MET=true
  MISSING_DEP=""

  for DEP in $DEPS; do
    # Use ls glob (not [ -f glob ]) for safe wildcard existence check
    if ! ls "${SESSION_DIR}/deliverables/${DEP}-"* > /dev/null 2>&1; then
      ALL_DEPS_MET=false
      MISSING_DEP="$DEP"
      break
    fi
  done

  if [ "$ALL_DEPS_MET" = true ]; then
    echo "  Wave 3: $TASK_ID dependencies met — launch specialist via Task tool dispatch"
  elif [ "$WAVE3_PASS" -ge "$MAX_WAVE3_PASSES" ]; then
    echo "  ESCALATION: $TASK_ID — dependencies unresolved after $MAX_WAVE3_PASSES passes"
    echo "$TASK_ID|ESCALATED|$(date -u +%Y-%m-%dT%H:%M:%SZ)|missing:$MISSING_DEP" >> "$SESSION_DIR/escalations.txt"
  else
    echo "  Wave 3: $TASK_ID still blocked on $MISSING_DEP (pass $WAVE3_PASS/$MAX_WAVE3_PASSES)"
  fi
done

# Present escalations to user if any occurred
if [ -f "$SESSION_DIR/escalations.txt" ] && [ -s "$SESSION_DIR/escalations.txt" ]; then
  ESCALATED_COUNT=$(wc -l < "$SESSION_DIR/escalations.txt" | tr -d ' ')
  echo ""
  echo "HARD ESCALATION: $ESCALATED_COUNT task(s) have unresolvable dependencies."
  echo "Options: (A) Skip blocked tasks, proceed with completed work  (B) Return to Phase 3  (C) Abort"
fi
```

## Step 3: Completion Tracking

Specialists return their results to you directly, so tracking is event-driven — there is no
background poll to run. As each specialist returns:

1. Append `{TASK_ID}|COMPLETED|<UTC timestamp>` to `$SESSION_DIR/task-status.txt` and drop the
   task from `$SESSION_DIR/running-agents.txt`.
2. Confirm the expected deliverable exists under `$SESSION_DIR/deliverables/` (Gate 4b expects
   >100 words).
3. Re-run the Wave 3 dependency pass above and launch anything it unblocked.

If a specialist returns without a usable deliverable or reports that it is blocked, apply the
intervention options from timeout-config.md — narrow the scope, substitute a specialist, or
escalate to the user — rather than relaunching it unchanged.

When every non-escalated task has a deliverable, proceed to Step 4.

## Step 4: Quality Gate 4a - Specialist Completion Check

Validate specialist outputs before proceeding:

```bash
echo "Quality Gate 4a: Specialist Completion Check"

# Count critical vs. implementation specialists
CRITICAL_COUNT=$(grep -cE "mathematician|statistician" "$SESSION_DIR/task-assignments.txt" || echo 0)
IMPL_COUNT=$(grep -cE "senior-developer|junior-developer" "$SESSION_DIR/task-assignments.txt" || echo 0)
TOTAL_COUNT=$((CRITICAL_COUNT + IMPL_COUNT))

# Count completed tasks
COMPLETED_COUNT=$(grep -c "COMPLETED" "$SESSION_DIR/task-status.txt" || echo 0)

echo "Completion status: $COMPLETED_COUNT / $TOTAL_COUNT tasks"

# Decision table
if [ "$COMPLETED_COUNT" -eq "$TOTAL_COUNT" ]; then
  echo "✅ Gate 4a: PASS (100% completion)"
elif [ "$COMPLETED_COUNT" -ge $((TOTAL_COUNT * 3 / 4)) ]; then
  echo "⚠️  Gate 4a: CONDITIONAL PASS (75%+ completion)"
  echo "Note: $((TOTAL_COUNT - COMPLETED_COUNT)) task(s) incomplete"

  # Check if critical specialists completed
  CRITICAL_COMPLETED=$(grep -E "TASK-.*mathematician|TASK-.*statistician" "$SESSION_DIR/task-status.txt" | grep -c "COMPLETED" || echo 0)

  if [ "$CRITICAL_COMPLETED" -eq "$CRITICAL_COUNT" ]; then
    echo "→ All critical specialists completed. Proceeding with note."
  else
    echo "→ Critical specialists incomplete. RETRY required."
    exit 1
  fi
else
  echo "❌ Gate 4a: FAIL (<75% completion)"
  echo "→ RETRY required"
  exit 1
fi
```

**Quality Gate 4b: Implementation Validation**:
- Type: Automated (output validation)
- Criteria:
  - [ ] All specialist outputs exist and are >100 words
  - [ ] No critical blocking issues flagged
  - [ ] Handoffs validate against schema (validate-handoff.py)
  - [ ] Acceptance criteria met per task
- Pass Condition: All criteria checked OR 75%+ with critical specialists complete
- Fail Action: Retry incomplete tasks or escalate to user

**Handoff Validation** (Phase 4 → Phase 5):
```bash
# Validate all code handoffs from Phase 4 (may be multiple task handoffs)
echo "Validating Phase 4 code handoffs..."

VALIDATION_ERRORS=0

for HANDOFF_FILE in "${SESSION_DIR}/handoffs/phase4-"*"-handoff-"*.yaml; do
  [ ! -f "$HANDOFF_FILE" ] && continue

  # Determine handoff type based on filename
  if echo "$HANDOFF_FILE" | grep -q "math-handoff"; then
    HANDOFF_TYPE="math_handoff"
  elif echo "$HANDOFF_FILE" | grep -q "stats-handoff"; then
    HANDOFF_TYPE="stats_handoff"
  elif echo "$HANDOFF_FILE" | grep -q "code-handoff"; then
    HANDOFF_TYPE="code_handoff"
  else
    echo "⚠️  Unknown handoff type: $HANDOFF_FILE"
    continue
  fi

  echo "  Validating $(basename "$HANDOFF_FILE") as $HANDOFF_TYPE..."

  python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
    "$HANDOFF_FILE" \
    "$HANDOFF_TYPE"

  if [ $? -ne 0 ]; then
    echo "  ❌ Validation failed for $(basename "$HANDOFF_FILE")"
    VALIDATION_ERRORS=$((VALIDATION_ERRORS + 1))
  else
    echo "  ✅ Validated successfully"
  fi
done

if [ "$VALIDATION_ERRORS" -gt 0 ]; then
  echo ""
  echo "❌ Phase 4 handoff validation FAILED ($VALIDATION_ERRORS error(s))"
  read -p "Fix issues and retry? [Y/n]: " RETRY_CHOICE

  if [ "$RETRY_CHOICE" != "n" ] && [ "$RETRY_CHOICE" != "N" ]; then
    echo "Returning to Phase 4..."
    exit 1
  else
    echo "⚠️  Override: Proceeding with validation errors (documented)"
    jq --arg count "$VALIDATION_ERRORS" '.phase4_handoff_override = true | .phase4_validation_errors = ($count | tonumber)' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ All Phase 4 handoffs validated successfully"
fi
```

**Phase Transition**: Phase 4 complete -> Quality Gate 4 -> PROCEED to Phase 5: Code Review and Testing
