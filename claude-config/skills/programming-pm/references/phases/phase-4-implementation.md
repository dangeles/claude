# Phase 4: Implementation — Full Detail

This reference holds the complete implementation walkthrough for Phase 4 of the
programming-pm workflow. SKILL.md keeps the phase summary, gate criteria, and a
pointer here. Read this file when you are actively executing Phase 4.

Before starting Phase 4: Read `~/.claude/programming-pm-sessions/{workflow-id}/state.yaml`. Confirm Phases 0-3 are complete.

**Objective**: Implement architecture with specialist agents in parallel.

**Mode-based execution**:
- **SIMPLE**: Sequential execution (one specialist at a time)
- **STANDARD/EXTENDED**: Wave-based parallel execution (waves at T=0s, T=30s, T=60s)

**Steps**:

## Step 0: Architecture Pre-flight Validation

**Skip condition**: SIMPLE mode projects with <= 1 component skip Step 0 and proceed directly to Step 1.

Before any task decomposition, validate the architecture handoff for implementability.

**Owner**: programming-pm dispatches to senior-developer via Task tool.

**Context size check**: Before filling the dispatch:
- If architecture handoff YAML is <= 200 lines: paste verbatim
- If > 200 lines: pass file path and instruct senior-developer to read it directly

**Dispatch** (using template below): Send architecture handoff + requirements summary to senior-developer for pre-flight review.

**senior-developer validates**:
1. All component interfaces are fully specified (inputs have types, outputs have types)
2. Component dependencies are traversable without cycles (check by tracing dependency chains)
3. Technology choices are recognizable and compatible
4. Implementation order in the handoff is internally consistent

**senior-developer returns** a validation report (write to `{SESSION_DIR}/deliverables/phase4-preflight.yaml`):
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

**Timeout**: 15 minutes. If Step 0 Task tool invocation fails or times out: log "Step 0 pre-flight unavailable" and proceed to Step 1 with a warning.

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

## Step 1b: Junior-Developer Evaluation (Mandatory for STANDARD/EXTENDED mode; optional for SIMPLE mode unless task count >= 6)

The grep-based assignment in Step 1a is a first-pass heuristic. Step 1b is a mandatory second-pass evaluation that reads and may update `task-assignments.txt`. Even if Step 1a assigned no tasks to junior-developer, Step 1b may reassign some.

**Evaluation criteria** -- evaluate junior-developer for a task if ANY of these are true:
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

**SIMPLE mode**: Skip waves, execute sequentially.

**STANDARD/EXTENDED mode**: Execute in waves with stagger.

```bash
# Check execution mode
if [ "$PROGRAMMING_PM_MODE" = "SIMPLE" ]; then
  echo "SIMPLE mode: Sequential execution"

  # Execute tasks one at a time
  while IFS='|' read -r TASK_ID COMPONENT SPECIALIST DEPS; do
    echo "Executing $TASK_ID ($COMPONENT) with $SPECIALIST..."

    # Invoke specialist (synchronous)
    # Record start time for timeout monitoring
    START_TIME=$(date +%s)

    # ... invoke specialist ...

    # Wait for completion
  done < "$SESSION_DIR/task-assignments.txt"

else
  echo "STANDARD/EXTENDED mode: Wave-based parallel execution"
fi
```

## Wave-Based Specialist Launch

Launch specialists in three waves to respect dependency ordering. Track all running agents to prevent double-launches.

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

**Wave 3 -- Bounded dependency retry protocol**:

```bash
MAX_WAVE3_RETRIES=3
WAVE3_RETRY_INTERVAL=120  # seconds

# For each Wave 3 task with unsatisfied dependencies
for TASK_ID in $(cat "$SESSION_DIR/wave3-pending.txt" 2>/dev/null); do
  RETRY_COUNT=0
  DEPS_SATISFIED=false

  while [ "$RETRY_COUNT" -lt "$MAX_WAVE3_RETRIES" ] && [ "$DEPS_SATISFIED" = false ]; do
    RETRY_COUNT=$((RETRY_COUNT + 1))
    MISSING_DEP=""

    # Read dependencies for this task
    DEPS=$(grep "^$TASK_ID|" "$SESSION_DIR/task-assignments.txt" | cut -d'|' -f4 | tr ',' ' ')
    ALL_DEPS_MET=true

    for DEP in $DEPS; do
      # Use ls glob (not [ -f glob ]) for safe wildcard existence check
      if ! ls "${SESSION_DIR}/deliverables/${DEP}-"* > /dev/null 2>&1; then
        ALL_DEPS_MET=false
        MISSING_DEP="$DEP"
        break
      fi
    done

    if [ "$ALL_DEPS_MET" = true ]; then
      DEPS_SATISFIED=true
      # Launch specialist via Task tool dispatch
    else
      echo "  Wave 3: $TASK_ID retry $RETRY_COUNT/$MAX_WAVE3_RETRIES (waiting on: $MISSING_DEP)"
      sleep "$WAVE3_RETRY_INTERVAL"
    fi
  done

  if [ "$DEPS_SATISFIED" = false ]; then
    echo "  ESCALATION: $TASK_ID — dependencies unresolved after $MAX_WAVE3_RETRIES retries"
    echo "$TASK_ID|ESCALATED|$(date -u +%Y-%m-%dT%H:%M:%SZ)|missing:$MISSING_DEP" >> "$SESSION_DIR/escalations.txt"
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

Monitor all agents via output file existence. When all non-escalated tasks show deliverables, proceed to quality gate.

## Step 3: Progress Monitoring

Monitor specialist outputs using file-based tracking:

```bash
# Progress monitoring loop
echo "Monitoring specialist progress..."

TIMEOUT_THRESHOLD=7200  # 2 hours (STANDARD mode)
if [ "$PROGRAMMING_PM_MODE" = "EXTENDED" ]; then
  TIMEOUT_THRESHOLD=14400  # 4 hours (EXTENDED mode)
fi

# Derive cycle count from TIMEOUT_THRESHOLD (single source of truth for timeout duration)
MAX_MONITORING_CYCLES=$(( TIMEOUT_THRESHOLD / 60 ))  # TIMEOUT_THRESHOLD defined in timeout-config.md
MONITORING_CYCLE=0

while [ "$MONITORING_CYCLE" -lt "$MAX_MONITORING_CYCLES" ]; do
  MONITORING_CYCLE=$((MONITORING_CYCLE + 1))

  # Check running agents
  RUNNING_COUNT=$(wc -l < "$SESSION_DIR/running-agents.txt" 2>/dev/null || echo 0)

  if [ "$RUNNING_COUNT" -eq 0 ]; then
    echo "✅ All specialists completed"
    break
  fi

  # Check each running agent
  while IFS='|' read -r TASK_ID AGENT_PID; do
    # Check if process still running
    if ! ps -p "$AGENT_PID" > /dev/null 2>&1; then
      echo "  $TASK_ID completed (PID $AGENT_PID exited)"

      # Mark as completed
      echo "$TASK_ID|COMPLETED|$(date -u +"%Y-%m-%dT%H:%M:%SZ")" >> "$SESSION_DIR/task-status.txt"

      # Remove from running list
      grep -v "^$TASK_ID|" "$SESSION_DIR/running-agents.txt" > "$SESSION_DIR/running-agents.txt.tmp"
      mv "$SESSION_DIR/running-agents.txt.tmp" "$SESSION_DIR/running-agents.txt"
    else
      # Check for timeout
      START_TIME=$(grep "^$TASK_ID|" "$SESSION_DIR/task-start-times.txt" | cut -d'|' -f2)
      CURRENT_TIME=$(date +%s)
      ELAPSED=$((CURRENT_TIME - START_TIME))

      if [ "$ELAPSED" -gt "$TIMEOUT_THRESHOLD" ]; then
        echo "  ⚠️  $TASK_ID TIMEOUT (elapsed: ${ELAPSED}s, threshold: ${TIMEOUT_THRESHOLD}s)"

        # Timeout intervention (see timeout-config.md)
        # Option: Extend deadline, narrow scope, substitute specialist, or escalate
      fi
    fi
  done < "$SESSION_DIR/running-agents.txt"

  # Check progress file outputs (must be >100 words)
  for TASK_ID in $(awk -F'|' '{print $1}' "$SESSION_DIR/task-assignments.txt"); do
    PROGRESS_FILE="/tmp/progress-${TASK_ID}.md"

    if [ -f "$PROGRESS_FILE" ]; then
      WORD_COUNT=$(wc -w < "$PROGRESS_FILE")

      if [ "$WORD_COUNT" -ge 100 ]; then
        echo "  $TASK_ID progress OK ($WORD_COUNT words)"
      else
        echo "  $TASK_ID progress insufficient ($WORD_COUNT words, min 100)"
      fi
    fi
  done

  # Sleep before next check
  sleep 60  # Check every minute
done

# Hard abort if loop exhausted with tasks still running
if [ "${RUNNING_COUNT:-0}" -gt 0 ]; then
  echo "HARD ABORT: Monitoring exhausted ($MAX_MONITORING_CYCLES cycles × 60s = $TIMEOUT_THRESHOLD s)"
  echo "Still running: $RUNNING_COUNT task(s)"
  echo "Options: (A) Proceed with completed tasks  (B) Extend monitoring (+60 cycles)  (C) Abort workflow"
fi
```

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
