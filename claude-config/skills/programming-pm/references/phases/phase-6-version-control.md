# Phase 6: Version Control Integration — Full Detail

This reference holds the complete VCS integration walkthrough (Steps 1-5) for
Phase 6 of the programming-pm workflow. SKILL.md keeps the phase summary, gate
criteria, optional shortcuts/advisory notes, and a pointer here. Read this file
when actively executing Phase 6.

Before starting Phase 6: Read `~/.claude/programming-pm-sessions/{workflow-id}/state.yaml`. Confirm Phases 0-5 are complete.

**Objective**: Integrate changes with sync-config.py and version control.

## Optional: Commit Commands Shortcut

For straightforward projects where no branching complexity is involved, the `commit-commands` plugin provides a single-action alternative: `/commit-push-pr` handles commit + push + PR creation in one step. Use this instead of the manual Steps 2–3 below when the git strategy is simple.

## Optional: Git Strategy Advisory

Before proceeding with version control integration, you MAY invoke `git-strategy-advisor`
via Task tool in post-work mode to get scope-adaptive git recommendations:

**Invocation** (via Task tool):
```
Use git-strategy-advisor to determine git strategy for completed work.

mode: post-work
```

The advisor analyzes actual changes (files, lines, directories) and recommends branch
strategy, branch naming, push timing, and PR creation.

**Conflict resolution**: If the advisor's recommendation differs from Phase 6 Step 3's
existing logic (which creates a feature branch when on main), Phase 6 Step 3 takes
precedence unconditionally. Present the advisor's recommendation as an informational
note in the completion summary (e.g., "Note: git-strategy-advisor suggests direct-commit
for this trivial change"). The orchestrator proceeds with its default behavior.

**Response handling**: Read the advisor's `summary` field for the human-readable
recommendation. Optionally read `strategy.branch.action` to note whether the advisor
agrees with the default strategy. Include the summary in the Phase 6 completion report.

**Confidence handling**: If the advisor returns confidence "none" (e.g., no git repository
found), silently skip the git strategy section. If confidence is "low", present the
recommendation with a caveat noting limited accuracy.

This is **advisory only**. If `git-strategy-advisor` is not available or returns an error,
proceed with existing Phase 6 logic unchanged.

**Steps**:

## Step 1: Pre-Merge Validation

```bash
# Check sync-config.py availability
SYNC_CONFIG_PATH="/Users/davidangelesalbores/repos/claude/sync-config.py"

if [ -f "$SYNC_CONFIG_PATH" ]; then
  echo "Using sync-config.py for VCS integration"
  USE_SYNC_CONFIG=true
else
  echo "⚠️  sync-config.py not found. Falling back to direct git commands."
  USE_SYNC_CONFIG=false
fi

# If using sync-config.py, check status
if [ "$USE_SYNC_CONFIG" = true ]; then
  echo "Checking sync-config.py status..."

  "$SYNC_CONFIG_PATH" status

  if [ $? -ne 0 ]; then
    echo "⚠️  Sync status check failed. Proceeding with caution."
  fi
fi

# Validate all handoffs one final time
echo "Final handoff validation before merge..."

FINAL_VALIDATION_ERRORS=0

for HANDOFF_FILE in "${SESSION_DIR}/handoffs/"*.yaml; do
  [ ! -f "$HANDOFF_FILE" ] && continue

  HANDOFF_NAME=$(basename "$HANDOFF_FILE" .yaml)

  # Determine handoff type from filename
  if echo "$HANDOFF_NAME" | grep -q "session-handoff"; then
    HANDOFF_TYPE="session_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "requirements-handoff"; then
    HANDOFF_TYPE="requirements_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "premortem-handoff"; then
    HANDOFF_TYPE="premortem_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "architecture-handoff"; then
    HANDOFF_TYPE="architecture_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "math-handoff"; then
    HANDOFF_TYPE="math_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "stats-handoff"; then
    HANDOFF_TYPE="stats_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "code-handoff"; then
    HANDOFF_TYPE="code_handoff"
  elif echo "$HANDOFF_NAME" | grep -q "review-handoff"; then
    HANDOFF_TYPE="review_handoff"
  else
    echo "⚠️  Unknown handoff type: $HANDOFF_NAME"
    continue
  fi

  python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
    "$HANDOFF_FILE" \
    "$HANDOFF_TYPE" > /dev/null 2>&1

  if [ $? -ne 0 ]; then
    FINAL_VALIDATION_ERRORS=$((FINAL_VALIDATION_ERRORS + 1))
  fi
done

if [ "$FINAL_VALIDATION_ERRORS" -gt 0 ]; then
  echo "⚠️  $FINAL_VALIDATION_ERRORS handoff(s) still have validation issues"
  echo "Review overrides in session-state.json"
fi

# Dry-run to detect conflicts (if using sync-config.py)
if [ "$USE_SYNC_CONFIG" = true ]; then
  echo "Running sync dry-run to detect conflicts..."

  "$SYNC_CONFIG_PATH" push --dry-run

  if [ $? -ne 0 ]; then
    echo "❌ Dry-run detected conflicts or issues"
    echo "Review and resolve before proceeding."
    read -p "Continue anyway? [y/N]: " CONTINUE_CHOICE

    if [ "$CONTINUE_CHOICE" != "y" ] && [ "$CONTINUE_CHOICE" != "Y" ]; then
      echo "Aborting merge. Resolve conflicts first."
      exit 1
    fi
  else
    echo "✅ Dry-run successful (no conflicts detected)"
  fi
fi
```

## Step 2: Quality Gate 6 Validation

```bash
# Run Quality Gate 6 validation script
echo "Running Quality Gate 6 validation..."

"${SKILL_DIR}/scripts/validate-gate.sh" 6 \
  "${SESSION_DIR}/handoffs/phase5-review-handoff.yaml" \
  "${SESSION_DIR}"

if [ $? -ne 0 ]; then
  echo "❌ Quality Gate 6 FAILED"
  read -p "Override and proceed? [y/N]: " OVERRIDE_CHOICE

  if [ "$OVERRIDE_CHOICE" != "y" ] && [ "$OVERRIDE_CHOICE" != "Y" ]; then
    echo "Aborting. Fix Quality Gate 6 issues first."
    exit 1
  else
    echo "⚠️  GATE 6 OVERRIDE (logged)"
    jq '.gate6_override = true | .gate6_override_timestamp = now | .gate6_override_user = env.USER' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ Quality Gate 6 PASSED"
fi
```

**Quality Gate 6 Criteria**:
- Type: Automated (VCS checks)
- Criteria:
  - [ ] All previous gates passed (or overrides documented)
  - [ ] No merge conflicts (verified by dry-run)
  - [ ] Review approved (from phase5-review-handoff.yaml)
  - [ ] Deliverable location documented
  - [ ] Files staged (if in git repo)
- Override: Repository admin can force merge (logged for audit)

## Step 3: Commit and Sync

```bash
# Create feature branch if needed
if git rev-parse --git-dir > /dev/null 2>&1; then
  CURRENT_BRANCH=$(git branch --show-current)

  if [ "$CURRENT_BRANCH" = "main" ] || [ "$CURRENT_BRANCH" = "master" ]; then
    echo "⚠️  On main branch. Creating feature branch..."

    # Generate branch name from requirements
    BRANCH_NAME="feature/programming-pm-$(date +%Y%m%d-%H%M%S)"

    git checkout -b "$BRANCH_NAME"
    echo "Created branch: $BRANCH_NAME"
  fi
fi

# Stage specific files (NEVER git add . or git add -A)
echo "Staging specific files..."

# Get changed files from handoff
if command -v yq &> /dev/null; then
  CHANGED_FILES=$(yq eval '.handoff.changes.files_changed[].path' \
    "${SESSION_DIR}/handoffs/phase5-review-handoff.yaml" 2>/dev/null)

  if [ -n "$CHANGED_FILES" ]; then
    echo "$CHANGED_FILES" | while read -r FILE; do
      if [ -f "$FILE" ]; then
        git add "$FILE"
        echo "  Staged: $FILE"
      fi
    done
  else
    echo "⚠️  No files listed in review handoff. Manual staging required."
  fi
else
  echo "⚠️  yq not found. Manual staging required."
  git status
fi

# Create commit with conventional format
PROBLEM_STATEMENT=$(yq eval '.handoff.requirements.problem_statement' \
  "${SESSION_DIR}/handoffs/phase1-requirements-handoff.yaml" 2>/dev/null | head -n1)

COMMIT_MESSAGE="feat(programming-pm): ${PROBLEM_STATEMENT}

Implemented via programming-pm workflow ($(date -u +"%Y-%m-%d"))

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"

git commit -m "$COMMIT_MESSAGE"

if [ $? -ne 0 ]; then
  echo "❌ Commit failed"
  exit 1
else
  echo "✅ Commit created successfully"
fi

# Sync to ~/.claude/ (if using sync-config.py)
if [ "$USE_SYNC_CONFIG" = true ]; then
  echo "Syncing changes to ~/.claude/..."

  "$SYNC_CONFIG_PATH" push

  if [ $? -ne 0 ]; then
    echo "❌ sync-config.py push failed"
    echo "Changes committed to git but not synced to ~/.claude/"
    echo "Run manually: $SYNC_CONFIG_PATH push"
  else
    echo "✅ Changes synced to ~/.claude/ successfully"
  fi
else
  echo "⚠️  Skipping sync-config.py push (not available)"
fi
```

## Step 4: Create Planning Journal Entry

```bash
# Create planning journal entry documenting the workflow execution
if [ "$USE_SYNC_CONFIG" = true ]; then
  echo "Creating planning journal entry..."

  # Extract brief description from problem statement (first 60 chars)
  BRIEF_DESC=$(echo "$PROBLEM_STATEMENT" | cut -c1-60)

  "$SYNC_CONFIG_PATH" plan --title "$BRIEF_DESC"

  if [ $? -ne 0 ]; then
    echo "⚠️  Failed to create planning journal entry"
    echo "Create manually with: $SYNC_CONFIG_PATH plan --title '...'"
  else
    echo "Document the following in the journal entry:"
    echo "  - Objective: $PROBLEM_STATEMENT"
    echo "  - Specialists used: $(awk -F'|' '{print $3}' "${SESSION_DIR}/task-assignments.txt" | sort -u | tr '\n' ', ')"
    echo "  - Files changed: $(git diff --name-only HEAD~1 HEAD | wc -l | tr -d ' ')"
    echo "  - Testing: All quality gates passed"
    echo "  - Outcome: Success"
    echo ""
    read -p "Press Enter after documenting in journal..."
  fi
fi
```

## Step 5: Session Cleanup

```bash
# Mark session as completed
jq '.status = "completed" | .completion_timestamp = now' \
  "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"

# Determine if session should be deleted or preserved
SESSION_SUCCESSFUL=true

# Check for any gate overrides
if jq -e '.phase0_handoff_override or .phase1_handoff_override or .phase2_handoff_override or .phase3_handoff_override or .phase4_handoff_override or .phase5_handoff_override or .gate6_override' \
  "${SESSION_DIR}/session-state.json" > /dev/null 2>&1; then
  SESSION_SUCCESSFUL=false
  echo "⚠️  Session had overrides. Preserving for review."
fi

# Check for validation errors
if [ "$FINAL_VALIDATION_ERRORS" -gt 0 ]; then
  SESSION_SUCCESSFUL=false
  echo "⚠️  Session had validation errors. Preserving for review."
fi

# Delete or preserve session directory
if [ "$SESSION_SUCCESSFUL" = true ]; then
  echo "Session completed successfully. Cleaning up..."
  echo "Session directory: ${SESSION_DIR}"
  read -p "Delete session directory? [Y/n]: " DELETE_CHOICE

  if [ "$DELETE_CHOICE" != "n" ] && [ "$DELETE_CHOICE" != "N" ]; then
    rm -rf "${SESSION_DIR}"
    echo "✅ Session directory deleted"
  else
    echo "Session directory preserved: ${SESSION_DIR}"
  fi
else
  echo "⚠️  Session preserved for debugging: ${SESSION_DIR}"
  echo "Review session-state.json for overrides and validation errors."
fi

echo ""
echo "================================================"
echo "  Programming-PM Workflow Complete"
echo "================================================"
echo ""
```

**Post-Merge Verification**:
After sync, prompt user to verify deliverable meets expectations. If issues found, create follow-up task (not rollback unless critical).
