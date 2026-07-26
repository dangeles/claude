# Handoff Validation Gates

Validation scripts run at each phase transition of the programming-pm workflow. SKILL.md
names the handoff file and schema type per transition and points here for the script and
its failure/override handling. Run these yourself with the Bash tool.

Two transitions have their scripts in their phase references instead:
Phase 0 → 1 in `phases/phase-0-session-setup.md`, Phase 4 → 5 in
`phases/phase-4-implementation.md`.

## Contents

- Phase 1 → Phase 2 (requirements_handoff)
- Phase 2 → Phase 3 (premortem_handoff)
- Phase 3 → Phase 4 (architecture_handoff)
- Phase 5 → Phase 6 (review_handoff)

## Phase 1 → Phase 2

```bash
# Validate requirements handoff before mode selection
python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
  "${SESSION_DIR}/handoffs/phase1-requirements-handoff.yaml" \
  "requirements_handoff"

if [ $? -ne 0 ]; then
  echo "❌ Phase 1 handoff validation FAILED"
  echo "Options:"
  echo "  (A) Return to requirements-analyst to fix issues"
  echo "  (B) Override with documented gaps"
  read -p "Choice [A/b]: " OVERRIDE_CHOICE

  if [ "$OVERRIDE_CHOICE" != "b" ] && [ "$OVERRIDE_CHOICE" != "B" ]; then
    echo "Returning to requirements-analyst..."
    exit 1
  else
    echo "⚠️  Override: Proceeding with gaps (documented)"
    jq '.phase1_handoff_override = true' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ Phase 1 handoff validated successfully"
fi
```

## Phase 2 → Phase 3

```bash
# Validate pre-mortem handoff
python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
  "${SESSION_DIR}/handoffs/phase2-premortem-handoff.yaml" \
  "premortem_handoff"

if [ $? -ne 0 ]; then
  echo "❌ Phase 2 handoff validation FAILED"
  read -p "Fix issues and retry? [Y/n]: " RETRY_CHOICE

  if [ "$RETRY_CHOICE" != "n" ] && [ "$RETRY_CHOICE" != "N" ]; then
    exit 1
  else
    echo "⚠️  Override: Proceeding with validation gaps"
    jq '.phase2_handoff_override = true' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ Phase 2 handoff validated successfully"
fi
```

## Phase 3 → Phase 4

```bash
# Validate architecture handoff before implementation
python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
  "${SESSION_DIR}/handoffs/phase3-architecture-handoff.yaml" \
  "architecture_handoff"

if [ $? -ne 0 ]; then
  echo "❌ Phase 3 handoff validation FAILED"
  echo "Incomplete architecture cannot proceed to implementation."
  read -p "Return to systems-architect? [Y/n]: " RETRY_CHOICE

  if [ "$RETRY_CHOICE" != "n" ] && [ "$RETRY_CHOICE" != "N" ]; then
    exit 1
  else
    echo "⚠️  Override: Proceeding with incomplete architecture (HIGH RISK)"
    jq '.phase3_handoff_override = true | .phase3_override_risk = "HIGH"' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ Phase 3 handoff validated successfully"
fi
```

## Phase 5 → Phase 6

```bash
# Validate review handoff before VCS integration
python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
  "${SESSION_DIR}/handoffs/phase5-review-handoff.yaml" \
  "review_handoff"

if [ $? -ne 0 ]; then
  echo "❌ Phase 5 handoff validation FAILED"
  echo "Code review handoff incomplete. Cannot proceed to merge."
  read -p "Return to code review? [Y/n]: " RETRY_CHOICE

  if [ "$RETRY_CHOICE" != "n" ] && [ "$RETRY_CHOICE" != "N" ]; then
    exit 1
  else
    echo "⚠️  Override: Proceeding without complete review (CRITICAL RISK)"
    jq '.phase5_handoff_override = true | .phase5_override_risk = "CRITICAL"' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ Phase 5 handoff validated successfully"
fi
```
