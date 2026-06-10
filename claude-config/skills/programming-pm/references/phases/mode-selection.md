# Mode Selection (After Phase 1) — Full Detail

This reference holds the complete Mode Selection walkthrough (the 5-step bash
procedure) for the programming-pm workflow. SKILL.md keeps the mode summary,
tier triggers, and a pointer here. Read this file when actively running mode
selection after Quality Gate 1.

**Objective**: Select workflow execution mode based on project complexity.

**Trigger**: After Quality Gate 1 passes (requirements approved)

**Three execution modes**:
- **SIMPLE** (~1-2 hrs): Single component, no stats/math, <5 implementation tasks
- **STANDARD** (~4-6 hrs): Multi-component (2-5), optional stats/math, 5-15 tasks (default)
- **EXTENDED** (~8-12 hrs): >5 components OR both stats+math OR >15 tasks OR architectural complexity

**Steps**:

## Step 1: Run Complexity Detection

```bash
# Source the detection function
source "${SKILL_DIR}/references/mode-selection-criteria.md"

# Run detection on requirements handoff
REQUIREMENTS_FILE="${SESSION_DIR}/handoffs/phase1-requirements-handoff.yaml"
DETECTION_RESULT=$(detect_tier "$REQUIREMENTS_FILE")

# Parse results
DETECTED_TIER=$(echo "$DETECTION_RESULT" | cut -d'|' -f1)
CONFIDENCE=$(echo "$DETECTION_RESULT" | cut -d'|' -f2)
REASON=$(echo "$DETECTION_RESULT" | cut -d'|' -f3)
```

**Triggers for each tier**:

**EXTENDED**:
- Component count >5
- Requires BOTH statistics AND mathematics
- Task count >15
- Architectural complexity keywords: "distributed system", "microservices", "event-driven", "real-time processing"
- User explicit request: "extended analysis" or "comprehensive review"

**SIMPLE**:
- Single component AND no stats/math
- Utility script with <5 tasks
- Data pipeline (ETL) with single component

**STANDARD** (default):
- Multiple components (2-5)
- Single specialization (stats OR math, not both)
- Moderate task count (5-15)
- Standard patterns: "web API", "CLI tool", "data analysis", "visualization"

## Step 2: Display Mode Selection Prompt

```bash
echo ""
echo "================================================"
echo "  Mode Selection (After Quality Gate 1)"
echo "================================================"
echo ""
echo "Detected tier: $DETECTED_TIER (confidence: $CONFIDENCE)"
echo "Reason: $REASON"
echo ""
echo "Mode descriptions:"
echo "  SIMPLE (1-2 hrs): Single component, no stats/math, <5 tasks"
echo "  STANDARD (4-6 hrs): Multi-component, optional stats/math, 5-15 tasks (default)"
echo "  EXTENDED (8-12 hrs): >5 components OR both stats+math OR >15 tasks"
echo ""
```

## Step 3: User Override Confirmation

```bash
# High-confidence: allow override with 60s timeout
if [ "$CONFIDENCE" = "high" ]; then
  read -t 60 -p "Proceed with $DETECTED_TIER mode? [Y/n]: " USER_CHOICE

  if [ $? -ne 0 ]; then
    # Timeout - proceed with detected tier
    echo "No response (timeout 60s). Proceeding with: $DETECTED_TIER"
    SELECTED_TIER="$DETECTED_TIER"
  elif [ "$USER_CHOICE" = "n" ] || [ "$USER_CHOICE" = "N" ]; then
    # User wants override
    read -p "Select mode (1=SIMPLE, 2=STANDARD, 3=EXTENDED): " MODE_OVERRIDE
    case "$MODE_OVERRIDE" in
      1) SELECTED_TIER="SIMPLE" ;;
      2) SELECTED_TIER="STANDARD" ;;
      3) SELECTED_TIER="EXTENDED" ;;
      *) SELECTED_TIER="$DETECTED_TIER" ;;
    esac

    # Risky override confirmation
    if [ "$DETECTED_TIER" != "SIMPLE" ] && [ "$SELECTED_TIER" = "SIMPLE" ]; then
      echo "⚠️  WARNING: Selecting SIMPLE when $DETECTED_TIER recommended."
      read -p "Confirm risky override? [y/N]: " RISKY_CONFIRM
      if [ "$RISKY_CONFIRM" != "y" ]; then
        SELECTED_TIER="$DETECTED_TIER"
      fi
    fi
  else
    SELECTED_TIER="$DETECTED_TIER"
  fi
else
  # Medium/low confidence: require user confirmation
  read -p "Select mode (1=SIMPLE, 2=STANDARD, 3=EXTENDED) [default: $DETECTED_TIER]: " USER_CHOICE
  case "$USER_CHOICE" in
    1) SELECTED_TIER="SIMPLE" ;;
    2) SELECTED_TIER="STANDARD" ;;
    3) SELECTED_TIER="EXTENDED" ;;
    "") SELECTED_TIER="$DETECTED_TIER" ;;
    *) SELECTED_TIER="STANDARD" ;;  # Safest default
  esac
fi
```

## Step 4: Record Mode Selection

```bash
# Create mode-selection.json
cat > "$SESSION_DIR/mode-selection.json" <<EOF
{
  "detected_tier": "$DETECTED_TIER",
  "confidence": "$CONFIDENCE",
  "reason": "$REASON",
  "selected_tier": "$SELECTED_TIER",
  "override": $([ "$DETECTED_TIER" != "$SELECTED_TIER" ] && echo "true" || echo "false"),
  "timestamp": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
}
EOF

# Update session-state.json
if command -v jq &> /dev/null; then
  jq --arg tier "$SELECTED_TIER" '.mode = $tier' \
    "$SESSION_DIR/session-state.json" > "$SESSION_DIR/session-state.json.tmp"
  mv "$SESSION_DIR/session-state.json.tmp" "$SESSION_DIR/session-state.json"
fi

# Export for workflow use
export PROGRAMMING_PM_MODE="$SELECTED_TIER"

echo "✅ Mode selected: $SELECTED_TIER"
```

## Step 5: Mode-Based Branching

```bash
# Workflow branching based on selected mode
if [ "$PROGRAMMING_PM_MODE" = "SIMPLE" ]; then
  echo "→ SIMPLE mode: Sequential execution, automated checks only"
  SKIP_EXTENDED_ANALYSIS=true
  PARALLEL_EXECUTION=false
elif [ "$PROGRAMMING_PM_MODE" = "EXTENDED" ]; then
  echo "→ EXTENDED mode: Wave-based parallel execution, extended reviews"
  SKIP_EXTENDED_ANALYSIS=false
  PARALLEL_EXECUTION=true
  EXTENDED_TIMEOUTS=true
else
  echo "→ STANDARD mode: Wave-based parallel execution, standard checks"
  SKIP_EXTENDED_ANALYSIS=false
  PARALLEL_EXECUTION=true
  EXTENDED_TIMEOUTS=false
fi
```

**Backwards Compatibility**:
```bash
# For sessions without mode-selection.json, default to STANDARD
if [ ! -f "$SESSION_DIR/mode-selection.json" ]; then
  echo "⚠️  Legacy session (no mode selection). Defaulting to STANDARD."
  export PROGRAMMING_PM_MODE="STANDARD"
fi
```

**Phase Transition**: Phase 1 complete -> Quality Gate 1 (user approval required) -> PROCEED to Phase 2: Pre-Mortem and Risk Assessment
