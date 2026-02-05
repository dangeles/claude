# Mode Selection Criteria for programming-pm

## Purpose

Determine workflow execution mode (Simple, Standard, Extended) based on project complexity.

**Primary goal**: Route projects to appropriate execution mode:
- **SIMPLE**: Fast path for single-component, no stats/math projects (~1-2 hours)
- **STANDARD**: Full workflow for typical multi-component projects (~4-6 hours, default)
- **EXTENDED**: Standard + extended analysis for complex architectural projects (~8-12 hours)

**Target metrics**:
- False positive rate: <5% (simple projects incorrectly flagged as extended)
- False negative rate: <5% (extended projects incorrectly flagged as simple)
- User override rate: <20% (indicates good auto-detection)

---

## Detection Criteria

### High-Confidence Extended (Always Triggers Extended Mode)

1. **Many components** (>5 components)
   - Component count exceeds 5
   - Rationale: Multi-component architectures require extended coordination and testing

2. **Requires both statistics AND mathematics**
   - Requirements mention: statistical methods AND algorithm design
   - Rationale: Dual specialization increases complexity and coordination overhead

3. **Large task count** (>15 implementation tasks)
   - Task decomposition yields >15 tasks
   - Rationale: High task count indicates architectural complexity

4. **Architectural complexity keywords**
   - Requirements contain: "distributed system", "microservices", "event-driven", "real-time processing"
   - Rationale: These patterns require extended architectural assessment

5. **Explicit user request**
   - Requirements contain: "extended analysis" or "comprehensive review"
   - Rationale: User explicitly wants thorough assessment

---

### High-Confidence Simple (Always Triggers Simple Mode)

1. **Single component, no specialization**
   - Architecture has 1 component
   - AND: No statistical methods
   - AND: No algorithm design
   - Rationale: Straightforward implementation with minimal coordination

2. **Utility script or tool**
   - Requirements mention: "utility", "script", "helper", "tool"
   - AND: <5 implementation tasks
   - Rationale: Small-scale tools don't require full workflow

3. **Data pipeline with clear structure**
   - Requirements mention: "ETL", "data pipeline", "batch processing"
   - AND: Single-component architecture
   - Rationale: Well-understood pattern with minimal architectural decisions

---

### High-Confidence Standard (Default Mode)

1. **Multiple components (2-5)**
   - Component count between 2 and 5
   - Rationale: Typical multi-component project

2. **Single specialization**
   - Requires EITHER statistics OR mathematics, not both
   - Rationale: Single specialist reduces coordination overhead

3. **Moderate task count (5-15 tasks)**
   - Task decomposition yields 5-15 tasks
   - Rationale: Balanced between simple and extended

4. **Standard patterns**
   - Requirements mention: "web API", "CLI tool", "data analysis", "visualization"
   - AND: Not triggering simple or extended criteria
   - Rationale: Well-understood patterns fit standard workflow

---

## Detection Function

```bash
# Three-tier detection function (POSIX-compatible)
# Analyzes requirements handoff to determine project tier
detect_tier() {
  local REQUIREMENTS_FILE="$1"  # Path to phase1-requirements-handoff.yaml
  local TIER="STANDARD"
  local CONFIDENCE="low"
  local REASON=""

  # Fail-safe: Check requirements file exists
  if [ ! -f "$REQUIREMENTS_FILE" ] || [ ! -s "$REQUIREMENTS_FILE" ]; then
    echo "STANDARD|error|Requirements file unreadable, defaulting to STANDARD"
    return
  fi

  # Extract metrics from requirements handoff
  # Assumes yq is available for YAML parsing
  if ! command -v yq &> /dev/null; then
    echo "STANDARD|error|yq not found, defaulting to STANDARD"
    return
  fi

  # Count components (from problem statement or architecture implications)
  COMPONENT_COUNT=0
  if grep -qi "component" "$REQUIREMENTS_FILE"; then
    COMPONENT_COUNT=$(grep -oi "component" "$REQUIREMENTS_FILE" | wc -l | tr -d ' ')
  fi

  # Detect statistical methods requirement
  REQUIRES_STATS=false
  if grep -qiE "statistic|hypothesis test|regression|ANOVA|chi-square|t-test|bayesian" "$REQUIREMENTS_FILE"; then
    REQUIRES_STATS=true
  fi

  # Detect mathematical/algorithm requirement
  REQUIRES_MATH=false
  if grep -qiE "algorithm|optimization|numerical|simulation|solver|matrix|complexity" "$REQUIREMENTS_FILE"; then
    REQUIRES_MATH=true
  fi

  # Estimate task count from scope
  ESTIMATED_TASKS=0
  if yq eval '.handoff.requirements.scope.in_scope' "$REQUIREMENTS_FILE" &>/dev/null; then
    ESTIMATED_TASKS=$(yq eval '.handoff.requirements.scope.in_scope | length' "$REQUIREMENTS_FILE" 2>/dev/null || echo 0)
  fi

  # Read problem statement for architectural keywords
  PROBLEM_STATEMENT=$(yq eval '.handoff.requirements.problem_statement' "$REQUIREMENTS_FILE" 2>/dev/null || echo "")

  # === EXTENDED DETECTION ===
  if [ "$COMPONENT_COUNT" -gt 5 ]; then
    TIER="EXTENDED"
    CONFIDENCE="high"
    REASON="Many components (>5)"
  elif [ "$REQUIRES_STATS" = true ] && [ "$REQUIRES_MATH" = true ]; then
    TIER="EXTENDED"
    CONFIDENCE="high"
    REASON="Requires both statistics AND mathematics"
  elif [ "$ESTIMATED_TASKS" -gt 15 ]; then
    TIER="EXTENDED"
    CONFIDENCE="high"
    REASON="Large task count (>15)"
  elif echo "$PROBLEM_STATEMENT" | grep -qiE "distributed system|microservices|event-driven|real-time processing"; then
    TIER="EXTENDED"
    CONFIDENCE="high"
    REASON="Architectural complexity keywords detected"
  elif grep -qi "extended analysis\|comprehensive review" "$REQUIREMENTS_FILE"; then
    TIER="EXTENDED"
    CONFIDENCE="high"
    REASON="User explicitly requested extended mode"
  fi

  # === SIMPLE DETECTION ===
  if [ "$CONFIDENCE" != "high" ]; then
    if [ "$COMPONENT_COUNT" -eq 1 ] && [ "$REQUIRES_STATS" = false ] && [ "$REQUIRES_MATH" = false ]; then
      TIER="SIMPLE"
      CONFIDENCE="high"
      REASON="Single component, no specialization"
    elif echo "$PROBLEM_STATEMENT" | grep -qiE "utility|script|helper|tool" && [ "$ESTIMATED_TASKS" -lt 5 ]; then
      TIER="SIMPLE"
      CONFIDENCE="high"
      REASON="Utility script or tool"
    elif echo "$PROBLEM_STATEMENT" | grep -qiE "ETL|data pipeline|batch processing" && [ "$COMPONENT_COUNT" -eq 1 ]; then
      TIER="SIMPLE"
      CONFIDENCE="high"
      REASON="Data pipeline with clear structure"
    fi
  fi

  # === STANDARD DETECTION (default) ===
  if [ "$CONFIDENCE" != "high" ]; then
    if [ "$COMPONENT_COUNT" -ge 2 ] && [ "$COMPONENT_COUNT" -le 5 ]; then
      TIER="STANDARD"
      CONFIDENCE="high"
      REASON="Multiple components (2-5)"
    elif [ "$REQUIRES_STATS" = true ] || [ "$REQUIRES_MATH" = true ]; then
      if [ "$REQUIRES_STATS" != "$REQUIRES_MATH" ]; then
        TIER="STANDARD"
        CONFIDENCE="high"
        REASON="Single specialization (statistics OR mathematics)"
      fi
    elif [ "$ESTIMATED_TASKS" -ge 5 ] && [ "$ESTIMATED_TASKS" -le 15 ]; then
      TIER="STANDARD"
      CONFIDENCE="high"
      REASON="Moderate task count (5-15)"
    elif echo "$PROBLEM_STATEMENT" | grep -qiE "web API|CLI tool|data analysis|visualization"; then
      TIER="STANDARD"
      CONFIDENCE="medium"
      REASON="Standard pattern detected"
    fi
  fi

  # === DEFAULT for unclear ===
  if [ "$CONFIDENCE" = "low" ]; then
    TIER="STANDARD"
    CONFIDENCE="low"
    REASON="Unclear project scope, defaulting to STANDARD (safest option)"
  fi

  echo "$TIER|$CONFIDENCE|$REASON"
}
```

---

## User Confirmation Logic

```bash
# Mode selection with user override capability
select_mode() {
  local REQUIREMENTS_FILE="$1"
  local SESSION_DIR="$2"

  DETECTION_RESULT=$(detect_tier "$REQUIREMENTS_FILE")
  DETECTED_TIER=$(echo "$DETECTION_RESULT" | cut -d'|' -f1)
  CONFIDENCE=$(echo "$DETECTION_RESULT" | cut -d'|' -f2)
  REASON=$(echo "$DETECTION_RESULT" | cut -d'|' -f3)

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

  # High-confidence: proceed with detected tier but allow override
  if [ "$CONFIDENCE" = "high" ]; then
    echo "High confidence in detected tier."
    echo ""
    read -t 60 -p "Proceed with $DETECTED_TIER mode? [Y/n]: " USER_CHOICE

    if [ $? -ne 0 ]; then
      # Timeout (60s)
      echo ""
      echo "No response (timeout 60s). Proceeding with detected tier: $DETECTED_TIER"
      SELECTED_TIER="$DETECTED_TIER"
    elif [ "$USER_CHOICE" = "n" ] || [ "$USER_CHOICE" = "N" ]; then
      # User wants override
      echo ""
      read -p "Select mode (1=SIMPLE, 2=STANDARD, 3=EXTENDED): " MODE_OVERRIDE
      case "$MODE_OVERRIDE" in
        1)
          SELECTED_TIER="SIMPLE"
          ;;
        2)
          SELECTED_TIER="STANDARD"
          ;;
        3)
          SELECTED_TIER="EXTENDED"
          ;;
        *)
          echo "Invalid choice. Proceeding with detected tier: $DETECTED_TIER"
          SELECTED_TIER="$DETECTED_TIER"
          ;;
      esac

      # Risky override confirmation (e.g., SIMPLE when STANDARD/EXTENDED recommended)
      if [ "$DETECTED_TIER" != "SIMPLE" ] && [ "$SELECTED_TIER" = "SIMPLE" ]; then
        echo ""
        echo "⚠️  WARNING: Selecting SIMPLE mode when $DETECTED_TIER was recommended."
        echo "This may result in insufficient analysis and coordination."
        echo ""
        read -p "Confirm risky override? [y/N]: " RISKY_CONFIRM
        if [ "$RISKY_CONFIRM" != "y" ] && [ "$RISKY_CONFIRM" != "Y" ]; then
          echo "Override cancelled. Proceeding with detected tier: $DETECTED_TIER"
          SELECTED_TIER="$DETECTED_TIER"
        else
          echo "Risky override confirmed. Proceeding with SIMPLE mode (documented)."
        fi
      fi
    else
      # User accepted detected tier
      SELECTED_TIER="$DETECTED_TIER"
    fi
  else
    # Medium/low confidence: require user confirmation
    echo "Heuristic confidence is $CONFIDENCE. User confirmation required."
    echo ""
    read -p "Select mode (1=SIMPLE, 2=STANDARD, 3=EXTENDED) [default: $DETECTED_TIER]: " USER_CHOICE

    case "$USER_CHOICE" in
      1)
        SELECTED_TIER="SIMPLE"
        ;;
      2)
        SELECTED_TIER="STANDARD"
        ;;
      3)
        SELECTED_TIER="EXTENDED"
        ;;
      "")
        SELECTED_TIER="$DETECTED_TIER"
        echo "Using default: $DETECTED_TIER"
        ;;
      *)
        echo "Invalid choice. Defaulting to STANDARD (safest option)."
        SELECTED_TIER="STANDARD"
        ;;
    esac
  fi

  # Record mode selection
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

  # Update session state
  if [ -f "$SESSION_DIR/session-state.json" ]; then
    # Use jq if available, otherwise append
    if command -v jq &> /dev/null; then
      jq --arg tier "$SELECTED_TIER" '.mode = $tier' "$SESSION_DIR/session-state.json" > "$SESSION_DIR/session-state.json.tmp"
      mv "$SESSION_DIR/session-state.json.tmp" "$SESSION_DIR/session-state.json"
    else
      # Fallback: create minimal state file
      echo "{\"mode\": \"$SELECTED_TIER\"}" > "$SESSION_DIR/session-state.json"
    fi
  fi

  echo ""
  echo "✅ Mode selected: $SELECTED_TIER"
  echo ""
  echo "================================================"

  # Export for use in workflow
  export PROGRAMMING_PM_MODE="$SELECTED_TIER"
}
```

---

## Test Cases

### Test Case 1: Simple - Single Component Utility
**Input**:
- Problem statement: "Create a utility script to parse CSV files"
- Components: 1
- Requires stats: No
- Requires math: No
- Estimated tasks: 3

**Expected**: SIMPLE (high confidence)

---

### Test Case 2: Standard - Multi-Component Web API
**Input**:
- Problem statement: "Build a REST API with database integration"
- Components: 3 (API layer, database layer, auth layer)
- Requires stats: No
- Requires math: No
- Estimated tasks: 8

**Expected**: STANDARD (high confidence)

---

### Test Case 3: Extended - Statistical Analysis + Algorithms
**Input**:
- Problem statement: "Develop statistical analysis pipeline with optimization algorithms"
- Components: 4
- Requires stats: Yes (hypothesis testing, regression)
- Requires math: Yes (optimization, numerical methods)
- Estimated tasks: 12

**Expected**: EXTENDED (high confidence, dual specialization)

---

### Test Case 4: Extended - Many Components
**Input**:
- Problem statement: "Build microservices architecture with 7 services"
- Components: 7
- Requires stats: No
- Requires math: No
- Estimated tasks: 20

**Expected**: EXTENDED (high confidence, >5 components)

---

### Test Case 5: Standard - Statistics Only
**Input**:
- Problem statement: "Implement A/B test analysis framework"
- Components: 2
- Requires stats: Yes (t-tests, chi-square)
- Requires math: No
- Estimated tasks: 6

**Expected**: STANDARD (high confidence, single specialization)

---

### Test Case 6: Simple - ETL Pipeline
**Input**:
- Problem statement: "Create ETL pipeline for data warehouse ingestion"
- Components: 1
- Requires stats: No
- Requires math: No
- Estimated tasks: 4

**Expected**: SIMPLE (high confidence, data pipeline pattern)

---

### Test Case 7: Extended - Distributed System
**Input**:
- Problem statement: "Build real-time event-driven processing system"
- Components: 4
- Requires stats: No
- Requires math: No
- Estimated tasks: 10

**Expected**: EXTENDED (high confidence, architectural complexity)

---

### Test Case 8: Standard - Data Visualization
**Input**:
- Problem statement: "Create interactive data visualization dashboard"
- Components: 3 (frontend, backend, data layer)
- Requires stats: No
- Requires math: No
- Estimated tasks: 9

**Expected**: STANDARD (medium confidence, standard pattern)

---

### Test Case 9: Standard (Default) - Unclear Scope
**Input**:
- Problem statement: "Improve data processing workflow"
- Components: Unknown
- Requires stats: Maybe
- Requires math: Maybe
- Estimated tasks: Unknown

**Expected**: STANDARD (low confidence, default to safest)

---

### Test Case 10: Extended - User Request
**Input**:
- Problem statement: "Build CLI tool with extended analysis and comprehensive review"
- Components: 2
- Requires stats: No
- Requires math: No
- Estimated tasks: 7

**Expected**: EXTENDED (high confidence, user explicit request)

---

## Mode-Based Workflow Branching

After mode selection, the workflow branches:

### SIMPLE Mode
- **Skip**: Phase 2.5 (extended risk analysis)
- **Specialist invocation**: Sequential (not parallel)
- **Quality gates**: Automated checks only (skip manual reviews where possible)
- **Estimated duration**: 1-2 hours

### STANDARD Mode (Default)
- **Include**: All standard phases (0-6)
- **Specialist invocation**: Wave-based parallel execution
- **Quality gates**: Full automated + manual checks
- **Estimated duration**: 4-6 hours

### EXTENDED Mode
- **Include**: All phases + Phase 2.5 (extended pre-mortem)
- **Specialist invocation**: Wave-based parallel with extended timeouts
- **Quality gates**: Full checks + additional architectural review
- **Code review**: senior-developer reviews ALL code (including senior-developer outputs)
- **Estimated duration**: 8-12 hours

---

## Integration with Session State

Mode selection updates `session-state.json`:

```json
{
  "session_dir": "/tmp/programming-pm-session-...",
  "mode": "STANDARD",
  "phase": 1,
  "status": "active",
  "start_time": "2025-02-05T10:00:00Z",
  "mode_selection": {
    "detected": "STANDARD",
    "confidence": "high",
    "selected": "STANDARD",
    "override": false,
    "timestamp": "2025-02-05T10:15:00Z"
  }
}
```

Workflow logic checks mode before executing phases:

```bash
MODE=$(jq -r '.mode' "$SESSION_DIR/session-state.json")

if [ "$MODE" = "SIMPLE" ]; then
  echo "SIMPLE mode: Skipping Phase 2.5"
  # Jump to Phase 3
elif [ "$MODE" = "EXTENDED" ]; then
  echo "EXTENDED mode: Including Phase 2.5 with extended analysis"
  # Run all phases including 2.5
else
  echo "STANDARD mode: Full workflow"
  # Run all standard phases
fi
```

---

## Backwards Compatibility

For existing sessions without mode selection:

```bash
# Check if mode-selection.json exists
if [ ! -f "$SESSION_DIR/mode-selection.json" ]; then
  echo "⚠️  Legacy session detected (no mode selection)"
  echo "Defaulting to STANDARD mode for safety"
  export PROGRAMMING_PM_MODE="STANDARD"
else
  export PROGRAMMING_PM_MODE=$(jq -r '.selected_tier' "$SESSION_DIR/mode-selection.json")
fi
```

---

## Metrics and Continuous Improvement

Track mode selection accuracy:

```bash
# Log mode selection for retrospective
cat >> "$SESSION_DIR/mode-selection-log.txt" <<EOF
$(date -u +"%Y-%m-%dT%H:%M:%SZ"),DETECTED:$DETECTED_TIER,SELECTED:$SELECTED_TIER,CONFIDENCE:$CONFIDENCE,REASON:$REASON
EOF
```

Periodically review logs to:
- Identify patterns in user overrides
- Adjust detection criteria if override rate >20%
- Improve confidence calibration
