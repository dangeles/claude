# Complexity Detection Criteria

## Purpose

Determine when to invoke Phase 2.5 (strategy-consultant agent) vs. skip directly to Phase 3.

**Goal**: Avoid overhead for simple changes while ensuring strategic review for complex changes.

**Target metrics**:
- False positive rate: <10% (simple changes incorrectly flagged as complex)
- False negative rate: <5% (complex changes incorrectly flagged as simple)

---

## Detection Criteria

### High-Confidence Complex (Always Triggers Phase 2.5)

1. **New skill creation**
   - Specification contains: "Create new skill"
   - Rationale: New skills are architectural additions requiring strategic assessment

2. **Many files affected** (>3 files)
   - File count exceeds 3
   - Rationale: Multi-file changes likely have architectural implications

3. **Large change** (>200 lines)
   - Total lines changed exceeds 200
   - Rationale: Large changes carry higher risk of architectural issues

4. **Explicit user request**
   - Specification contains: "strategic review" or "architectural assessment"
   - Rationale: User explicitly wants strategic perspective

---

### High-Confidence Simple (Always Skips Phase 2.5)

1. **Documentation-only changes**
   - Scope mentions: "documentation" or "typo" or "comment" or "example"
   - AND: Single file, ≤50 lines
   - Rationale: Documentation has no architectural implications

2. **Minor bug fix**
   - Scope mentions: "fix bug" or "fix typo" or "fix error"
   - AND: Single file, ≤50 lines
   - Rationale: Small fixes don't require strategic assessment

---

### Medium-Confidence (User Confirmation Required)

1. **Keywords detected but small change**
   - Contains: "agent" or "workflow" or "phase" or "quality gate"
   - BUT: ≤2 files AND ≤100 lines
   - Action: Prompt user with reasoning, default to skip

2. **Moderate size with unclear scope**
   - 2-3 files OR 100-200 lines
   - No clear documentation/bug-fix/architectural indicators
   - Action: Prompt user with reasoning, default to complex

---

## Detection Function

```bash
detect_complexity() {
  local SPEC_FILE="/tmp/skill-editor-session/refined-specification.md"
  local COMPLEX=false
  local CONFIDENCE="low"
  local REASON=""

  # Extract metrics from spec
  FILES_CHANGED=$(grep -c "File:" "$SPEC_FILE" 2>/dev/null || echo 0)
  LINES_CHANGED=$(grep -oP "Lines: \K[0-9]+" "$SPEC_FILE" 2>/dev/null | awk '{sum+=$1} END {print sum}')
  [ -z "$LINES_CHANGED" ] && LINES_CHANGED=0

  SCOPE=$(grep -A10 "^## Scope" "$SPEC_FILE")

  # High-confidence complex triggers
  if grep -qi "Create new skill" "$SPEC_FILE"; then
    COMPLEX=true
    CONFIDENCE="high"
    REASON="New skill creation"
  elif [ "$FILES_CHANGED" -gt 3 ]; then
    COMPLEX=true
    CONFIDENCE="high"
    REASON="Multiple files affected (>3)"
  elif [ "$LINES_CHANGED" -gt 200 ]; then
    COMPLEX=true
    CONFIDENCE="high"
    REASON="Large change (>200 lines)"
  elif grep -qi "strategic review\|architectural assessment" "$SPEC_FILE"; then
    COMPLEX=true
    CONFIDENCE="high"
    REASON="User explicitly requested strategic review"
  fi

  # High-confidence simple (override complex if both match)
  if [ "$CONFIDENCE" != "high" ]; then
    if echo "$SCOPE" | grep -qi "documentation\|typo\|comment\|example"; then
      if [ "$FILES_CHANGED" -le 1 ] && [ "$LINES_CHANGED" -le 50 ]; then
        COMPLEX=false
        CONFIDENCE="high"
        REASON="Documentation-only change"
      fi
    fi

    if echo "$SCOPE" | grep -qi "fix bug\|fix typo\|fix error"; then
      if [ "$FILES_CHANGED" -le 1 ] && [ "$LINES_CHANGED" -le 50 ]; then
        COMPLEX=false
        CONFIDENCE="high"
        REASON="Minor bug fix"
      fi
    fi
  fi

  # Medium-confidence detection
  if [ "$CONFIDENCE" != "high" ]; then
    if grep -qi "agent\|workflow\|phase\|quality gate\|multi-agent" "$SPEC_FILE"; then
      if [ "$FILES_CHANGED" -le 2 ] && [ "$LINES_CHANGED" -le 100 ]; then
        COMPLEX=false
        CONFIDENCE="medium"
        REASON="Keywords detected but change is small (user confirmation recommended)"
      else
        COMPLEX=true
        CONFIDENCE="medium"
        REASON="Workflow/agent keywords with moderate change size"
      fi
    fi
  fi

  # Default for unclear cases
  if [ "$CONFIDENCE" = "low" ]; then
    if [ "$FILES_CHANGED" -ge 2 ] || [ "$LINES_CHANGED" -ge 100 ]; then
      COMPLEX=true
      CONFIDENCE="low"
      REASON="Moderate size with unclear scope (user confirmation recommended)"
    else
      COMPLEX=false
      CONFIDENCE="medium"
      REASON="Small change with unclear scope"
    fi
  fi

  echo "$COMPLEX|$CONFIDENCE|$REASON"
}
```

---

## User Confirmation Logic

```bash
DETECTION_RESULT=$(detect_complexity)
IS_COMPLEX=$(echo "$DETECTION_RESULT" | cut -d'|' -f1)
CONFIDENCE=$(echo "$DETECTION_RESULT" | cut -d'|' -f2)
REASON=$(echo "$DETECTION_RESULT" | cut -d'|' -f3)

echo "Complexity detection: $IS_COMPLEX (confidence: $CONFIDENCE)"
echo "Reason: $REASON"
echo ""

# High-confidence: proceed automatically
if [ "$CONFIDENCE" = "high" ]; then
  if [ "$IS_COMPLEX" = "true" ]; then
    echo "→ Complex change detected: Launching Phase 2.5 (strategy consultant)"
  else
    echo "→ Simple change detected: Skipping Phase 2.5 (proceeding directly to Phase 3)"
  fi
else
  # Medium/low confidence: ask user
  echo "Heuristic confidence is $CONFIDENCE. User confirmation recommended."
  echo ""
  read -p "Do you want strategic review (Phase 2.5)?
  (Y) Yes - run Phase 2.5 (adds 10-30 min, strategic architectural assessment)
  (N) No - skip Phase 2.5 (faster, proceed directly to synthesis)
Choice [Y/n]: " USER_OVERRIDE

  if [ "$USER_OVERRIDE" = "n" ] || [ "$USER_OVERRIDE" = "N" ]; then
    IS_COMPLEX=false
    echo "→ User override: Skipping Phase 2.5"
  else
    IS_COMPLEX=true
    echo "→ User confirmed: Running Phase 2.5"
  fi
fi

# Record decision in session state
jq -n \
  --argjson complex "$IS_COMPLEX" \
  --arg confidence "$CONFIDENCE" \
  --arg reason "$REASON" \
  '{
    complexity_detected: $complex,
    confidence: $confidence,
    reason: $reason,
    timestamp: (now | strftime("%Y-%m-%dT%H:%M:%SZ"))
  }' \
  > /tmp/skill-editor-session/complexity-detection.json
```

---

## Test Cases

### Test Case 1: Simple Change (Should Skip)

**Specification**:
```markdown
## Objective
Fix typo in skill-editor SKILL.md documentation

## Scope
- Edit 1 file: claude-config/skills/skill-editor/SKILL.md
- Lines changed: 1 line
- Change: "strategig" → "strategic"
```

**Expected Result**:
- IS_COMPLEX=false
- CONFIDENCE=high
- REASON="Documentation-only change"
- Phase 2.5: SKIP

---

### Test Case 2: New Skill Creation (Should Trigger)

**Specification**:
```markdown
## Objective
Create new skill for code review automation

## Scope
- Create new skill: claude-config/skills/code-reviewer/SKILL.md
- Create agent: claude-config/agents/code-reviewer.md
- Create references: claude-config/skills/code-reviewer/references/
```

**Expected Result**:
- IS_COMPLEX=true
- CONFIDENCE=high
- REASON="New skill creation"
- Phase 2.5: TRIGGER

---

### Test Case 3: Major Refactoring (Should Trigger)

**Specification**:
```markdown
## Objective
Refactor skill-editor to use hierarchical agent coordination

## Scope
- Modify 4 files: SKILL.md, 3 agent files
- Lines changed: ~250 lines
- Add new workflow phases
```

**Expected Result**:
- IS_COMPLEX=true
- CONFIDENCE=high
- REASON="Multiple files affected (>3)" OR "Large change (>200 lines)"
- Phase 2.5: TRIGGER

---

### Test Case 4: Borderline (User Confirmation)

**Specification**:
```markdown
## Objective
Update skill-editor agent coordination pattern

## Scope
- Modify 2 files: SKILL.md, decision-synthesizer.md
- Lines changed: ~80 lines
- Change workflow description
```

**Expected Result**:
- IS_COMPLEX=false
- CONFIDENCE=medium
- REASON="Keywords detected but change is small (user confirmation recommended)"
- Action: PROMPT USER

---

### Test Case 5: Documentation Update (Should Skip)

**Specification**:
```markdown
## Objective
Update skill-editor README with usage examples

## Scope
- Modify 1 file: README.md
- Lines changed: 30 lines
- Add 2 examples
```

**Expected Result**:
- IS_COMPLEX=false
- CONFIDENCE=high
- REASON="Documentation-only change"
- Phase 2.5: SKIP

---

## Monitoring and Refinement

**Track metrics over time**:
- False positive rate: (Simple changes that triggered Phase 2.5) / (Total simple changes)
- False negative rate: (Complex changes that skipped Phase 2.5) / (Total complex changes)
- User override rate: (User overrides) / (Total borderline cases)

**Refinement process**:
1. After first 10 uses of Phase 2.5, review false positive/negative incidents
2. Adjust thresholds if false positive rate >10% or false negative rate >5%
3. Consider adding new high-confidence indicators based on patterns
4. Document adjustments in this file

**Target thresholds for adjustment**:
- If false positive >10%: Increase file count threshold (3 → 4) or line threshold (200 → 300)
- If false negative >5%: Decrease thresholds or add new keyword indicators
