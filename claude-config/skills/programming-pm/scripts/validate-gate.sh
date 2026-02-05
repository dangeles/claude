#!/usr/bin/env bash
#
# Quality Gate Validation Script for programming-pm workflow
#
# Validates each of 6 quality gates with automated checks.
# Called at each phase boundary to enforce quality criteria.
#
# Usage:
#   ./validate-gate.sh <gate_number> <handoff_path> [session_dir]
#
# Example:
#   ./validate-gate.sh 1 /tmp/session/handoffs/phase1-requirements-handoff.yaml /tmp/session
#
# Exit codes:
#   0: Gate PASSED
#   1: Gate FAILED
#   2: Invalid usage
#

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <gate_number> <handoff_path> [session_dir]" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  $0 1 /tmp/session/handoffs/phase1-requirements-handoff.yaml" >&2
    echo "  $0 6 /tmp/session/handoffs/phase5-review-handoff.yaml /tmp/session" >&2
    exit 2
fi

GATE_NUMBER="$1"
HANDOFF_PATH="$2"
SESSION_DIR="${3:-}"

# Validate gate number
if ! [[ "$GATE_NUMBER" =~ ^[1-6]$ ]]; then
    echo -e "${RED}❌ Invalid gate number: $GATE_NUMBER (must be 1-6)${NC}" >&2
    exit 2
fi

# Check if handoff file exists
if [ ! -f "$HANDOFF_PATH" ]; then
    echo -e "${RED}❌ Handoff file not found: $HANDOFF_PATH${NC}" >&2
    exit 1
fi

# Helper functions
check_yaml_field() {
    local file="$1"
    local field="$2"
    local expected_value="${3:-}"

    if command -v yq &> /dev/null; then
        if [ -n "$expected_value" ]; then
            result=$(yq eval "$field" "$file" 2>/dev/null || echo "null")
            [ "$result" = "$expected_value" ]
        else
            result=$(yq eval "$field" "$file" 2>/dev/null || echo "null")
            [ "$result" != "null" ] && [ "$result" != "" ]
        fi
    else
        # Fallback to grep if yq not available
        grep -q "$field" "$file"
    fi
}

count_yaml_array() {
    local file="$1"
    local field="$2"

    if command -v yq &> /dev/null; then
        yq eval "$field | length" "$file" 2>/dev/null || echo "0"
    else
        # Fallback: count lines indented under field
        grep -A 100 "$field" "$file" | grep -c "^  - " || echo "0"
    fi
}

# Gate 1: Requirements Approval
validate_gate_1() {
    echo -e "${YELLOW}Validating Quality Gate 1: Requirements Approval${NC}"

    local errors=0

    # Check 1: Problem statement exists and is specific
    if check_yaml_field "$HANDOFF_PATH" "handoff.requirements.problem_statement"; then
        problem_statement=$(yq '.handoff.requirements.problem_statement' "$HANDOFF_PATH" 2>/dev/null || echo "")

        # Check for vague terms
        if echo "$problem_statement" | grep -qiE '\b(better|faster|improve|enhance)\b' && ! echo "$problem_statement" | grep -qE '[0-9]+'; then
            echo -e "${RED}  ✗ Problem statement contains vague terms without measurable targets${NC}"
            errors=$((errors + 1))
        else
            echo -e "${GREEN}  ✓ Problem statement is specific${NC}"
        fi
    else
        echo -e "${RED}  ✗ Missing problem statement${NC}"
        errors=$((errors + 1))
    fi

    # Check 2: Success criteria are measurable
    success_criteria_count=$(count_yaml_array "$HANDOFF_PATH" "handoff.requirements.success_criteria")
    if [ "$success_criteria_count" -gt 0 ]; then
        echo -e "${GREEN}  ✓ Success criteria defined ($success_criteria_count items)${NC}"
    else
        echo -e "${RED}  ✗ No success criteria defined${NC}"
        errors=$((errors + 1))
    fi

    # Check 3: Scope boundaries (IN/OUT) explicitly defined
    if check_yaml_field "$HANDOFF_PATH" "handoff.requirements.scope.in_scope" && \
       check_yaml_field "$HANDOFF_PATH" "handoff.requirements.scope.out_of_scope"; then
        echo -e "${GREEN}  ✓ Scope boundaries defined (in_scope and out_of_scope)${NC}"
    else
        echo -e "${RED}  ✗ Missing scope boundaries${NC}"
        errors=$((errors + 1))
    fi

    # Check 4: Dependencies identified
    if check_yaml_field "$HANDOFF_PATH" "handoff.requirements.dependencies"; then
        deps_count=$(count_yaml_array "$HANDOFF_PATH" "handoff.requirements.dependencies")
        echo -e "${GREEN}  ✓ Dependencies field present ($deps_count items)${NC}"
    else
        echo -e "${RED}  ✗ Missing dependencies field${NC}"
        errors=$((errors + 1))
    fi

    return $errors
}

# Gate 2: Pre-Mortem Completion
validate_gate_2() {
    echo -e "${YELLOW}Validating Quality Gate 2: Pre-Mortem Completion${NC}"

    local errors=0

    # Check 1: At least 3 risks identified
    risk_count=$(count_yaml_array "$HANDOFF_PATH" "handoff.risk_summary")
    if [ "$risk_count" -ge 3 ]; then
        echo -e "${GREEN}  ✓ At least 3 risks identified ($risk_count total)${NC}"
    else
        echo -e "${RED}  ✗ Insufficient risks identified ($risk_count, minimum 3)${NC}"
        errors=$((errors + 1))
    fi

    # Check 2: Each risk has likelihood and impact ratings
    # This requires iterating through risk_summary array
    if command -v yq &> /dev/null; then
        for i in $(seq 0 $((risk_count - 1))); do
            if ! check_yaml_field "$HANDOFF_PATH" "handoff.risk_summary[$i].score"; then
                echo -e "${RED}  ✗ Risk $i missing score${NC}"
                errors=$((errors + 1))
            fi
        done
        echo -e "${GREEN}  ✓ All risks have scores${NC}"
    else
        echo -e "${YELLOW}  ⚠ Cannot validate risk scores (yq not available)${NC}"
    fi

    # Check 3: Each risk has disposition
    if command -v yq &> /dev/null; then
        for i in $(seq 0 $((risk_count - 1))); do
            if ! check_yaml_field "$HANDOFF_PATH" "handoff.risk_summary[$i].disposition"; then
                echo -e "${RED}  ✗ Risk $i missing disposition${NC}"
                errors=$((errors + 1))
            fi
        done
        echo -e "${GREEN}  ✓ All risks have disposition${NC}"
    else
        echo -e "${YELLOW}  ⚠ Cannot validate dispositions (yq not available)${NC}"
    fi

    # Check 4: Critical risks (score >= 15) have mitigation plans
    if command -v yq &> /dev/null; then
        critical_count=0
        mitigated_count=0

        for i in $(seq 0 $((risk_count - 1))); do
            score=$(yq eval ".handoff.risk_summary[$i].score" "$HANDOFF_PATH" 2>/dev/null || echo "0")
            disposition=$(yq eval ".handoff.risk_summary[$i].disposition" "$HANDOFF_PATH" 2>/dev/null || echo "")

            if [ "$score" -ge 15 ]; then
                critical_count=$((critical_count + 1))
                if [ "$disposition" = "mitigate" ] && check_yaml_field "$HANDOFF_PATH" "handoff.risk_summary[$i].mitigation"; then
                    mitigated_count=$((mitigated_count + 1))
                fi
            fi
        done

        if [ "$critical_count" -eq 0 ]; then
            echo -e "${GREEN}  ✓ No critical risks (score >= 15)${NC}"
        elif [ "$mitigated_count" -eq "$critical_count" ]; then
            echo -e "${GREEN}  ✓ All critical risks have mitigation plans${NC}"
        else
            echo -e "${RED}  ✗ Critical risks without mitigation: $((critical_count - mitigated_count))${NC}"
            errors=$((errors + 1))
        fi
    fi

    return $errors
}

# Gate 3: Architecture Approval
validate_gate_3() {
    echo -e "${YELLOW}Validating Quality Gate 3: Architecture Approval${NC}"

    local errors=0

    # Check 1: All components identified with responsibilities
    component_count=$(count_yaml_array "$HANDOFF_PATH" "handoff.components")
    if [ "$component_count" -gt 0 ]; then
        echo -e "${GREEN}  ✓ Components identified ($component_count total)${NC}"

        # Verify each component has required fields
        if command -v yq &> /dev/null; then
            for i in $(seq 0 $((component_count - 1))); do
                if ! check_yaml_field "$HANDOFF_PATH" "handoff.components[$i].responsibility"; then
                    echo -e "${RED}  ✗ Component $i missing responsibility${NC}"
                    errors=$((errors + 1))
                fi
            done
        fi
    else
        echo -e "${RED}  ✗ No components identified${NC}"
        errors=$((errors + 1))
    fi

    # Check 2: Data flow documented
    if check_yaml_field "$HANDOFF_PATH" "handoff.data_flow.description"; then
        echo -e "${GREEN}  ✓ Data flow documented${NC}"
    else
        echo -e "${RED}  ✗ Missing data flow description${NC}"
        errors=$((errors + 1))
    fi

    # Check 3: Technology choices justified
    tech_count=$(count_yaml_array "$HANDOFF_PATH" "handoff.technology_choices")
    if [ "$tech_count" -gt 0 ]; then
        echo -e "${GREEN}  ✓ Technology choices documented ($tech_count items)${NC}"
    else
        echo -e "${RED}  ✗ No technology choices documented${NC}"
        errors=$((errors + 1))
    fi

    # Check 4: Component interfaces defined
    if command -v yq &> /dev/null; then
        interfaces_ok=true
        for i in $(seq 0 $((component_count - 1))); do
            if ! check_yaml_field "$HANDOFF_PATH" "handoff.components[$i].interfaces.inputs"; then
                echo -e "${RED}  ✗ Component $i missing input interfaces${NC}"
                interfaces_ok=false
                errors=$((errors + 1))
            fi
            if ! check_yaml_field "$HANDOFF_PATH" "handoff.components[$i].interfaces.outputs"; then
                echo -e "${RED}  ✗ Component $i missing output interfaces${NC}"
                interfaces_ok=false
                errors=$((errors + 1))
            fi
        done
        if $interfaces_ok; then
            echo -e "${GREEN}  ✓ Component interfaces defined${NC}"
        fi
    fi

    # Check 5: Testing strategy outlined
    if check_yaml_field "$HANDOFF_PATH" "handoff.testing_strategy.unit" && \
       check_yaml_field "$HANDOFF_PATH" "handoff.testing_strategy.integration"; then
        echo -e "${GREEN}  ✓ Testing strategy outlined${NC}"
    else
        echo -e "${RED}  ✗ Incomplete testing strategy${NC}"
        errors=$((errors + 1))
    fi

    return $errors
}

# Gate 4: Implementation Validation
validate_gate_4() {
    echo -e "${YELLOW}Validating Quality Gate 4: Implementation Validation${NC}"

    local errors=0

    # Check 1: Specialist outputs exist
    if [ -z "$SESSION_DIR" ]; then
        echo -e "${YELLOW}  ⚠ No session directory provided, skipping file checks${NC}"
        return 0
    fi

    # Check for code handoff files
    code_handoffs=$(find "$SESSION_DIR/handoffs" -name "phase4-code-handoff-*.yaml" 2>/dev/null | wc -l)
    if [ "$code_handoffs" -gt 0 ]; then
        echo -e "${GREEN}  ✓ Code handoffs present ($code_handoffs files)${NC}"
    else
        echo -e "${RED}  ✗ No code handoffs found${NC}"
        errors=$((errors + 1))
    fi

    # Check 2: Outputs meet minimum length criteria (>100 words)
    if command -v yq &> /dev/null && [ -f "$HANDOFF_PATH" ]; then
        if check_yaml_field "$HANDOFF_PATH" "handoff.summary"; then
            summary=$(yq '.handoff.summary' "$HANDOFF_PATH" 2>/dev/null || echo "")
            word_count=$(echo "$summary" | wc -w | tr -d ' ')

            if [ "$word_count" -ge 100 ]; then
                echo -e "${GREEN}  ✓ Summary meets minimum length ($word_count words)${NC}"
            else
                echo -e "${RED}  ✗ Summary too short ($word_count words, minimum 100)${NC}"
                errors=$((errors + 1))
            fi
        fi
    fi

    # Check 3: No critical blocking issues
    if check_yaml_field "$HANDOFF_PATH" "handoff.quality.status"; then
        status=$(yq '.handoff.quality.status' "$HANDOFF_PATH" 2>/dev/null || echo "")
        if [ "$status" = "complete" ] || [ "$status" = "partial" ]; then
            echo -e "${GREEN}  ✓ Quality status acceptable ($status)${NC}"
        else
            echo -e "${RED}  ✗ Invalid quality status: $status${NC}"
            errors=$((errors + 1))
        fi
    fi

    # Check 4: Acceptance criteria met
    if check_yaml_field "$HANDOFF_PATH" "handoff.context.task_id"; then
        echo -e "${GREEN}  ✓ Task properly identified${NC}"
    else
        echo -e "${RED}  ✗ Missing task identification${NC}"
        errors=$((errors + 1))
    fi

    return $errors
}

# Gate 5: Code Review Approval
validate_gate_5() {
    echo -e "${YELLOW}Validating Quality Gate 5: Code Review Approval${NC}"

    local errors=0

    # Check automated checks
    if check_yaml_field "$HANDOFF_PATH" "handoff.automated_checks"; then
        # Check ruff
        if check_yaml_field "$HANDOFF_PATH" "handoff.automated_checks.ruff" "pass"; then
            echo -e "${GREEN}  ✓ Ruff checks passed${NC}"
        else
            echo -e "${RED}  ✗ Ruff checks failed${NC}"
            errors=$((errors + 1))
        fi

        # Check mypy
        if check_yaml_field "$HANDOFF_PATH" "handoff.automated_checks.mypy" "pass"; then
            echo -e "${GREEN}  ✓ Mypy checks passed${NC}"
        else
            echo -e "${RED}  ✗ Mypy checks failed${NC}"
            errors=$((errors + 1))
        fi

        # Check tests
        if check_yaml_field "$HANDOFF_PATH" "handoff.automated_checks.tests" "pass"; then
            echo -e "${GREEN}  ✓ Tests passed${NC}"
        else
            echo -e "${RED}  ✗ Tests failed${NC}"
            errors=$((errors + 1))
        fi

        # Check coverage
        if command -v yq &> /dev/null; then
            coverage=$(yq '.handoff.automated_checks.coverage' "$HANDOFF_PATH" 2>/dev/null || echo "0")
            if awk "BEGIN {exit !($coverage >= 80)}"; then
                echo -e "${GREEN}  ✓ Coverage meets threshold (${coverage}%)${NC}"
            else
                echo -e "${RED}  ✗ Coverage below threshold (${coverage}%, minimum 80%)${NC}"
                errors=$((errors + 1))
            fi
        fi
    else
        echo -e "${RED}  ✗ Missing automated_checks field${NC}"
        errors=$((errors + 1))
    fi

    # Check manual review
    if check_yaml_field "$HANDOFF_PATH" "handoff.manual_review"; then
        if command -v yq &> /dev/null; then
            for aspect in code_quality documentation testing architecture; do
                status=$(yq eval ".handoff.manual_review.$aspect" "$HANDOFF_PATH" 2>/dev/null || echo "")
                if [ "$status" = "pass" ]; then
                    echo -e "${GREEN}  ✓ Manual review: $aspect passed${NC}"
                else
                    echo -e "${RED}  ✗ Manual review: $aspect has issues${NC}"
                    errors=$((errors + 1))
                fi
            done
        fi
    else
        echo -e "${RED}  ✗ Missing manual_review field${NC}"
        errors=$((errors + 1))
    fi

    # Check approval
    if check_yaml_field "$HANDOFF_PATH" "handoff.approval.approved" "true"; then
        echo -e "${GREEN}  ✓ Review approved${NC}"
    else
        echo -e "${RED}  ✗ Review not approved${NC}"
        errors=$((errors + 1))
    fi

    return $errors
}

# Gate 6: PR Merge
validate_gate_6() {
    echo -e "${YELLOW}Validating Quality Gate 6: VCS Integration${NC}"

    local errors=0

    # Check 1: All previous gates passed (check handoff history)
    if check_yaml_field "$HANDOFF_PATH" "handoff.quality.status" "complete"; then
        echo -e "${GREEN}  ✓ Quality status complete${NC}"
    else
        echo -e "${YELLOW}  ⚠ Quality status not complete${NC}"
    fi

    # Check 2: No merge conflicts (check git status if in repo)
    if [ -n "$SESSION_DIR" ] && command -v git &> /dev/null; then
        # Navigate to repo root
        if git rev-parse --git-dir > /dev/null 2>&1; then
            if git diff --name-only --diff-filter=U | grep -q .; then
                echo -e "${RED}  ✗ Merge conflicts detected${NC}"
                errors=$((errors + 1))
            else
                echo -e "${GREEN}  ✓ No merge conflicts${NC}"
            fi

            # Check for unstaged changes
            if git diff --quiet && git diff --cached --quiet; then
                echo -e "${GREEN}  ✓ Working directory clean${NC}"
            else
                echo -e "${YELLOW}  ⚠ Uncommitted changes present${NC}"
            fi
        else
            echo -e "${YELLOW}  ⚠ Not in git repository${NC}"
        fi
    fi

    # Check 3: Handoff contains review approval
    if check_yaml_field "$HANDOFF_PATH" "handoff.review.status" "approved"; then
        echo -e "${GREEN}  ✓ Review approved${NC}"
    else
        echo -e "${RED}  ✗ Review not approved${NC}"
        errors=$((errors + 1))
    fi

    # Check 4: PR description requirements met
    if check_yaml_field "$HANDOFF_PATH" "handoff.deliverable.location"; then
        echo -e "${GREEN}  ✓ Deliverable location documented${NC}"
    else
        echo -e "${RED}  ✗ Missing deliverable location${NC}"
        errors=$((errors + 1))
    fi

    # Check 5: Only staged files (if in git repo)
    if [ -n "$SESSION_DIR" ] && command -v git &> /dev/null; then
        if git rev-parse --git-dir > /dev/null 2>&1; then
            staged_count=$(git diff --cached --name-only | wc -l | tr -d ' ')
            if [ "$staged_count" -gt 0 ]; then
                echo -e "${GREEN}  ✓ Files staged for commit ($staged_count files)${NC}"
            else
                echo -e "${YELLOW}  ⚠ No files staged${NC}"
            fi
        fi
    fi

    return $errors
}

# Main validation logic
main() {
    echo ""
    echo "================================================"
    echo "  Quality Gate $GATE_NUMBER Validation"
    echo "================================================"
    echo ""

    case "$GATE_NUMBER" in
        1)
            validate_gate_1
            ;;
        2)
            validate_gate_2
            ;;
        3)
            validate_gate_3
            ;;
        4)
            validate_gate_4
            ;;
        5)
            validate_gate_5
            ;;
        6)
            validate_gate_6
            ;;
        *)
            echo -e "${RED}❌ Unknown gate number: $GATE_NUMBER${NC}" >&2
            exit 2
            ;;
    esac

    result=$?

    echo ""
    echo "================================================"

    if [ $result -eq 0 ]; then
        echo -e "${GREEN}✅ Quality Gate $GATE_NUMBER: PASSED${NC}"
        echo "================================================"
        return 0
    else
        echo -e "${RED}❌ Quality Gate $GATE_NUMBER: FAILED${NC}"
        echo -e "${RED}   $result check(s) failed${NC}"
        echo "================================================"
        return 1
    fi
}

main
