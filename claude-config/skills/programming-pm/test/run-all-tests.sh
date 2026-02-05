#!/usr/bin/env bash
#
# Automated test runner for programming-pm
#

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(dirname "$SCRIPT_DIR")"

echo ""
echo "================================================"
echo "  Programming-PM Test Suite"
echo "================================================"
echo ""

TESTS_PASSED=0
TESTS_FAILED=0
TESTS_SKIPPED=0

# Test 1: Handoff Validation - Valid Session
echo -e "${BLUE}Test 1: Handoff Validation - Valid Session${NC}"
if python3 "$SKILL_DIR/scripts/validate-handoff.py" \
   "$SCRIPT_DIR/fixtures/valid-session-handoff.yaml" \
   "session_handoff" > /dev/null 2>&1; then
  echo -e "${GREEN}✅ PASSED${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
else
  echo -e "${RED}❌ FAILED${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Test 2: Handoff Validation - Invalid Session (should fail)
echo -e "${BLUE}Test 2: Handoff Validation - Invalid Session (expect failure)${NC}"
if python3 "$SKILL_DIR/scripts/validate-handoff.py" \
   "$SCRIPT_DIR/fixtures/invalid-session-handoff.yaml" \
   "session_handoff" > /dev/null 2>&1; then
  echo -e "${RED}❌ FAILED (should have rejected invalid handoff)${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
else
  echo -e "${GREEN}✅ PASSED (correctly rejected invalid handoff)${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
fi
echo ""

# Test 3: Handoff Validation - Valid Requirements
echo -e "${BLUE}Test 3: Handoff Validation - Valid Requirements${NC}"
if python3 "$SKILL_DIR/scripts/validate-handoff.py" \
   "$SCRIPT_DIR/fixtures/valid-requirements-handoff.yaml" \
   "requirements_handoff" > /dev/null 2>&1; then
  echo -e "${GREEN}✅ PASSED${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
else
  echo -e "${RED}❌ FAILED${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Test 4: Mode Detection
echo -e "${BLUE}Test 4: Mode Detection${NC}"
if bash "$SCRIPT_DIR/test-mode-detection.sh" > /dev/null 2>&1; then
  echo -e "${GREEN}✅ PASSED${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
else
  echo -e "${RED}❌ FAILED${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Test 5: Quality Gate Validation (requires yq)
echo -e "${BLUE}Test 5: Quality Gate Validation${NC}"
if command -v yq &> /dev/null; then
  if bash "$SKILL_DIR/scripts/validate-gate.sh" 1 \
     "$SCRIPT_DIR/fixtures/valid-requirements-handoff.yaml" > /dev/null 2>&1; then
    echo -e "${GREEN}✅ PASSED${NC}"
    TESTS_PASSED=$((TESTS_PASSED + 1))
  else
    echo -e "${RED}❌ FAILED${NC}"
    TESTS_FAILED=$((TESTS_FAILED + 1))
  fi
else
  echo -e "${YELLOW}⚠️  SKIPPED (yq not installed)${NC}"
  TESTS_SKIPPED=$((TESTS_SKIPPED + 1))
fi
echo ""

# Test 6: Validate script executability
echo -e "${BLUE}Test 6: Script Executability${NC}"
if [ -x "$SKILL_DIR/scripts/validate-handoff.py" ] && \
   [ -x "$SKILL_DIR/scripts/validate-gate.sh" ]; then
  echo -e "${GREEN}✅ PASSED (scripts are executable)${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
else
  echo -e "${RED}❌ FAILED (scripts not executable)${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Test 7: Check required files exist
echo -e "${BLUE}Test 7: Required Files Exist${NC}"
REQUIRED_FILES=(
  "$SKILL_DIR/scripts/validate-handoff.py"
  "$SKILL_DIR/scripts/validate-gate.sh"
  "$SKILL_DIR/references/mode-selection-criteria.md"
  "$SKILL_DIR/references/handoff-schema.md"
)

FILES_OK=true
for FILE in "${REQUIRED_FILES[@]}"; do
  if [ ! -f "$FILE" ]; then
    echo -e "${RED}  Missing: $(basename "$FILE")${NC}"
    FILES_OK=false
  fi
done

if $FILES_OK; then
  echo -e "${GREEN}✅ PASSED (all required files present)${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
else
  echo -e "${RED}❌ FAILED (missing files)${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Test 8: Check SKILL.md integration points
echo -e "${BLUE}Test 8: SKILL.md Integration Points${NC}"
INTEGRATION_CHECKS=(
  "Mode Selection"
  "Handoff Validation"
  "Wave-Based Parallel Execution"
  "sync-config.py"
  "Quality Gate 4a"
)

INTEGRATION_OK=true
for CHECK in "${INTEGRATION_CHECKS[@]}"; do
  if ! grep -q "$CHECK" "$SKILL_DIR/SKILL.md" 2>/dev/null; then
    echo -e "${RED}  Missing: $CHECK${NC}"
    INTEGRATION_OK=false
  fi
done

if $INTEGRATION_OK; then
  echo -e "${GREEN}✅ PASSED (all integration points present in SKILL.md)${NC}"
  TESTS_PASSED=$((TESTS_PASSED + 1))
else
  echo -e "${RED}❌ FAILED (missing integration points)${NC}"
  TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# Summary
echo "================================================"
echo "  Test Summary"
echo "================================================"
echo ""
echo -e "Tests Passed:  ${GREEN}$TESTS_PASSED${NC}"
echo -e "Tests Failed:  ${RED}$TESTS_FAILED${NC}"
echo -e "Tests Skipped: ${YELLOW}$TESTS_SKIPPED${NC}"
echo ""

TOTAL_TESTS=$((TESTS_PASSED + TESTS_FAILED + TESTS_SKIPPED))
echo "Total: $TOTAL_TESTS tests"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
  echo -e "${GREEN}✅ All tests passed!${NC}"
  exit 0
else
  echo -e "${RED}❌ Some tests failed${NC}"
  exit 1
fi
