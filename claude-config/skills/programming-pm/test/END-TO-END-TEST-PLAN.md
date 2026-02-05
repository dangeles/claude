# End-to-End Test Plan for programming-pm

## Overview

This document outlines the end-to-end testing strategy for the programming-pm workflow with Priority 1 improvements implemented.

**Test Goal**: Verify that all 5 improvements work together in real workflow execution.

---

## Test Prerequisites

### Environment Setup

1. **Install dependencies**:
   ```bash
   # Install yq for YAML parsing (required for quality gates)
   brew install yq  # macOS
   # OR
   pip install yq   # Python version

   # Verify installation
   yq --version
   ```

2. **Verify sync-config.py availability**:
   ```bash
   ls -la ~/repos/claude/sync-config.py
   # Should exist at this path for full Phase 6 testing
   ```

3. **Ensure programming-pm skill is loadable**:
   ```bash
   # In Claude Code, verify skill loads
   # /help should list programming-pm
   ```

---

## Test Scenario 1: Simple Project (SIMPLE Mode)

### Project Specification
**Objective**: Create a CSV parsing utility script

**Characteristics**:
- Single component
- No statistical methods
- No algorithm design
- Estimated: 3-4 tasks
- Expected mode: SIMPLE

### Test Steps

1. **Invoke programming-pm**:
   ```
   I need to create a utility script that parses CSV files and validates data formats.
   Requirements:
   - Read CSV with configurable delimiter
   - Validate column types (string, int, float, date)
   - Output validation report
   - Handle malformed rows gracefully
   ```

2. **Verify Phase 0 (Archival Setup)**:
   - [ ] Session directory created in /tmp/
   - [ ] archival-guidelines-summary.md created
   - [ ] session-handoff.yaml validated
   - [ ] No errors, proceed to Phase 1

3. **Verify Phase 1 (Requirements)**:
   - [ ] requirements-analyst invoked
   - [ ] Requirements document created
   - [ ] Quality Gate 1 automated checks pass
   - [ ] requirements-handoff.yaml validated
   - [ ] User approval received

4. **Verify Mode Selection**:
   - [ ] Complexity detection runs
   - [ ] Detected tier: SIMPLE (high confidence)
   - [ ] Reason: "Single component, no specialization"
   - [ ] 60s timeout prompt displayed
   - [ ] Mode selected: SIMPLE
   - [ ] mode-selection.json created
   - [ ] session-state.json updated with mode

5. **Verify Phase 2 (Pre-Mortem)**:
   - [ ] At least 3 risks identified
   - [ ] Quality Gate 2 passes
   - [ ] premortem-handoff.yaml validated

6. **Verify Phase 3 (Architecture)**:
   - [ ] systems-architect invoked
   - [ ] Single component identified
   - [ ] Quality Gate 3 passes
   - [ ] architecture-handoff.yaml validated

7. **Verify Phase 4 (Implementation - Sequential)**:
   - [ ] SIMPLE mode: Sequential execution confirmed
   - [ ] Task decomposition: 1 task assigned to senior-developer
   - [ ] NO parallel waves (SIMPLE mode characteristic)
   - [ ] specialist completes implementation
   - [ ] code-handoff.yaml created and validated
   - [ ] Quality Gate 4a/4b pass

8. **Verify Phase 5 (Code Review)**:
   - [ ] Automated checks run (ruff, mypy, tests)
   - [ ] Code review completed
   - [ ] Quality Gate 5 passes
   - [ ] review-handoff.yaml validated

9. **Verify Phase 6 (VCS Integration)**:
   - [ ] sync-config.py status checked
   - [ ] All handoffs validated (final check)
   - [ ] Quality Gate 6 passes
   - [ ] Feature branch created
   - [ ] Specific files staged (NOT git add .)
   - [ ] Commit created with conventional format
   - [ ] sync-config.py push executed
   - [ ] Planning journal entry created
   - [ ] Session directory deleted (success)

### Success Criteria
- All phases complete without errors
- Mode selection correctly identifies SIMPLE
- Sequential execution (no parallel waves)
- All handoffs valid
- All quality gates pass
- Session cleaned up on success

---

## Test Scenario 2: Standard Project (STANDARD Mode)

### Project Specification
**Objective**: Build REST API with database integration

**Characteristics**:
- 3 components (API layer, database layer, auth layer)
- No statistical methods (but may need statistician for A/B testing)
- No algorithm design
- Estimated: 8-10 tasks
- Expected mode: STANDARD

### Test Steps

1. **Invoke programming-pm**:
   ```
   Build a REST API for user management with the following:
   - User CRUD endpoints (create, read, update, delete)
   - JWT authentication
   - PostgreSQL database
   - API documentation (OpenAPI/Swagger)
   - Unit and integration tests
   - Docker deployment
   ```

2. **Verify Phase 0-1**: (same as Scenario 1)

3. **Verify Mode Selection**:
   - [ ] Detected tier: STANDARD (high confidence)
   - [ ] Reason: "Multiple components (2-5)"
   - [ ] Mode selected: STANDARD
   - [ ] Parallel execution enabled

4. **Verify Phase 2-3**: (same as Scenario 1)

5. **Verify Phase 4 (Implementation - Parallel)**:
   - [ ] STANDARD mode: Wave-based execution confirmed
   - [ ] Task decomposition: 3 tasks (one per component)
   - [ ] **Wave 1 (T=0s)**: No critical specialists (no math/stats)
   - [ ] **Wave 2 (T=30s)**: 3 senior-developer agents launched in parallel
   - [ ] **Wave 3 (T=60s)**: Dependent tasks (if any)
   - [ ] Progress monitoring active (checking every 60s)
   - [ ] All specialists complete
   - [ ] 3 code-handoff.yaml files created and validated
   - [ ] Quality Gate 4a: 100% completion
   - [ ] Quality Gate 4b passes

6. **Verify Phase 5-6**: (same as Scenario 1)

### Success Criteria
- Mode selection correctly identifies STANDARD
- **Parallel execution in Wave 2** (key test of parallel improvements)
- Multiple code handoffs validated
- All quality gates pass

---

## Test Scenario 3: Extended Project (EXTENDED Mode)

### Project Specification
**Objective**: Statistical analysis pipeline with optimization

**Characteristics**:
- 4 components (data ingestion, statistical analysis, optimization, reporting)
- **Requires BOTH statistics AND mathematics** (dual specialization trigger)
- Estimated: 12-15 tasks
- Expected mode: EXTENDED

### Test Steps

1. **Invoke programming-pm**:
   ```
   Build a statistical analysis pipeline with optimization:
   - Data ingestion from multiple sources
   - Hypothesis testing (t-tests, chi-square)
   - Regression analysis (linear, logistic)
   - Optimization algorithms (gradient descent, simulated annealing)
   - Interactive visualizations
   - Automated reporting
   ```

2. **Verify Phase 0-1**: (same as Scenario 1)

3. **Verify Mode Selection**:
   - [ ] Detected tier: EXTENDED (high confidence)
   - [ ] Reason: "Requires both statistics AND mathematics"
   - [ ] Mode selected: EXTENDED
   - [ ] Parallel execution enabled with extended timeouts

4. **Verify Phase 2-3**: (same as Scenario 1)

5. **Verify Phase 4 (Implementation - Parallel with Critical Specialists)**:
   - [ ] EXTENDED mode: Wave-based execution confirmed
   - [ ] Task decomposition: 4 tasks
   - [ ] **Wave 1 (T=0s)**: mathematician AND statistician launched in parallel (CRITICAL TEST)
   - [ ] **Wave 2 (T=30s)**: 2 senior-developer agents for independent tasks
   - [ ] **Wave 3 (T=60s)**: Dependent tasks
   - [ ] Extended timeouts active (4 hours vs. 2 hours for STANDARD)
   - [ ] All specialists complete
   - [ ] Multiple handoff types: math-handoff.yaml, stats-handoff.yaml, code-handoff.yaml
   - [ ] Quality Gate 4a: Verify critical specialist completion
   - [ ] Quality Gate 4b passes

6. **Verify Phase 5**:
   - [ ] EXTENDED mode: senior-developer reviews ALL code (including senior outputs)
   - [ ] Extended review process

7. **Verify Phase 6**: (same as Scenario 1)

### Success Criteria
- Mode selection correctly identifies EXTENDED
- **Wave 1 launches BOTH mathematician AND statistician** (key test)
- Extended timeouts applied
- Multiple handoff types validated
- Enhanced code review process

---

## Test Scenario 4: Validation Failure Handling

### Project Specification
**Objective**: Test handoff validation blocking behavior

**Test Approach**: Manually inject invalid handoff mid-workflow

### Test Steps

1. **Start workflow** (use Simple or Standard project)

2. **After Phase 1, manually corrupt requirements-handoff.yaml**:
   ```bash
   # Remove required field
   sed -i '' '/problem_statement:/d' ${SESSION_DIR}/handoffs/phase1-requirements-handoff.yaml
   ```

3. **Verify Phase 1→2 boundary**:
   - [ ] Handoff validation runs
   - [ ] Validation FAILS with clear error message
   - [ ] Error: "Missing required field: problem_statement"
   - [ ] User prompted: (A) Fix and retry, (B) Override
   - [ ] Choose (A): Workflow blocks, returns to Phase 1
   - [ ] Fix handoff, retry
   - [ ] Validation PASSES, workflow continues

4. **After Phase 3, test override mechanism**:
   ```bash
   # Remove a less critical field
   sed -i '' '/testing_strategy:/d' ${SESSION_DIR}/handoffs/phase3-architecture-handoff.yaml
   ```

5. **Verify Phase 3→4 boundary**:
   - [ ] Handoff validation runs
   - [ ] Validation FAILS
   - [ ] User prompted: (A) Fix, (B) Override
   - [ ] Choose (B): Override accepted
   - [ ] Override logged to session-state.json
   - [ ] phase3_handoff_override = true
   - [ ] phase3_override_risk = "HIGH"
   - [ ] Workflow continues with warning

6. **Verify Phase 6 (Session Cleanup)**:
   - [ ] Session NOT deleted (has overrides)
   - [ ] Message: "Session preserved for debugging"
   - [ ] session-state.json contains override records

### Success Criteria
- Validation blocks on invalid handoffs
- Clear error messages displayed
- Override mechanism works
- Overrides logged to session state
- Session preserved when overrides present

---

## Test Scenario 5: Specialist Timeout Handling

### Project Specification
**Objective**: Test timeout detection and intervention

**Test Approach**: Simulate specialist timeout (requires manual simulation)

### Test Steps

1. **Start workflow** (use Standard project)

2. **During Phase 4, simulate timeout**:
   ```bash
   # While specialist is running, manually mark start time as 3 hours ago
   # This triggers timeout threshold (2 hours for STANDARD mode)

   # In task-start-times.txt
   echo "TASK-001|$(date -v-3H +%s)" > ${SESSION_DIR}/task-start-times.txt
   ```

3. **Verify timeout detection**:
   - [ ] Progress monitoring detects timeout
   - [ ] Warning displayed: "TASK-001 TIMEOUT (elapsed: 10800s, threshold: 7200s)"
   - [ ] Timeout intervention prompt displayed

4. **Test intervention options**:
   - Option 1: Extend deadline (+30 min, +1 hour)
   - Option 2: Narrow scope (reduce requirements)
   - Option 3: Substitute specialist
   - Option 4: Escalate to user

5. **Choose intervention and verify**:
   - [ ] Intervention applied
   - [ ] Decision logged to session state
   - [ ] Workflow continues or retries

### Success Criteria
- Timeout detected correctly
- Intervention options presented
- Decision logged
- Graceful degradation (no crash)

---

## Test Scenario 6: Backwards Compatibility

### Project Specification
**Objective**: Test legacy session handling

**Test Approach**: Simulate old session without mode-selection.json

### Test Steps

1. **Create legacy session directory**:
   ```bash
   mkdir -p /tmp/legacy-session-test

   # Create session-state.json WITHOUT mode field
   cat > /tmp/legacy-session-test/session-state.json <<EOF
   {
     "session_dir": "/tmp/legacy-session-test",
     "phase": 1,
     "status": "active"
   }
   EOF

   # DO NOT create mode-selection.json
   ```

2. **Resume workflow with legacy session**:
   - [ ] Workflow detects missing mode-selection.json
   - [ ] Warning displayed: "Legacy session (no mode selection)"
   - [ ] Defaults to STANDARD mode
   - [ ] PROGRAMMING_PM_MODE = "STANDARD"
   - [ ] Workflow continues normally

### Success Criteria
- Legacy session detected
- Defaults to STANDARD (safest)
- No crashes or errors
- Workflow continues

---

## Automated Test Suite

### Create automated test runner

```bash
#!/usr/bin/env bash
# test/run-all-tests.sh

set -euo pipefail

echo "Running programming-pm test suite..."
echo ""

# Test 1: Handoff validation
echo "Test 1: Handoff Validation"
python3 scripts/validate-handoff.py test/fixtures/valid-session-handoff.yaml session_handoff
python3 scripts/validate-handoff.py test/fixtures/valid-requirements-handoff.yaml requirements_handoff
echo "✅ Handoff validation tests passed"
echo ""

# Test 2: Mode detection
echo "Test 2: Mode Detection"
bash test/test-mode-detection.sh
echo "✅ Mode detection tests passed"
echo ""

# Test 3: Quality gates (requires yq)
if command -v yq &> /dev/null; then
  echo "Test 3: Quality Gate Validation"
  bash scripts/validate-gate.sh 1 test/fixtures/valid-requirements-handoff.yaml
  echo "✅ Quality gate tests passed"
else
  echo "⚠️  Test 3: Quality Gate Validation SKIPPED (yq not installed)"
fi
echo ""

echo "================================================"
echo "Automated test suite complete"
echo "================================================"
```

---

## Manual Testing Checklist

Use this checklist when running end-to-end tests:

### Pre-Flight Checks
- [ ] yq installed and available
- [ ] sync-config.py available at expected path
- [ ] programming-pm skill loads in Claude Code
- [ ] Test project requirements prepared

### Phase 0 (Archival Setup)
- [ ] Session directory created
- [ ] Archival guidelines extracted
- [ ] session-handoff.yaml created
- [ ] Handoff validates successfully
- [ ] No errors

### Phase 1 (Requirements)
- [ ] requirements-analyst invoked
- [ ] Requirements document created
- [ ] Quality Gate 1 passes
- [ ] requirements-handoff.yaml validates
- [ ] User approval received

### Mode Selection
- [ ] Complexity detection runs
- [ ] Tier detected correctly (SIMPLE/STANDARD/EXTENDED)
- [ ] Confidence level appropriate
- [ ] User prompt displayed with 60s timeout
- [ ] Mode selection recorded
- [ ] mode-selection.json created
- [ ] session-state.json updated

### Phase 2 (Pre-Mortem)
- [ ] Risks identified (>=3)
- [ ] Risk scores calculated
- [ ] Dispositions assigned
- [ ] Quality Gate 2 passes
- [ ] premortem-handoff.yaml validates

### Phase 3 (Architecture)
- [ ] systems-architect invoked
- [ ] Components identified
- [ ] Data flow documented
- [ ] Technology choices justified
- [ ] Quality Gate 3 passes
- [ ] architecture-handoff.yaml validates

### Phase 4 (Implementation)
- [ ] Task decomposition complete
- [ ] Specialists assigned correctly
- [ ] Mode-based execution (sequential vs. parallel)
- [ ] **SIMPLE**: Sequential execution
- [ ] **STANDARD**: Wave 2 parallel execution
- [ ] **EXTENDED**: Wave 1 critical specialists + Wave 2 parallel
- [ ] Progress monitoring active
- [ ] Handoffs created (math/stats/code as needed)
- [ ] All handoffs validate
- [ ] Quality Gate 4a/4b pass

### Phase 5 (Code Review)
- [ ] Automated checks run
- [ ] Code review completed
- [ ] Quality Gate 5 passes
- [ ] review-handoff.yaml validates

### Phase 6 (VCS Integration)
- [ ] sync-config.py status checked
- [ ] All handoffs validated (final)
- [ ] Quality Gate 6 passes
- [ ] Feature branch created
- [ ] Specific files staged
- [ ] Commit created
- [ ] sync-config.py push executed
- [ ] Planning journal entry created
- [ ] Session cleanup (delete or preserve based on overrides)

---

## Success Metrics

### Quantitative Metrics

1. **Test Pass Rate**: Target 100% for automated tests
2. **Mode Detection Accuracy**:
   - False positive rate: <5%
   - False negative rate: <5%
   - User override rate: <20%
3. **Handoff Validation Accuracy**: 100% (all invalid handoffs rejected)
4. **Quality Gate Pass Rate**: Target >90% (some failures expected for testing)

### Qualitative Metrics

1. **User Experience**:
   - Clear error messages
   - Helpful prompts
   - Reasonable timeouts
2. **Graceful Degradation**:
   - Fallback to STANDARD if uncertain
   - Fallback to git if sync-config.py unavailable
   - Warnings for missing dependencies (yq)
3. **Session Management**:
   - Clean on success
   - Preserved on errors for debugging

---

## Known Limitations

1. **yq dependency**: Quality gate validation requires yq for full functionality
   - Fallback: grep-based parsing (limited)
   - Recommendation: Install yq for production use

2. **Specialist invocation**: Full end-to-end testing requires actual agent invocations
   - This document provides test plan
   - Actual execution best done with real projects

3. **Timeout testing**: Requires manual simulation
   - Hard to automate timeout scenarios
   - Best tested in real long-running projects

---

## Conclusion

This test plan provides comprehensive coverage of all Priority 1 improvements:
1. ✅ Mode selection system
2. ✅ Parallel execution framework
3. ✅ Handoff validation at all boundaries
4. ✅ Quality gate automation
5. ✅ sync-config.py integration

**Recommendation**: Execute Scenarios 1-3 with real projects to verify end-to-end behavior.

**Automated testing**: Use `test/run-all-tests.sh` for regression testing.

**Manual testing**: Use checklist above for comprehensive verification.
