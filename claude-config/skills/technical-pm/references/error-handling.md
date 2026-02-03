# Error Handling Protocol

## Overview

When skills fail during orchestrated workflows, technical-pm follows this protocol to ensure graceful degradation, preserve completed work, and provide clear recovery options.

## Error Categories

| Category | Examples | Default Action | Recovery Options |
|----------|----------|----------------|------------------|
| **Transient** | Timeout, rate limit, network error | Retry 2x with backoff | Skip after retries |
| **Skill Error** | Invalid output, crash, assertion failure | Halt and report | Retry, skip, abort |
| **Validation Failure** | Handoff invalid, output missing | Halt immediately | Regenerate, fix, abort |
| **Quality Failure** | Output doesn't meet criteria | Route to user review | Accept, retry, abort |

## Error Response Protocol

### For Sequential Pipelines

When skill N fails in a sequential pipeline:

1. **Preserve outputs**: Skills 1 through N-1 outputs are preserved
2. **Save partial**: If skill N has partial output, save with `_partial` suffix
3. **Report clearly**:

```
Pipeline Status: HALTED at skill 3/5

Completed successfully:
- [x] researcher: docs/literature/review.md
- [x] synthesizer: docs/literature/synthesis.md

Failed:
- [ ] devils-advocate: TIMEOUT after 15 minutes
      Last progress: "Reviewing section 2 of 4"

Not started:
- [ ] fact-checker
- [ ] editor

Completed outputs preserved. Options:
(A) Retry devils-advocate with extended timeout (20 min)
(B) Retry with narrowed scope (review executive summary only)
(C) Skip devils-advocate, continue with fact-checker
(D) Abort pipeline, keep completed outputs
```

### For Parallel Execution

When one of multiple parallel tasks fails:

1. **Wait for others**: Let remaining tasks complete (don't abort prematurely)
2. **Collect all results**: Note which succeeded, which failed
3. **Preserve successful**: Never discard successful task outputs
4. **Report comprehensively**:

```
Parallel Execution Status: PARTIAL SUCCESS (1/2 completed)

Completed:
- [x] researcher: docs/literature/oxygen-review.md

Failed:
- [ ] calculator: TIMEOUT at 40% progress
      Partial output: scratchpad/calculator/membrane-calc_partial.md

Researcher output is complete and preserved.
Calculator partial output available for review.

Options:
(A) Proceed with researcher output only (skip calculator)
(B) Retry calculator with extended timeout
(C) Retry calculator with narrower scope
(D) Abort, keep researcher output
```

## Error Detection

After EVERY skill invocation, validate:

```python
# Pseudo-code for error detection
result = invoke_skill(skill_name, params)

# Check 1: Explicit error status
if result.status == "error":
    handle_skill_error(skill_name, result.error)
    return

# Check 2: Empty or missing output
if not result.output or len(result.output) < MIN_OUTPUT_LENGTH[skill_name]:
    handle_empty_output(skill_name, result)
    return

# Check 3: Output validation against skill criteria
if not validate_skill_output(skill_name, result.output):
    handle_validation_failure(skill_name, result.output)
    return

# Check 4: Handoff document validity
if not validate_handoff(result.handoff):
    handle_handoff_error(skill_name, result.handoff)
    return

# All checks pass - proceed
proceed_to_next_stage(result)
```

## Minimum Output Lengths

| Skill | Minimum Output | Required Content |
|-------|---------------|------------------|
| researcher | 500 chars | At least 1 citation |
| synthesizer | 300 chars | Executive summary section |
| devils-advocate | 200 chars | "Strong points" AND "Challenges" sections |
| fact-checker | 100 chars | "Verified" or "Issues found" statement |
| editor | Original length | Style checklist present |
| calculator | 100 chars | Result with units |

## Retry Protocol

For transient failures, use exponential backoff:

| Attempt | Delay | Action |
|---------|-------|--------|
| 1 | Immediate | First try |
| 2 | 30 seconds | First retry |
| 3 | 2 minutes | Second retry |
| 4 | N/A | Escalate to user |

**NEVER retry automatically**:
- Skills that write to files (may duplicate)
- Skills that have side effects
- Skills that failed validation (not transient)

## Cascading Failure Prevention

If more than 2 skills fail in the same workflow:
1. PAUSE workflow immediately
2. Do NOT continue retrying
3. Escalate to user with full context

```
Multiple failures detected in workflow.

Failed skills:
1. researcher: Timeout
2. calculator: Validation error

This may indicate a systemic issue. Automatic retry disabled.

Options:
(A) Review failures and retry individually
(B) Abort workflow, preserve any outputs
(C) Start fresh with modified goal
```

## Error Context Preservation

Every error report includes:

```yaml
error_context:
  error_id: "err-abc123"
  error_type: SKILL_TIMEOUT | SKILL_ERROR | VALIDATION_FAILURE | QUALITY_FAILURE
  skill: "researcher"
  task_id: "T1"
  workflow_id: "workflow-xyz"

  message: "Skill timed out after 15 minutes"
  timestamp: "2026-02-03T14:15:00Z"

  context:
    last_progress: "Reading paper 5/8"
    elapsed_time: "15m"
    expected_duration: "10m"

  preserved_outputs:
    - "scratchpad/researcher/paper-notes_partial.md"

  recovery_options:
    - id: A
      action: retry
      description: "Retry with extended timeout (20m)"
    - id: B
      action: skip
      description: "Skip this skill, continue workflow"
    - id: C
      action: abort
      description: "Abort workflow, preserve outputs"
```

## Skill-Specific Error Handling

### Researcher Errors

| Error | Likely Cause | Recommended Action |
|-------|--------------|-------------------|
| Timeout | Scope too broad | Narrow to specific papers/topics |
| Empty output | No relevant sources found | Broaden search terms or accept gap |
| Low confidence | Conflicting sources | Flag for user review |

### Calculator Errors

| Error | Likely Cause | Recommended Action |
|-------|--------------|-------------------|
| Validation failure | Missing units | Retry with explicit unit requirements |
| Inconsistent results | Calculation error | Show work, request verification |
| Timeout | Complex computation | Break into smaller calculations |

### Synthesizer Errors

| Error | Likely Cause | Recommended Action |
|-------|--------------|-------------------|
| Missing sections | Incomplete input | Check researcher output completeness |
| Contradictions | Source conflicts | Flag for user resolution |
| Timeout | Too much input | Prioritize key sources |

## User Communication

### Error Notification Format

```
ERROR in workflow [{workflow_id}]

Skill: {skill_name}
Error type: {error_type}
Message: {error_message}

What happened:
{brief_explanation}

What was preserved:
- {list of preserved outputs}

What was lost:
- {list of lost/partial outputs, if any}

Options:
(A) {option_a}
(B) {option_b}
(C) {option_c}

Recommendation: {your_recommendation_with_rationale}
```

### Severity Indicators

Use visual indicators for severity:

- **CRITICAL** (P0): Use all-caps, immediate action required
- **ERROR** (P1): Standard error format, options provided
- **WARNING** (P2): Informational, workflow continues

## Logging

Log all errors for debugging:

```yaml
# Append to /tmp/workflow-errors-{workflow_id}.log
- timestamp: "2026-02-03T14:15:00Z"
  error_id: "err-abc123"
  skill: researcher
  error_type: SKILL_TIMEOUT
  message: "Skill timed out after 15 minutes"
  resolution: "User selected: retry with narrower scope"
  outcome: "SUCCESS on retry"
```

Logs enable:
- Post-workflow analysis
- Pattern detection (which skills fail most?)
- Estimation improvement (are timeouts too aggressive?)
