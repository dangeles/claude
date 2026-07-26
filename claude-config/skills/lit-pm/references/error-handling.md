# Error Handling

Compensation logic, recovery, and escalation for the lit-pm pipeline.

## Contents

- Compensation Matrix
- Repeated-Failure Handling
- Atomic State Writes
- Interrupt Handling and Resume
- Escalation
- Error Logging
- Failure Mode Reference

---

## Compensation Matrix

When a stage fails or needs rolling back:

| Stage | Forward Action | Compensation Action | Trigger |
|-------|----------------|---------------------|---------|
| 1 | Create scope document | Delete scope, reset workflow | User abort |
| 2 | Collect reviews | Delete reviews, retry discovery | Convergence < 2 after retry |
| 3 | Generate outline | Delete outline, return to Stage 2 | User rejects 2x |
| 4 | Write introduction | Delete intro, regenerate from outline | Inconsistent with approved outline |
| 5 | Write sections | Archive failed section, reassign to different agent | Fails validation 3x |
| 6a | Validate section | Return to writer with specific issues | Validation fails |
| 6b | Comprehensive check | Flag issues for manual review | Deep check cannot complete |
| 6c | DA reviews section | Archive review state, pass with uncertainty | Agent crash, persistent disagreement |
| 7 | Synthesize document | Archive synthesis, return to sections with feedback | Cross-section inconsistency detected |
| 7.5 | DA reviews synthesis | Archive review, proceed to Stage 8 with warning | Trigger evaluation failed |
| 8 | Polish document | Re-run editor on specific sections | Major issues found |

When compensation triggers: log the reason, archive the current state rather than deleting it, execute the compensation action, append to `workflow_state.compensation_history[]`, and notify the user if checkpoints are affected.

```yaml
compensation_record:
  timestamp: ISO8601
  stage: integer
  reason: string
  action_taken: string
  state_archived: string  # path to archived state
  user_notified: boolean
```

Stage 6c specifics: on agent crash, retry once, then proceed with an uncertainty note; notify the user on high-stakes work. If two sections in a row end in unresolved disagreement, escalate with the options proceed-with-uncertainty, extend the review, or abort.

Stage 7.5 specifics: if `addition_percentage` cannot be calculated, fall back to complexity tier — run 7.5 for high-stakes, otherwise skip it with a warning.

---

## Repeated-Failure Handling

Rather than retrying indefinitely, each stage has a failure ceiling and a fallback:

| Stage | Ceiling | Fallback when reached |
|-------|---------|-----------------------|
| 2 | All agents fail, or convergence never achieved | Proceed with partial reviews |
| 5 | Validation fails 3x, or the agent crashes | Proceed without the section, flag the gap |
| 6c | 2 consecutive section failures or thesis identification failures | Escalate to user with options: skip 6c to synthesis, review sections manually, or abort |

Failures reset once a stage succeeds.

---

## Atomic State Writes

To prevent workflow state corruption on interrupt or system failure:

1. Write state to `{workflow_id}.state.tmp`
2. Validate: parse the temp file, confirm `workflow_id` present, `stage_current` is an integer, artifact paths exist
3. Rename `{workflow_id}.state.tmp` -> `{workflow_id}.state.yaml` (`rename()` is atomic on POSIX)
4. Copy to `{workflow_id}.state.bak`

If `.state.yaml` fails to parse, try `.state.bak`. If that is also invalid, scan `{scratchpad}/lit-pm/{workflow_id}/` for the most recent `stage_*_complete.yaml` and offer the user: resume from the last complete stage, resume from a stage they pick, or restart preserving only the scope.

---

## Interrupt Handling and Resume

On interrupt: log it, send STOP to running Task agents, collect whatever partial state each returns, write state atomically, and tell the user how to resume.

During parallel section writing (Stage 5), classify each section as COMPLETE (finished and validated), PARTIAL (in progress with a checkpoint available), or NOT_STARTED, and record it:

```yaml
sections_complete: [list of finished sections]
sections_partial:
  section_N:
    completion_percentage: X%
    checkpoint: "description of last milestone"
    scratchpad_path: "/path/to/partial.md"
sections_not_started: [list]
```

Then report the counts and `Resume with: lit-pm --resume {workflow_id}`.

Stage 6c interrupts classify sections as COMPLETE, IN_PROGRESS (challenge document exists), or NOT_STARTED. Stage 7.5 interrupts record `exchanges_completed`, current challenges, and partial resolution.

On resume: load `{workflow_id}.state.yaml`, show a resume summary, and offer to continue from the interrupted point (default), restart a specific stage or section, skip to synthesis with completed sections only (noting the gaps), or view current artifacts.

---

## Escalation

| Condition | Escalation Message |
|-----------|-------------------|
| Stage cannot complete | "Stage {N} could not complete. Options: retry, skip, abort" |
| Failure ceiling reached | "Stage {N} failed {X} times. Proceeding with partial results." |
| Compensation triggered | "Rolling back Stage {N}. Reason: {reason}" |
| Quality floor violated | "Cannot proceed - minimum quality not met. Issue: {issue}" |
| Resource limit reached | "Agent queue full. {X} agents waiting." |

```yaml
escalation:
  type: enum  # timeout | failure_ceiling | compensation | quality_floor | resource
  stage: integer
  timestamp: ISO8601
  message: string
  options:
    - label: string
      action: string
      recommended: boolean
  default_action: string
```

Compensation and quality-floor escalations wait for the user; they do not auto-proceed.

---

## Error Logging

```yaml
error_log:
  timestamp: ISO8601
  workflow_id: string
  stage: integer
  error_type: enum
  severity: enum  # info | warning | error | critical
  message: string
  context:
    agent: string | null
    section: string | null
    attempt: integer
  resolution:
    action_taken: string
    successful: boolean
```

Write logs for a completed workflow to `{workflow_id}.log.yaml`; archive after 30 days.

---

## Failure Mode Reference

| Failure Mode | Detection | Response | Recovery |
|--------------|-----------|----------|----------|
| Agent stalls | No deliverable returned | Log, escalate | Reassign or skip |
| Agent crash | Task tool error | Log, retry once | Reassign |
| YAML parse error | Validation fails | Log, reject handoff | Return to producer |
| State corruption | Load fails | Try backup | Scan artifacts |
| Convergence never achieved | Stage 2 ceiling reached | Proceed with partial | Broaden search |
| Infinite validation loop | Cycle count >= 3 | Escape hatch | User decision |
| Resource exhaustion | Queue full | Escalate | Wait or abort |
| Network failure | WebSearch fails | Retry 3x | Proceed without papers |
