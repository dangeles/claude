# Error Handling

Compensation logic, circuit breakers, and recovery procedures for the lit-pm pipeline.

---

## Compensation Matrix (Saga Pattern)

When a stage fails or needs to be rolled back, use these compensation actions:

| Stage | Forward Action | Compensation Action | Trigger |
|-------|----------------|---------------------|---------|
| 1 | Create scope document | Delete scope, reset workflow | User abort |
| 2 | Collect reviews | Delete reviews, retry discovery | Convergence < 2 after retry |
| 3 | Generate outline | Delete outline, return to Stage 2 | User rejects 2x |
| 4 | Write introduction | Delete intro, regenerate from outline | Inconsistent with approved outline |
| 5 | Write sections | Archive failed section, reassign to different agent | Fails validation 3x |
| 6a | Validate section | Return to writer with specific issues | Validation fails |
| 6b | Comprehensive check | Flag issues for manual review | Deep check timeout |
| 7 | Synthesize document | Archive synthesis, return to sections with feedback | Cross-section inconsistency detected |
| 8 | Polish document | Re-run editor on specific sections | Major issues found |

### Compensation Protocol

When compensation is triggered:
1. Log compensation reason with timestamp
2. Archive current state (do NOT delete)
3. Execute compensation action
4. Update workflow_state.compensation_history[]
5. Notify user if compensation affects checkpoints

```yaml
compensation_record:
  timestamp: ISO8601
  stage: integer
  reason: string
  action_taken: string
  state_archived: string  # path to archived state
  user_notified: boolean
```

---

## Circuit Breaker Configuration

Circuit breakers prevent infinite loops and resource exhaustion.

### Stage 5 Circuit Breaker

```yaml
circuit_breaker:
  stage_5:
    failure_threshold: 3        # Open after 3 failures
    timeout_seconds: 21600      # 6 hours
    half_open_after: 1800       # Try again after 30 min cooldown
    action_when_open: proceed_without_section

    failure_definition:
      - validation_fails_3x
      - timeout_exceeded
      - agent_crash
```

### Stage 2 Circuit Breaker

```yaml
circuit_breaker:
  stage_2:
    failure_threshold: 2        # Open after 2 failures
    timeout_seconds: 7200       # 2 hours
    half_open_after: 900        # Try again after 15 min
    action_when_open: proceed_with_partial_reviews

    failure_definition:
      - all_agents_fail
      - convergence_never_achieved
      - timeout_exceeded
```

### Circuit Breaker States

- **CLOSED**: Normal operation, requests proceed
- **OPEN**: Stage failed too many times, skip or escalate
- **HALF-OPEN**: After cooldown, try one request to test recovery

### State Transitions

```
CLOSED --[failure_threshold reached]--> OPEN
OPEN --[half_open_after elapsed]--> HALF-OPEN
HALF-OPEN --[success]--> CLOSED
HALF-OPEN --[failure]--> OPEN
```

---

## Atomic State Writes Protocol

To prevent workflow state corruption on Ctrl+C or system failure:

### Write Protocol

1. Write state to `{workflow_id}.state.tmp`
2. Validate: Parse temp file, verify all required fields present
3. Atomic rename: `{workflow_id}.state.tmp` -> `{workflow_id}.state.yaml`
4. Backup: Copy to `{workflow_id}.state.bak`

### Implementation

```yaml
atomic_write:
  steps:
    1_write_temp:
      path: "{workflow_id}.state.tmp"
      content: serialized_state

    2_validate:
      action: parse_yaml
      verify:
        - workflow_id present
        - stage_current is integer
        - artifacts paths exist

    3_atomic_rename:
      from: "{workflow_id}.state.tmp"
      to: "{workflow_id}.state.yaml"
      note: "rename() is atomic on POSIX systems"

    4_backup:
      copy_to: "{workflow_id}.state.bak"
```

### Recovery on Corruption

If `.state.yaml` is invalid:
1. Try `.state.bak`
2. If .bak also invalid, scan artifacts directory for latest valid checkpoint
3. Present user with recovery options:
   - Resume from last checkpoint
   - Resume from specific stage
   - Restart workflow

```yaml
recovery_protocol:
  on_corruption:
    1_try_backup:
      source: "{workflow_id}.state.bak"
      action: validate_and_load

    2_scan_artifacts:
      directory: "{scratchpad}/lit-pm/{workflow_id}/"
      find: "stage_*_complete.yaml"
      select: most_recent

    3_user_options:
      present:
        - "Resume from Stage {last_complete + 1}"
        - "Resume from Stage {user_selected}"
        - "Restart workflow (preserve scope only)"
```

---

## Interrupt Handling (Ctrl+C)

### Global Interrupt Handler

```yaml
on_interrupt:
  1_log: "Interrupt received at {timestamp}"
  2_signal_agents: "Send STOP to all running Task agents"
  3_wait_graceful: 30 seconds
  4_collect_states: "Gather partial state from each agent"
  5_write_state: "Atomic write with partial state"
  6_notify_user: "Workflow paused. Resume with --resume {workflow_id}"
```

### Stage 5 Parallel Execution Interrupt

When interrupt occurs during parallel section writing:

```yaml
on_interrupt:
  1_signal_agents: Send STOP to all running agents
  2_wait_graceful: 30 seconds for partial state save
  3_collect_states:
    - COMPLETE: Section finished, validated
    - PARTIAL: Section in progress, checkpoint available
    - NOT_STARTED: Section queued, no progress
  4_write_state:
    sections_complete: [list of finished sections]
    sections_partial:
      section_N:
        completion_percentage: X%
        checkpoint: "description of last milestone"
        scratchpad_path: "/path/to/partial.md"
    sections_not_started: [list]
  5_notify_user: |
    Workflow paused.
    - {X} sections complete
    - {Y} sections partial (checkpointed)
    - {Z} sections not started

    Resume with: lit-pm --resume {workflow_id}
```

### Resume Behavior

On resume:
1. Load workflow state from `{workflow_id}.state.yaml`
2. Display resume summary
3. Offer options:
   - Continue from interrupted point (default)
   - Restart specific stage
   - View current artifacts

```yaml
resume_options:
  - option: "Continue from Section {N} (in-progress)"
    action: resume_section_from_checkpoint
  - option: "Restart Section {N} from outline"
    action: restart_section
  - option: "Skip to synthesis with completed sections only"
    action: proceed_with_partial
    warning: "Not recommended - may have gaps"
```

---

## Escalation Procedures

### When to Escalate to User

| Condition | Escalation Message |
|-----------|-------------------|
| Stage timeout exceeded | "Stage {N} exceeded timeout ({X} min). Options: extend, skip, abort" |
| Circuit breaker opens | "Stage {N} failed {X} times. Proceeding with partial results." |
| Compensation triggered | "Rolling back Stage {N}. Reason: {reason}" |
| Quality floor violated | "Cannot proceed - minimum quality not met. Issue: {issue}" |
| Resource limit reached | "Agent queue timeout. {X} agents waiting." |

### Escalation Format

```yaml
escalation:
  type: enum  # timeout | circuit_breaker | compensation | quality_floor | resource
  stage: integer
  timestamp: ISO8601
  message: string
  options:
    - label: string
      action: string
      recommended: boolean
  default_action: string
  auto_proceed_after: integer | null  # seconds, or null for no auto
```

### User Response Handling

If user doesn't respond within timeout:
1. Log "User response timeout"
2. Take default action (if safe)
3. Continue workflow with note in state

```yaml
response_timeout:
  wait_seconds: 300  # 5 minutes for interactive
  default_action_policy:
    timeout: "extend_by_50%"
    circuit_breaker: "proceed_with_partial"
    compensation: "wait_for_user"  # Do NOT auto-proceed
    quality_floor: "wait_for_user"  # Do NOT auto-proceed
```

---

## Error Logging

### Log Format

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

### Log Retention

- Current workflow: Keep all logs in memory
- Completed workflow: Write to `{workflow_id}.log.yaml`
- Archive after 30 days

---

## Recovery Testing Checklist

Before deploying lit-pm, verify these recovery scenarios work:

- [ ] Ctrl+C during Stage 2 parallel execution -> resumes correctly
- [ ] Ctrl+C during Stage 5 parallel section writing -> partial sections preserved
- [ ] Stage 5 section fails validation 3x -> circuit breaker opens, proceeds without section
- [ ] Stage 3 user rejects 2x -> compensation triggers, returns to Stage 2
- [ ] State file corrupted -> recovery from backup works
- [ ] Agent timeout -> timeout action executes correctly
- [ ] Resource limit reached -> queue behaves correctly
- [ ] Full workflow resume after system restart

---

## Failure Mode Reference

| Failure Mode | Detection | Response | Recovery |
|--------------|-----------|----------|----------|
| Agent timeout | Timer exceeds limit | Log, escalate | Reassign or skip |
| Agent crash | Task tool error | Log, retry once | Reassign |
| YAML parse error | Validation fails | Log, reject handoff | Return to producer |
| State corruption | Load fails | Try backup | Scan artifacts |
| Convergence never achieved | Stage 2 timeout | Open circuit | Proceed with partial |
| Infinite validation loop | Cycle count >= 3 | Escape hatch | User decision |
| Resource exhaustion | Queue timeout | Escalate | Wait or abort |
| Network failure | WebSearch fails | Retry 3x | Proceed without papers |
