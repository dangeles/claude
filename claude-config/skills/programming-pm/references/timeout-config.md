# Timeout Configuration

Reference for timeout thresholds and intervention protocols in programming-pm orchestrated projects.

---

## Overview

Timeouts prevent stuck workflows and enable proactive intervention. The programming-pm monitors progress and intervenes when specialists exceed expected durations.

**Principle**: Better to intervene early than wait for complete stall.

---

## Per-Phase Timeouts

Default timeouts for each workflow phase (Phase 0-6).

| Phase | Default Timeout | Warning Threshold | Exceeded Action |
|-------|-----------------|-------------------|-----------------|
| Phase 0: Archival Setup | 5 min | 3 min | ABORT - cannot proceed |
| Phase 1: Requirements | 45 min | 30 min | Escalate to user |
| Phase 2: Pre-mortem | 30 min | 20 min | Proceed with available risks |
| Phase 3: Architecture | 90 min | 60 min | Escalate to user |
| Phase 4: Implementation | Varies by task | Per-task | Per-task escalation |
| Phase 5: Code Review | 30 min/task | 20 min | Skip to next, return later |
| Phase 6: VCS Integration | 15 min | 10 min | Manual intervention |

### Phase-Specific Notes

**Phase 0 (Archival Setup)**:
- Critical phase - workflow cannot proceed without session directory
- ABORT on timeout (no recovery possible)
- If CLAUDE.md not found, use defaults and continue (not a timeout)
- Session directory creation failure is the only timeout trigger
- Creates `/tmp/programming-pm-session-{timestamp}-{pid}/` for session isolation

**Phase 1 (Requirements)**:
- Complex projects may need longer
- User can extend if requirements discussion is productive
- Timeout indicates scope may be too broad

**Phase 2 (Pre-mortem)**:
- Faster phase, focus on critical risks
- If timing out, proceed with identified risks
- Can add risks during implementation if discovered

**Phase 3 (Architecture)**:
- Most variable phase
- Complex systems need more time
- Timeout may indicate need for simplification

**Phase 4 (Implementation)**:
- Per-task timeouts (see specialist section)
- Parallel tasks tracked independently
- Timeout intervention per task, not phase

**Phase 5 (Code Review)**:
- Per-task timeout
- If reviewer unavailable, queue for later
- Don't block other tasks

**Phase 6 (VCS)**:
- Should be quick if previous phases clean
- Timeout usually indicates CI issues
- May need manual intervention

---

## Per-Specialist Timeouts

Default timeouts for specialist tasks.

| Specialist | Default Timeout | Warning Threshold | Exceeded Action |
|------------|-----------------|-------------------|-----------------|
| mathematician | 60 min | 45 min | Substitute with senior-developer |
| statistician | 60 min | 45 min | Proceed without validation (flag) |
| senior-developer | 120 min/task | 90 min | Break task smaller |
| junior-developer | 90 min/task | 60 min | Senior-developer intervention |
| copilot | 30 min | 20 min | Skip assistance, manual review |

### Specialist-Specific Notes

**mathematician**:
- If timing out, problem may be too complex
- Substitution: senior-developer implements with "unverified complexity" flag
- Document for future review

**statistician**:
- If timing out, proceed with implementation
- Flag output as "statistical validation pending"
- Schedule follow-up validation

**senior-developer**:
- Most common timeout scenario
- Usually indicates task scope too broad
- Break into smaller tasks, reassign

**junior-developer**:
- Expected to complete faster with supervision
- Timeout triggers senior-developer takeover
- Document for retrospective (was task too complex?)

---

## Progress Monitoring

### Progress File Format

Each specialist updates a progress file every 15 minutes during active work.

**Location**: `/tmp/progress-{task-id}.md`

```markdown
# Progress: {Task ID}

**Status**: In Progress | Complete | Blocked
**Last Update**: YYYY-MM-DD HH:MM:SS
**Completion**: X%

## Latest Milestone
- Completed: [What was done since last update]
- Current: [What is in progress now]

## Next Steps
1. [Next item] (estimated time)
2. [Next item] (estimated time)

## Blockers
- [If any blockers exist]

## Time Status
- Elapsed: X min
- Estimated remaining: Y min
- Time limit: Z min
```

### Progress Check Frequency

| Phase | Check Frequency | Staleness Threshold |
|-------|-----------------|---------------------|
| Active implementation | Every 15 min | 30 min |
| Code review | Every 10 min | 20 min |
| Architecture | Every 15 min | 30 min |
| All others | Every 10 min | 20 min |

**Staleness**: If progress file not updated within staleness threshold, trigger intervention.

---

## Timeout Detection Protocol

### Detection Steps

1. **Check progress file** at configured frequency
2. **Compare timestamps**:
   - Last update vs. current time
   - Elapsed time vs. warning threshold
   - Elapsed time vs. timeout threshold
3. **Assess completion**:
   - Is completion percentage increasing?
   - Are next steps being completed?
4. **Evaluate blockers**:
   - Are blockers documented?
   - Have blockers changed?

### Detection Outcomes

| Condition | Status | Action |
|-----------|--------|--------|
| Progress updating, on track | Green | Continue monitoring |
| Progress updating, behind schedule | Yellow | Increase check frequency |
| Progress stale, approaching timeout | Orange | Send warning, offer options |
| Progress stale, exceeded timeout | Red | Trigger intervention |
| Blockers documented | Blue | Review blockers, assist |

---

## Timeout Intervention Protocol

When timeout detected, programming-pm executes this protocol.

### Step 1: Diagnose

```yaml
diagnosis:
  task_id: string
  specialist: string
  elapsed_time: int  # minutes
  timeout_threshold: int
  last_progress_update: ISO8601
  completion_percent: int
  documented_blockers: []
  progress_trend: "increasing" | "stalled" | "decreasing"
```

### Step 2: Present Options

Present to user with context:

```markdown
## Timeout Intervention: {Task ID}

**Specialist**: {specialist}
**Elapsed**: {elapsed_time} min (timeout: {timeout_threshold} min)
**Last Update**: {last_progress_update}
**Completion**: {completion_percent}%

### Diagnosis
{Analysis of what might be causing delay}

### Options

1. **Extend Deadline** (+30 min, +1 hour)
   - Use if progress is being made but slower than expected
   - Will continue monitoring

2. **Narrow Scope**
   - Reduce task requirements to achievable subset
   - Mark deferred items for follow-up

3. **Substitute Specialist**
   - Replace with alternative (e.g., senior-developer for mathematician)
   - May result in lower quality output

4. **Escalate to User**
   - Pause task for user guidance
   - User provides direction or takes over

### Recommendation
{programming-pm's suggested option with rationale}
```

### Step 3: Execute Decision

Based on user selection:

**Extend Deadline**:
```yaml
action: extend
new_timeout: {current + extension}
notes: "Extended by user decision"
continue_monitoring: true
```

**Narrow Scope**:
```yaml
action: narrow_scope
original_scope: [...]
reduced_scope: [...]
deferred_items: [...]
follow_up_issue: "ISSUE-XXX"
```

**Substitute Specialist**:
```yaml
action: substitute
original_specialist: string
replacement_specialist: string
handoff_notes: string
quality_flag: "unverified" | "substitute"
```

**Escalate**:
```yaml
action: escalate
paused_at: ISO8601
user_notified: true
awaiting: "user_guidance"
```

### Step 4: Log Decision

Add to exceptions log:

```yaml
exception_log:
  - timestamp: ISO8601
    task_id: string
    type: "timeout"
    diagnosis: {...}
    options_presented: [...]
    decision: string
    outcome: string  # filled in later
```

---

## Circuit Breaker Pattern

After repeated failures, stop automatic retries.

### Trigger Conditions

Circuit breaker opens after **3 consecutive failures** of the same type:
- 3 timeouts on same task
- 3 failed substitutions
- 3 rejected code reviews

### Circuit Breaker States

| State | Behavior |
|-------|----------|
| **Closed** | Normal operation, failures counted |
| **Open** | No automatic retries, user decision required |
| **Half-Open** | Testing if issue resolved (after user intervention) |

### Open Circuit Actions

1. **Stop retrying** automatically
2. **Alert user** with failure summary:
   ```markdown
   ## Circuit Breaker Open: {Task ID}

   **Failure Type**: {timeout | substitution | review}
   **Consecutive Failures**: 3

   ### Failure History
   1. {timestamp}: {failure description}
   2. {timestamp}: {failure description}
   3. {timestamp}: {failure description}

   ### Options
   1. **Retry with Changes** - User modifies task/approach
   2. **Skip Component** - Mark as out of scope, continue
   3. **Abort Workflow** - Stop project, reassess

   Cannot proceed without user decision.
   ```
3. **Require explicit decision** to continue

### Resetting Circuit Breaker

After user intervention and successful retry:
- State transitions to Half-Open
- If next attempt succeeds: Closed
- If next attempt fails: Open again

---

## Global Workflow Timeout

**Maximum workflow duration**: 72 hours (3 days)

### At Global Timeout

1. **Checkpoint current state**:
   ```yaml
   checkpoint:
     timestamp: ISO8601
     workflow_id: string
     phases_completed: []
     current_phase: int
     in_progress_tasks: []
     elapsed_hours: 72
   ```

2. **Prompt user**:
   ```markdown
   ## Workflow Timeout: {workflow_id}

   **Elapsed**: 72 hours
   **Phases Completed**: {list}
   **Current Phase**: {phase}

   ### Options
   1. **Continue** - Extend workflow, continue from checkpoint
   2. **Checkpoint and Pause** - Save state, resume later
   3. **Abort** - End workflow, preserve artifacts

   Note: State is automatically checkpointed.
   ```

3. **If no response within 24 hours**: Auto-checkpoint and pause

---

## Configuration Customization

Users can customize timeouts per project.

### In workflow state file

```yaml
workflow:
  id: string
  custom_timeouts:
    phases:
      phase_1: 60  # override default 45 min
      phase_3: 180  # complex architecture
    specialists:
      mathematician: 90  # complex algorithm
      senior_developer: 180  # large component
    global: 168  # 7 days for large project
```

### Via command line

```bash
# Extended mathematician timeout
programming-pm --timeout mathematician=90 "Complex algorithm project"

# Extended global timeout
programming-pm --global-timeout 168 "Large multi-phase project"
```
