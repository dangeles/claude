# Workflow State Management

## Overview

technical-pm maintains workflow state to enable:
- Resume after interruption (Ctrl+C, session timeout)
- Progress visibility during execution
- Debugging of failed workflows
- Partial result preservation

## State Machine

### States

| State | Description | Next States |
|-------|-------------|-------------|
| `INITIALIZING` | Parsing user goal, identifying required skills | PLANNING, FAILED |
| `PLANNING` | Building execution plan with dependencies | EXECUTING, FAILED |
| `EXECUTING` | Running skills (sequential or parallel) | WAITING, SYNTHESIZING, FAILED |
| `WAITING` | Blocked on skill completion or user input | EXECUTING, FAILED, PAUSED |
| `SYNTHESIZING` | Combining outputs from multiple skills | COMPLETED, FAILED |
| `COMPLETED` | All skills finished, final output ready | (terminal) |
| `FAILED` | Unrecoverable error occurred | (terminal) |
| `PAUSED` | User-requested pause for review | EXECUTING, COMPLETED |

### State Transitions

```
INITIALIZING -> PLANNING      (goal parsed successfully)
PLANNING -> EXECUTING         (execution plan created)
EXECUTING -> WAITING          (skill launched, awaiting completion)
WAITING -> EXECUTING          (skill completed, more skills to run)
WAITING -> SYNTHESIZING       (all skills complete, ready to combine)
SYNTHESIZING -> COMPLETED     (final output generated)
Any -> FAILED                 (unrecoverable error)
Any -> PAUSED                 (user interrupt with Ctrl+C)
PAUSED -> EXECUTING           (user resumes workflow)
```

## State File Format

State is persisted to enable recovery:

**Location**: `/tmp/workflow-state-{workflow_id}.yaml`

```yaml
workflow:
  id: "workflow-abc123"
  status: executing
  started_at: "2026-02-03T14:00:00Z"
  updated_at: "2026-02-03T14:15:00Z"

goal:
  original: "Write comprehensive literature review on hepatocyte oxygenation"
  parsed_tasks: ["research", "synthesize", "review", "edit"]

execution_plan:
  - task_id: T1
    skill: researcher
    status: complete
    output: "docs/literature/hepatocyte/review-draft.md"
    started_at: "2026-02-03T14:01:00Z"
    completed_at: "2026-02-03T14:45:00Z"

  - task_id: T2
    skill: synthesizer
    status: in_progress
    partial_output: "scratchpad/synthesizer/draft-wip.md"
    progress: "50%"
    started_at: "2026-02-03T14:46:00Z"

  - task_id: T3
    skill: devils-advocate
    status: pending
    depends_on: [T2]

  - task_id: T4
    skill: editor
    status: pending
    depends_on: [T3]

outputs:
  T1:
    location: "docs/literature/hepatocyte/review-draft.md"
    handoff: "/tmp/handoff-T1-T2.yaml"
```

## State Persistence Protocol

1. **Write state** after EVERY state transition
2. **Atomic writes**: Write to `.tmp` file, then rename
3. **Backup**: Keep last 3 state files for debugging
4. **Cleanup**: Remove state files after COMPLETED or explicit abort

### Atomic Write Pattern

```python
# Pseudo-code for atomic state persistence
def save_state(workflow_id, state):
    state_file = f"/tmp/workflow-state-{workflow_id}.yaml"
    tmp_file = f"{state_file}.tmp"

    # Write to temporary file
    with open(tmp_file, 'w') as f:
        yaml.dump(state, f)

    # Atomic rename
    os.rename(tmp_file, state_file)

    # Maintain backup rotation
    rotate_backups(state_file, keep=3)
```

## Resume Protocol

On session start, check for existing state files:

```
Found interrupted workflow: workflow-abc123
Status: EXECUTING (50% complete)
Last activity: 30 minutes ago

Completed:
- [x] researcher: docs/literature/hepatocyte/review-draft.md

In progress:
- [ ] synthesizer: 50% complete (partial saved)

Pending:
- [ ] devils-advocate
- [ ] editor

Options:
(A) Resume from synthesizer
(B) Restart synthesizer from beginning
(C) Abort workflow (keep completed outputs)
(D) Start fresh (discard all progress)
```

### Resume Decision Tree

```
Is there an existing state file for this goal?
  No  -> Start fresh
  Yes -> Check status
         |
         +-> COMPLETED -> Report "Already done", show outputs
         +-> FAILED -> Show error, offer: retry failed step or start fresh
         +-> PAUSED -> Offer: resume or abort
         +-> EXECUTING/WAITING -> Session was interrupted
              |
              +-> Last activity <1 hour ago -> Likely user pause, offer resume
              +-> Last activity >1 hour ago -> Likely crash, offer resume/fresh
```

## Cancellation Handling (Ctrl+C / SIGINT)

When user interrupts:

1. Set status to `PAUSED`
2. Save current state immediately
3. Preserve all completed outputs
4. Save partial outputs with `_interrupted` suffix
5. Display recovery instructions

```
Workflow paused.

Completed outputs preserved:
- docs/literature/hepatocyte/review-draft.md

Partial outputs saved:
- scratchpad/synthesizer/draft-wip_interrupted.md (50%)

To resume: Invoke technical-pm with resume=true
To abort: Invoke technical-pm with abort=true
```

## Session Timeout Handling

If session times out mid-workflow:

1. State file persists in `/tmp/` (survives session)
2. Completed outputs in user directories persist
3. Partial outputs in scratchpad may be lost (volatile)

**Mitigation**: For long tasks, checkpoint partial work to user directory, not scratchpad.

## Progress Tracking

Update state file with progress indicators:

```yaml
execution_plan:
  - task_id: T2
    skill: synthesizer
    status: in_progress
    progress: "50%"
    progress_detail: "Synthesized sections 1-3 of 6"
    last_checkpoint: "2026-02-03T14:30:00Z"
    estimated_remaining: "15 minutes"
```

Progress enables:
- User visibility into long-running tasks
- Informed resume decisions (restart from 0% or 50%?)
- Timeout detection (no progress update in 30min)

## State File Discovery

On workflow start, check for existing state:

```bash
# Find existing workflow states
ls /tmp/workflow-state-*.yaml

# Parse goal from state to match against current request
for state_file in /tmp/workflow-state-*.yaml; do
    goal=$(yaml_extract goal.original $state_file)
    if similar_goal("$goal", "$current_request"); then
        offer_resume($state_file)
    fi
done
```

## Cleanup Protocol

After COMPLETED or explicit abort:

1. Archive state file to `~/.claude/workflow-archive/` (optional, for debugging)
2. Remove `/tmp/workflow-state-{id}.yaml`
3. Remove `/tmp/handoff-*.yaml` files for this workflow
4. Keep all output files in user directories

**Never auto-delete user outputs**. Only cleanup temporary orchestration files.
