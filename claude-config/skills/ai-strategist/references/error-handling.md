# Error Handling Protocol

This reference defines the error handling, compensation, and recovery protocols for the ai-strategist workflow.

## Saga-Style Compensation Table

Each phase has a defined forward action and compensation action. If a phase fails, the compensation action restores the system to a consistent state.

| Phase | Forward Action | Compensation (on failure) |
|---|---|---|
| 0 | Create session directory and initialize state | Remove session directory (`rm -rf`) |
| 1 | Scope refinement via requirements-analyst | No side effects (read-only interaction) |
| 2 | Launch parallel researcher agents | Cancel running agents, retain completed agent results |
| 3 | Strategic assessment via strategist | Retain Phase 2 outputs; retry once with simplified scope, then escalate |
| 4 | Roadmap synthesis (orchestrator-owned) | Retain Phase 3 outputs; retry once, then escalate |
| 5 | Adversarial review via devils-advocate | Deliver without review (flag deliverable as "unreviewed") |
| 6 | Editorial polish via editor | Deliver pre-polish version of deliverable |

## Circuit Breaker

Adapted from lit-pm and programming-pm patterns.

**Trigger**: 2 consecutive failures of the same agent type (e.g., researcher agent fails twice in a row).

**When circuit opens**:
1. Stop dispatching to that agent type
2. Log the failure pattern in workflow-state.yaml
3. Present options to user:
   - **Retry with different scope**: Narrow the agent's task and retry
   - **Skip agent**: Proceed without that agent's contribution (if minimum thresholds still achievable)
   - **Abort**: Stop the workflow, preserve session for later resume

**Circuit reset**: Only after user selects an option and the selected action succeeds.

## Graceful Cancellation Protocol

Adapted from archive-workflow SIGINT handling pattern.

When the user cancels the workflow (Ctrl+C or explicit abort):

1. **Complete current atomic operation**: Finish any in-progress file write or agent result collection
2. **Update workflow-state.yaml** with current progress (use atomic write protocol below)
3. **Preserve session directory**: Do not clean up -- session is needed for resume
4. **Display resume instructions**:
   ```
   Workflow paused at Phase {N} - {phase_name}.
   Completed phases: {list}
   To resume: ai-strategist --resume {session-id}
   Session directory: {session_path}
   ```
5. **Provide manual cleanup instructions** (if user wants to abandon):
   ```
   To abandon this session: rm -rf {session_path}
   ```

## Atomic State Writes

All workflow-state.yaml updates follow this protocol to prevent corruption:

1. Write updated content to `workflow-state.yaml.tmp`
2. Validate YAML syntax of the temp file
3. If valid: Rename `workflow-state.yaml.tmp` to `workflow-state.yaml` (atomic on POSIX)
4. Maintain single backup: Copy current to `workflow-state.yaml.bak` before rename

**Update frequency**: After each agent completion (not just at phase boundaries). This ensures resume can detect exactly which agents have finished.

**workflow-state.yaml schema**:

```yaml
workflow_id: ai-strategist-session-{timestamp}-{pid}
session_path: /tmp/ai-strategist-session-{timestamp}-{pid}/
invocation_mode: quarterly | deep_dive | event_triggered
current_phase: 2
current_phase_name: Parallel Research
started_at: "2026-01-15T10:00:00Z"
last_updated: "2026-01-15T11:30:00Z"
phases:
  phase_0:
    status: complete
    completed_at: "2026-01-15T10:02:00Z"
  phase_1:
    status: complete
    completed_at: "2026-01-15T10:25:00Z"
    scope_approved: true
  phase_2:
    status: in_progress
    agents:
      agent_1_mcp:
        status: complete
        tools_found: 8
        output_file: research/agent-1-mcp-servers.md
      agent_2_frameworks:
        status: complete
        tools_found: 7
        output_file: research/agent-2-ai-frameworks.md
      agent_3_scientific:
        status: running
        started_at: "2026-01-15T11:00:00Z"
      agent_4_community:
        status: pending
```

## State Recovery Protocol (--resume)

When invoked with `--resume {session-id}`:

1. **Locate session**: Find directory matching `/tmp/ai-strategist-session-{session-id}*`
2. **Read workflow-state.yaml**: Parse current state
3. **Staleness check**: If `last_updated` is more than 72 hours ago, warn user:
   ```
   This session was last active 4 days ago. Tool landscape data may be stale.
   Options: Continue from Phase {N} / Restart from Phase 0
   ```
4. **Display current state**: Show completed phases, current phase, and progress within current phase
5. **Validate output files**: Check that output files referenced in workflow-state.yaml still exist
6. **Handle missing session**: If session directory not found, ABORT with error (cannot recover without session artifacts)
7. **Handle partial Phase 2**: If some agents completed but not all, offer to:
   - Resume with available research (re-launch only incomplete agents)
   - Restart Phase 2 entirely
8. **User choice**: Continue from current phase OR restart from Phase 0

## Global Timeout Handling

| Threshold | Action |
|---|---|
| 5-hour mark | Display warning: "Global workflow timeout approaching (5/6 hours elapsed)" |
| 6-hour mark (Phase 0-3) | Abort with preserved session. Research artifacts are incomplete. |
| 6-hour mark (Phase 4+) | Deliver current state as partial result. Research and assessment are the most valuable artifacts -- deliver those even if roadmap/review/editorial are incomplete. |

On any global timeout:
1. Persist all intermediate files to session directory
2. Write partial-completion summary to `final/partial-completion-summary.md`
3. Retain session directory for potential resume
4. Inform user of what was completed and what remains

## Phase-Specific Error Patterns

### Phase 2: Agent Timeout

If an individual research agent exceeds its timeout (default 30 min):
1. Record timeout in workflow-state.yaml
2. Check if minimum agent threshold is still achievable (need 2+ agents for degraded mode)
3. If achievable: Proceed with available results, inform user of missing category
4. If not achievable: Pause and present options to user

### Phase 2: Agent Produces Empty Output

If an agent returns no tools:
1. Log the empty result
2. Check if the agent's prompt was appropriate (did we send the right search targets?)
3. Offer to retry once with broadened search terms
4. If retry also empty: Proceed without that agent's contribution

### Phase 3: Strategist Cannot Score Tools

If the strategist reports inability to score (missing information, ambiguous criteria):
1. Identify which tools lack sufficient information
2. Option A: Dispatch a targeted researcher agent to gather missing data
3. Option B: Score with available information, flagging low-confidence scores
4. Option C: Exclude under-documented tools from scoring (note exclusion in methodology)

### Phase 5: Methodology Challenge

If the devils-advocate challenges the scoring methodology itself (not just individual scores):
1. Do not auto-resolve -- this is a significant concern
2. Present to user with options:
   - Accept the limitation and document it in the deliverable
   - Re-run Phase 3 with modified methodology
   - Include the critique as an appendix to the deliverable
3. Log the user's decision in workflow-state.yaml
