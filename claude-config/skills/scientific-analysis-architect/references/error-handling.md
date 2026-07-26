# Error Handling

Failure handling, compensation, and escalation for the scientific-analysis-architect workflow.

## Contents

- Failure Handling Protocol
- Retry and Failure Ceilings
- Compensation Actions by Phase
- Escalation Prompt Templates
- Fan-Out Failure Handling
- Logging
- Graceful Degradation
- Recovery Practices

---

## Failure Handling Protocol

**Detect**: catch Task tool failures, and validate that each agent's output file exists and has the expected format.

**Contain**: isolate the failed agent from its parallel group, preserve partial outputs, and log the details to `{session_dir}/logs/errors.log`.

**Recover**: retry, then compensate per phase.

**Escalate**: use the AskUserQuestion templates below.

## Retry and Failure Ceilings

Retry a failed agent once with the same inputs, logging "Retry attempt 1 for {agent_name}". On a second consecutive failure, escalate to the user. Two attempts per agent per invocation; the count resets on a successful execution.

Wider ceilings:

| Metric | Threshold | Action |
|--------|-----------|--------|
| Consecutive failures per agent | 2 | Stop using that agent, escalate |
| Failures within one phase | 50% of agents | Abort the phase, escalate |
| Session-wide failures | 5 | Suggest aborting the workflow |

## Compensation Actions by Phase

**Phase 0 — initialization.** Session directory creation or output validation failure is terminal: clean up any partial session directory, log the error, and exit with a clear message. No retry.

**Phase 1 — birds-eye planning.** Save any partial `research-structure.md`, log the failure to session state, escalate. Options: retry Phase 1, user supplies a manual research structure, or abort.

**Phase 2 — subsection planning.** Retry each failed consultant once. If the statistician (critical) still fails, escalate with the statistical guidance gap and offer manual guidance, proceeding without, retrying, or aborting. If the mathematician or programmer (optional) still fails, proceed with a warning and log the gap in the analysis plans. For an analysis-planner failure, save partial chapter plans and allow per-chapter retry; escalate if more than 50% of chapters fail.

**Phase 3 — structure review.** Skip the review with a warning. The user can accept the unreviewed structure, retry the review, or abort. If proceeding, carry an "unreviewed" flag into Phase 4 and show the warning at the structure approval gate.

**Phase 4 — plan review.** Mark failed chapters "unreviewed" and proceed with the reviewed ones. The user can accept the partial set, retry the failed chapters, or abort.

**Phase 5 — document generation.** Preserve successfully generated documents, mark failed chapters in session state, and offer per-chapter regeneration, manual creation instructions, or proceeding to Phase 6 with a partial set.

**Phase 6 — statistical fact-checking.** Save interview progress and `corrections-manifest.json` with the decisions made so far, then offer to resume from the last concern, skip the remaining concerns, or pass without statistical review (logged as a warning).

**Phase 7 — audience document generation.** The three documents are independent; one failing does not affect the others. Retry a failed document once, then mark it skipped and continue with the rest (up to 2 retries per document). Preserve what succeeded and log failures to the session state errors array. If all three fail, escalate with Template 1 — retry all, proceed without audience documents, or abort. Phase 0-6 outputs are unaffected either way.

## Escalation Prompt Templates

### Template 1: Critical Agent Failure (Single)

```
[AGENT FAILURE] {agent_name} failed after 2 attempts

Phase: {phase_number} ({phase_name})
Error: {error_type} - {error_message}
Impact: {impact_description}

Options:
  (A) Retry {agent_name} (attempt 3)
  (B) Proceed without {agent_name} output (may affect quality)
  (C) Provide manual input for {missing_output}
  (D) Abort workflow

Enter choice [A/B/C/D]:
```

### Template 2: Multiple Agent Failures (Consolidated)

```
[MULTIPLE FAILURES] {count} agents failed in Phase {phase_number}

Failed agents:
{foreach agent in failed_agents}
  - {agent.name}: {agent.error_type}
{end}

Working agents:
{foreach agent in working_agents}
  - {agent.name}: completed successfully
{end}

Options:
  (A) Retry all failed agents
  (B) Proceed with {working_count} available outputs
  (C) Abort phase and return to Phase {phase_number - 1}
  (D) Abort workflow

Enter choice [A/B/C/D]:
```

### Template 3: Partial Phase Completion

```
[PARTIAL] Phase {phase_number} ({phase_name}) did not complete fully

Status:
  - Completed: {completed_items}
  - Failed: {failed_items}
  - Not started: {not_started_items}

Options:
  (A) Retry the incomplete items
  (B) Proceed with completed items only
  (C) Abort phase and troubleshoot
  (D) Abort workflow

Enter choice [A/B/C/D]:
```

### Template 4: Session-Wide Failure Threshold

```
[WORKFLOW STABILITY WARNING]

This session has encountered {failure_count} failures across phases.
Failure rate: {failure_rate}%

Recent failures:
{foreach failure in recent_failures limit=5}
  - Phase {failure.phase}: {failure.agent} - {failure.type}
{end}

Recommendation: Consider aborting and troubleshooting.

Options:
  (A) Continue workflow (accepting higher risk)
  (B) Save session state and pause for troubleshooting
  (C) Abort workflow and start fresh

Enter choice [A/B/C]:
```

## Fan-Out Failure Handling

When several agents run in parallel (Phase 2 consultants, Phase 4 and 5 per-chapter), wait for all of them to return before escalating, collect every failure into one report, present Template 2, and let the user decide once for the whole group.

| Agent | Criticality | Failure behavior |
|-------|-------------|------------------|
| statistician-consultant | Critical | Escalate if fails |
| mathematician-consultant | Optional | Warn and proceed |
| programmer-consultant | Optional | Warn and proceed |
| notebook-reviewer (per chapter) | Critical | Per-chapter escalation |
| notebook-generator (per chapter) | Critical | Per-chapter escalation |

If more than 50% of a fan-out fails, abort the current phase and escalate with the options retry all, return to the previous phase, or abort the workflow.

## Logging

`{session_dir}/logs/errors.log`, one JSON object per line:

```json
{
  "timestamp": "2026-02-04T14:30:22Z",
  "session_id": "session-20260204-143022-12345",
  "phase": 2,
  "agent": "statistician-consultant",
  "error_type": "no_output",
  "error_message": "Agent returned without writing consultation file",
  "retry_count": 2,
  "resolution": "user_escalation",
  "user_choice": "proceed_without"
}
```

Levels: **ERROR** for agent failure or validation failure, **WARNING** for optional agent failure or partial completion, **INFO** for retry attempts, **DEBUG** for agent input/output details in verbose mode.

## Graceful Degradation

Degradation paths, in order of decreasing completeness: full workflow; optional consultants missing but statistician succeeded; some chapters failed generation or review; no statistical review (user accepts the risk); no audience documents.

When a degraded workflow completes, summarize what is missing:

```
Workflow Complete (with degradation)

Warnings:
- Chapter 3 analysis documents generated without mathematician input
- Chapter 4 not statistically reviewed

Generated outputs:
- 3 of 4 chapters complete
- 8 of 11 planned analysis documents

Recommendation: Manually review Chapter 3 algorithms and Chapter 4 statistics.
```

## Recovery Practices

Retry on a single isolated failure, a transient error, or a first occurrence in the session. Abort when the same agent fails 3+ times, when more than 50% of a parallel group fails, or when the user asks.

Checkpoint session state before spawning parallel agents, before user approval gates, and before document generation, so that user-provided context, approval decisions, and statistical interview progress survive a failure.

When a session aborts, report:

```
Session Summary

Completed: Phases {completed_list}
Failed at: Phase {phase} - {agent}

Saved artifacts:
- {session_dir}/session-state.json
- {session_dir}/logs/errors.log
- {list of generated outputs}

To resume:
  /scientific-analysis-architect

To start fresh:
  rm -rf {session_dir} && /scientific-analysis-architect
```
