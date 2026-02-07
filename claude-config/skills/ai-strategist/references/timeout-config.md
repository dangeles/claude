# Timeout Configuration

This reference defines per-phase and mode-specific timeout settings for ai-strategist.

## Per-Phase Timeouts

| Phase | Default Timeout | Quarterly | Deep Dive | Event-Triggered | Exceeded Action |
|---|---|---|---|---|---|
| 0 (Archival) | 5 min | 5 min | 5 min | 5 min | ABORT |
| 1 (Scope) | 30 min | 30 min | 20 min | 15 min | Escalate to user |
| 2 (Research total) | 90 min | 90 min | 60 min | 45 min | Proceed with available (min 2 agents) |
| 2 (per agent) | 30 min | 30 min | 30 min | 20 min | Proceed without that agent |
| 2 (brainstorming-pm) | 45 min | 45 min | N/A | N/A | Proceed without (optional) |
| 3 (Assessment) | 60 min | 60 min | 45 min | 30 min | Escalate to user |
| 4 (Roadmap) | 45 min | 45 min | 30 min | 30 min | Escalate to user |
| 5 (Review) | 30 min | 30 min | 20 min | 20 min | Deliver without review (flag risk) |
| 6 (Editorial) | 30 min | 30 min | 20 min | 20 min | Deliver as-is |
| Global | 6 hours | 6 hours | 3 hours | 4 hours | Deliver partial + preserve session |

## Mode-Specific Parameters

| Parameter | Quarterly | Deep Dive | Event-Triggered |
|---|---|---|---|
| Agent count | 4 (all categories) | 1-2 (focused) | 2-3 (gap-specific) |
| Min tools evaluated | 15 | 5 (depth over breadth) | 10 |
| User checkpoints | Phase 1 only | Phase 1 + Phase 4 | Phase 1 only |
| Expected duration | 4-6 hours | 2-3 hours | 2-4 hours |
| brainstorming-pm | Included | Not included | Not included |

## Timeout Handling Behavior

### Phase 0 Timeout (ABORT)
Phase 0 is foundational -- if session setup fails, the workflow cannot proceed. On timeout, abort with error message and no session directory to clean up.

### Phase 2 Agent Timeout
Individual agent timeouts within Phase 2 are handled independently:
1. Record timeout in workflow-state.yaml
2. Mark that agent as "timed_out" in the status board
3. Check if remaining completed agents meet the minimum threshold
4. Continue if threshold met, pause for user decision if not

### Phase 2 Total Timeout
If the total Phase 2 timeout fires while agents are still running:
1. Cancel any remaining running agents
2. Collect results from completed agents
3. Evaluate against quality gate (degraded mode thresholds)
4. Proceed if degraded threshold met

### Global Timeout
See `error-handling.md` for global timeout protocol.

## brainstorming-pm Timeout Note

brainstorming-pm receives a longer per-agent timeout (45 min vs 30 min for standard researchers) because it is itself a Tier 1 orchestrator that runs a multi-stage pipeline internally. However, it is still non-blocking: the workflow proceeds without creative brainstorming output if the timeout fires.

brainstorming-pm is only dispatched in Quarterly mode by default. In Deep Dive and Event-Triggered modes, it is skipped to reduce overall duration.
