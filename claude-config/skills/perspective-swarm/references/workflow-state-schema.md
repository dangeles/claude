# Workflow State Schema

This document defines the YAML schema for `workflow-state.yaml`, enabling session persistence and recovery.

## Schema Definition

```yaml
# workflow-state.yaml
# Version: 1.0

workflow_id: string          # Format: swarm-{YYYYMMDD}-{HHMMSS}-{uuid4-8char}
created_at: ISO8601          # Session creation timestamp
last_updated: ISO8601        # Last state change
expires_at: ISO8601          # Auto-archive after this time (default: +24h)

original_prompt: string      # User's original input
reframed_challenge: string   # Neutrally framed challenge (after Stage 1)
problem_type: enum           # decision | creative | analytical | strategic

current_state: enum          # INITIALIZED | FRAMING | DIVERGING | CONVERGING | AWAITING_USER | COMPLETED | FAILED | ABORTED
global_timeout_at: ISO8601   # 45-minute safety ceiling

stages:
  stage_1:
    status: enum             # pending | in_progress | complete | failed
    started_at: ISO8601
    completed_at: ISO8601
    error: string            # If failed

  stage_2:
    status: enum             # pending | in_progress | complete | failed
    started_at: ISO8601
    deadline_at: ISO8601     # started_at + 15 minutes
    completed_at: ISO8601
    agents:
      optimist:
        status: enum         # pending | in_progress | complete | failed | timeout
        started_at: ISO8601
        completed_at: ISO8601
        output_file: string  # perspectives/optimist.md
        confidence: integer  # 1-5
        search_count: integer
        search_success: boolean
        error: string
        validation_warnings: [string]  # Populated if confidence was clamped
      critic:
        # Same structure
      analyst:
        # Same structure
      innovator:
        # Same structure
      pragmatist:
        # Same structure

  stage_3:
    status: enum
    started_at: ISO8601
    completed_at: ISO8601
    convergent_count: integer
    divergent_count: integer
    error: string

  stage_4:
    status: enum
    user_decision: enum      # accept | refine | handoff
    refinement_feedback: string
    handoff_target: string   # lit-pm

recovery:
  can_resume: boolean
  resume_point: string       # Human-readable description
  completed_perspectives: integer
  minimum_met: boolean       # >= 4 complete
```

## Example State File

```yaml
workflow_id: swarm-20260204-183000-a1b2c3d4
created_at: 2026-02-04T18:30:00Z
last_updated: 2026-02-04T18:45:23Z
expires_at: 2026-02-05T18:30:00Z

original_prompt: "Should we expand into the European market?"
reframed_challenge: "Evaluate market expansion opportunity: European market entry for B2B SaaS product"
problem_type: strategic

current_state: DIVERGING
global_timeout_at: 2026-02-04T19:15:00Z

stages:
  stage_1:
    status: complete
    started_at: 2026-02-04T18:30:00Z
    completed_at: 2026-02-04T18:32:15Z

  stage_2:
    status: in_progress
    started_at: 2026-02-04T18:32:15Z
    deadline_at: 2026-02-04T18:47:15Z
    agents:
      optimist:
        status: complete
        started_at: 2026-02-04T18:32:20Z
        completed_at: 2026-02-04T18:38:45Z
        output_file: perspectives/optimist.md
        confidence: 4
        search_count: 2
        search_success: true
        validation_warnings: []
      critic:
        status: in_progress
        started_at: 2026-02-04T18:32:22Z
      analyst:
        status: complete
        started_at: 2026-02-04T18:32:25Z
        completed_at: 2026-02-04T18:40:12Z
        output_file: perspectives/analyst.md
        confidence: 5
        search_count: 2
        search_success: true
        validation_warnings: []
      innovator:
        status: failed
        started_at: 2026-02-04T18:32:28Z
        error: "WebSearch service unavailable"
        search_success: false
        validation_warnings: []
      pragmatist:
        status: complete
        started_at: 2026-02-04T18:32:30Z
        completed_at: 2026-02-04T18:41:55Z
        output_file: perspectives/pragmatist.md
        confidence: 4
        search_count: 1
        search_success: true
        validation_warnings: []

  stage_3:
    status: pending

  stage_4:
    status: pending

recovery:
  can_resume: true
  resume_point: "Stage 2: Waiting for critic, innovator failed"
  completed_perspectives: 3
  minimum_met: false
```

## State Transitions

```
State Machine:

INITIALIZED ──────────────────────────────────────────────┐
     │                                                    │
     │ start session                                      │
     ▼                                                    │
  FRAMING ────────────────────────────────────────────────┤
     │                                                    │
     │ framing complete                                   │
     ▼                                                    │
 DIVERGING ───────────────────────────────────────────────┤
     │                                                    │
     │ min 4 agents complete                              │
     ▼                                                    │
CONVERGING ───────────────────────────────────────────────┤
     │                                                    │
     │ synthesis complete                                 │
     ▼                                                    │
AWAITING_USER ─────────────────────────────────────────── ┤
     │ accept      │ refine                               │
     │             │                                      │
     ▼             ▼                                      │
 COMPLETED    CONVERGING                                  │
                                                          │
     ◄────────────── user abort ──────────────────────────┤
   ABORTED                                                │
                                                          │
     ◄────────────── unrecoverable error ─────────────────┘
   FAILED
```

## Resume Protocol

On skill invocation:

1. Check for `workflow-state.yaml` files in `/tmp/swarm-session-*/`
2. Check for `.session.lock` files to detect active sessions
3. Filter for `expires_at > now` and `current_state not in [COMPLETED, FAILED, ABORTED]`
4. If found, prompt user:

```
Found incomplete session:
  ID: {workflow_id}
  Started: {created_at}
  Prompt: "{original_prompt}"
  State: {current_state} - {recovery.resume_point}

Options:
(A) Resume session
(B) Start fresh (archives previous)
(C) Cancel
```

5. On resume:
   - Verify `.session.lock` is not held by another process (check if stale > 30 min)
   - Load state file
   - Update `.session.lock` with current workflow
   - Skip completed stages
   - Re-execute in-progress or failed stages
   - For Stage 2: only run incomplete/failed agents
