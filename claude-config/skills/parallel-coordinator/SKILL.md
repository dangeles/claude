---
name: parallel-coordinator
last_updated: 2026-06-09
description: Deprecated as of 2026-06-09. Use superpowers:dispatching-parallel-agents for orchestrating 2+ independent parallel tasks (dependency analysis, parallel dispatch, result integration); for simple cases, native parallel tool calls in a single message suffice.
deprecated: true
deprecated_date: 2026-06-09
replaced_by: superpowers:dispatching-parallel-agents
---

# parallel-coordinator (deprecated)

This skill was deprecated on 2026-06-09 as a duplicate of `superpowers:dispatching-parallel-agents`. That skill covers the same ground: identifying independent tasks, performing dependency analysis, dispatching them in parallel, and integrating the results.

## What to use instead

| Need | Use |
|---|---|
| Orchestrate 2+ independent tasks with dependency analysis and result integration | `superpowers:dispatching-parallel-agents` |
| Fire a few independent operations at once (file reads, searches, API calls) | Native parallel tool calls in a single message — no skill needed |

## Rationale

- Full redundancy with `superpowers:dispatching-parallel-agents`, which already covers task decomposition, dependency verification, parallel dispatch, and synthesis.
- This was a "what-it-does" skill that mostly described a native Claude capability (batching independent tool calls in one message) rather than adding orchestration value.
