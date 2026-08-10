# Model Selection Guide

## Overview

Model selection for each component of the brainstorming-pm workflow. These are
**enforced, not advisory**: pass `model` explicitly on every dispatch. A component that
inherits silently runs at the session's tier, which for an Opus session means five
parallel perspective agents on Opus — the most expensive thing this workflow does.

Use the tier aliases `opus` / `sonnet` / `haiku` rather than pinned version IDs. This
document previously pinned the 4.5 family and went stale twice; aliases track the current
generation on their own.

The optimization axis is cost-efficiency: reach for the most capable tier only where
reasoning complexity demands it.

## Assignments

| Component | Model | Rationale | Fallback |
|-----------|-------|-----------|----------|
| Orchestrator (brainstorming-pm) | `inherit` | Framing, quality gates, conflict resolution, user communication | n/a — runs at session tier |
| Perspective Agents (Stage 2) | `sonnet` | Well-bounded structured output, 5x concurrent makes cost dominate | `opus` (higher quality, ~5x cost) |
| LLM Grouping (Stage 3) | `haiku` | Mechanical clustering of 5 short snippets under a JSON schema | `sonnet` if groupings are low-quality |
| Workflow Discovery | n/a | File-system scan, no LLM involved | n/a |
| Relevance Scoring | inline | Simple scoring, no separate call needed | n/a |

## How to pass it

The Agent/Task tool takes a `model` parameter, and agent definitions take a `model:`
frontmatter field (`inherit` / `sonnet` / `opus` / `haiku`). The tool parameter takes
precedence over frontmatter, so an orchestrator can override a general-purpose agent
per-dispatch:

```
Agent(model="sonnet", subagent_type="pov-perspective-analyst", prompt="...")
Agent(model="haiku", prompt="Cluster these 5 snippets, return only JSON matching ...")
```

When a dedicated agent definition exists and always wants the same tier, set it in that
agent's frontmatter instead and omit the parameter.

## Rationale Details

### Perspective Agents (`sonnet`)

The perspective generation task is well-bounded: ~2000 token target, structured output
format, 1-2 web searches per agent. Sonnet handles this reliably. Running five Opus agents
in parallel costs roughly 5x for marginal quality gain here — the structured output format
(key insight, evidence, confidence, blind spots) constrains the task enough that Opus's
additional reasoning depth hits diminishing returns.

### LLM Grouping (`haiku`)

The grouping task clusters 5 short text snippets (each 1-2 sentences) by thematic
similarity, well within Haiku's range. Schema enforcement guarantees structured output
regardless of tier, and the prompt is explicit ("Return ONLY valid JSON" plus a schema),
making the task mechanical rather than reasoning-heavy. Typical latency: <2 seconds.

### Orchestrator (`inherit`)

The orchestrator makes judgment calls throughout: framing quality assessment, conflict
resolution during synthesis, quality gate evaluation, user communication. These benefit
from the strongest tier available. It runs once per session as a single long-lived process,
so its cost premium is amortized across the whole 15-30 minute workflow — which is exactly
why it is the one component worth leaving at the session tier.

## Fallback Escalation Chain

For LLM grouping:

```
haiku (primary)
  -> retry haiku (1x, 5-second delay)
  -> sonnet (1x escalation)
  -> Jaccard keyword overlap (deterministic fallback)
```

### Fallback Triggers

Fallback is triggered per-component, not globally. See
`../../perspective-swarm/references/convergence-algorithm.md` (Step 2a) for the complete
list of fallback trigger conditions for LLM grouping.

For perspective agents: if an agent produces empty or malformed output, retry at the same
tier rather than escalating. Agent failures are handled by the Stage 2 minimum-agents
protocol (>= 4 of 5 required).

## References

- [brainstorming-pm/SKILL.md](../SKILL.md) - Orchestrator skill (model selection summary)
- [convergence-algorithm.md](../../perspective-swarm/references/convergence-algorithm.md) - LLM grouping fallback triggers
