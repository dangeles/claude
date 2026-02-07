# Model Selection Guide

## Overview

This document provides model selection guidance for each component of the brainstorming-pm workflow. These recommendations are **advisory, not enforced** -- the orchestrator functions correctly with any single model. The guidance exists to inform future optimization when the Task tool supports explicit model selection.

These recommendations apply to the Claude 4.5 model family. Future model releases may warrant re-evaluation. The primary optimization axis is cost-efficiency: use the most capable model only where reasoning complexity demands it.

## Current Behavior

As of 2026-02, the Claude Code Task tool does not support explicit model selection. All components inherit the orchestrator's model.

| Component | Model Used | Why |
|-----------|-----------|-----|
| All components | Orchestrator model (inherited) | Task tool does not support a `model` parameter; all spawned tasks run the same model as the caller |

This means the fallback chains described below reduce to: orchestrator model inline -> retry -> Jaccard.

## Target Architecture

When the Task tool gains model selection support, the following component-level recommendations apply:

| Component | Recommended Model | Rationale | Fallback |
|-----------|------------------|-----------|----------|
| Orchestrator (brainstorming-pm) | Claude Opus 4.5 | Best reasoning for framing, quality gates, conflict resolution | Claude Sonnet 4.5 (slightly reduced quality gate precision) |
| Perspective Agents (Stage 2) | Claude Sonnet 4.5 | Good balance of speed and quality for parallel agents; 5x concurrent makes cost-efficiency important | Claude Opus 4.5 (higher quality but slower, more expensive) |
| LLM Grouping (Stage 3) | Claude Haiku 4.5 | Fast, low-cost; sufficient for semantic clustering of 5 short texts; JSON schema enforcement reliable | Claude Sonnet 4.5 (if Haiku produces low-quality groupings) |
| Workflow Discovery | N/A (file-system scan) | No LLM involved | N/A |
| Relevance Scoring | Orchestrator model (inline) | Simple scoring logic, no separate call needed | N/A |

## Rationale Details

### Perspective Agents (Sonnet 4.5)

The perspective generation task is well-bounded: ~2000 token target, structured output format, 1-2 web searches per agent. Sonnet 4.5 handles this reliably. Running 5 Opus 4.5 agents in parallel would be approximately 5x more expensive with marginal quality improvement for this specific task type. The structured output format (key insight, evidence, confidence, blind spots) constrains the task sufficiently that the additional reasoning depth of Opus provides diminishing returns.

### LLM Grouping (Haiku 4.5)

The grouping task clusters 5 short text snippets (each 1-2 sentences) by thematic similarity. This is well within Haiku's capabilities. The JSON schema enforcement ensures structured output regardless of model tier. Typical expected latency: <2 seconds. The grouping prompt includes explicit "Return ONLY valid JSON" instructions and a schema definition, making the task mechanical rather than requiring deep reasoning.

### Orchestrator (Opus 4.5)

The orchestrator makes judgment calls throughout the workflow: framing quality assessment, conflict resolution during synthesis, quality gate evaluation, and user communication. These tasks benefit from the strongest reasoning model. The orchestrator runs once per session as a single long-lived process, so the cost premium of Opus over Sonnet is amortized across the entire 15-30 minute workflow.

## Fallback Escalation Chains

### Current (Task tool inherits model)

```
Orchestrator model inline (primary)
  -> retry with orchestrator model (1x, 5-second delay)
  -> Jaccard keyword overlap (deterministic fallback)
```

### Target (when Task tool supports model selection)

```
Haiku 4.5 (primary)
  -> retry Haiku 4.5 (1x, 5-second delay)
  -> Sonnet 4.5 (1x escalation)
  -> Jaccard keyword overlap (deterministic fallback)
```

### Fallback Triggers

Fallback is triggered per-component, not globally. See `../perspective-swarm/references/convergence-algorithm.md` (Step 2a) for the complete list of fallback trigger conditions for LLM grouping.

For perspective agents: if an agent produces empty or malformed output, retry with the same model (do not escalate to a higher-tier model). Agent failures are handled by the Stage 2 minimum-agents protocol (>= 4 of 5 required).

## Future: Task Tool Model Parameter

When model selection becomes available in the Task tool, update Stage 2 Task invocations to specify the recommended model:

```
# Placeholder syntax (not yet supported)
Task(model="claude-sonnet-4-5", prompt="...", ...)
```

And update Stage 3 grouping to use:

```
# Placeholder syntax (not yet supported)
Task(model="claude-haiku-4-5", prompt="...", ...)
```

Until then, all components use the orchestrator's inherited model, and the "Current" fallback chain applies.

## References

- [brainstorming-pm/SKILL.md](../SKILL.md) - Orchestrator skill (model selection summary)
- [convergence-algorithm.md](../../perspective-swarm/references/convergence-algorithm.md) - LLM grouping fallback triggers
