---
name: pov-expansion
description: Use when seeking analogous solutions from other domains, when stuck on a problem and need fresh perspectives, or when evaluating whether approaches from field X might apply to field Y. Requires structured problem statement.
---

# POV-Expansion: Cross-Domain Perspective Generation

Generate systematic cross-domain perspectives through boundary spanning research frameworks. Uses structure-mapping theory, convergence tracking, and transfer evaluation to identify analogous problems and transferable solutions.

## When to Use

Use this skill when:
- **Stuck on a problem** - Need fresh perspectives from outside your domain
- **Evaluating analogies** - Want to assess if approach from field X applies to field Y
- **Seeking breakthrough insights** - Ready for 2-3 hour deep exploration
- **Have structured problem** - Can articulate problem clearly enough for abstraction

## When NOT to Use

Do NOT use when:
- **Exploratory browsing** - No specific problem statement (use /research-pipeline instead)
- **Time-constrained** - Need answer in <1 hour (this workflow takes 2-3 hours minimum)
- **Well-understood problem** - Already know the solution domain (direct research faster)
- **Vague problem** - Cannot articulate problem clearly (refine first)

## Architecture

**Three-tier architecture**:

1. **Orchestrator**: pov-expansion-pm (coordinates workflow, handles exceptions)
2. **Specialist Agents**: 5 new agents for abstraction, perspective generation, synthesis, transfer evaluation
3. **Supporting Skills**: Reuses requirements-analyst, fact-checker, devils-advocate, editor, strategist

## Pre-flight Validation

Before invoking `/pov-expansion`, verify:
- [ ] Problem statement is specific (not "improve product" but "reduce onboarding time")
- [ ] Expected runtime: 2-3 hours (workflow cannot be rushed)
- [ ] Dependencies available: requirements-analyst, fact-checker, devils-advocate, editor, strategist skills exist

## The 11-Stage Pipeline

### Stage 1: Problem Refinement
**Owner**: requirements-analyst (reused)
**Duration**: 15-30 min
**Quality Gate**: Specific problem statement with measurable criteria

### Stage 2: Problem Abstraction
**Owner**: pov-abstractor-classifier
**Duration**: 15-30 min
**Quality Gate**: Rasmussen hierarchy populated (5/5 levels)

### Stage 3: Domain Classification
**Owner**: pov-abstractor-classifier
**Duration**: 20-45 min (includes user confirmation)
**Quality Gate**: 4 domains proposed (1+ Near, 1+ Mid, 1+ Far), user confirmed

### Stage 4: Parallel Perspective Generation
**Owner**: pov-perspective-analyst (x4 parallel)
**Duration**: 30-60 min (parallel)
**Quality Gate**: 3+ perspectives complete, each with structural analogies

### Stage 5: Parallel Verification
**Owner**: fact-checker (x2 parallel, reused)
**Duration**: 15-30 min (parallel)
**Quality Gate**: No factual errors in central claims

### Stage 6: Convergence Analysis
**Owner**: pov-synthesizer
**Duration**: 20-40 min
**Quality Gate**: Convergence floor met (1+ theme in 2+ reports)

### Stage 7: Transfer Evaluation
**Owner**: pov-transfer-evaluator
**Duration**: 20-40 min
**Quality Gate**: Transfer feasibility scored, Far field warnings added

### Stage 8: Master Synthesis
**Owner**: pov-synthesizer
**Duration**: 45-90 min
**Quality Gate**: All perspectives integrated, convergent insights highlighted

### Stage 9: Strategic Assessment
**Owner**: strategist (reused skill)
**Duration**: 20-40 min
**Quality Gate**: 3+ actionable recommendations with implementation paths

### Stage 10: Adversarial Review
**Owner**: devils-advocate (reused)
**Duration**: 20-40 min
**Quality Gate**: Critical issues identified, mitigations proposed

### Stage 11: Editorial Polish
**Owner**: editor (reused)
**Duration**: 15-30 min
**Quality Gate**: Consistent voice, no substantive errors

## Timeout Configuration

| Stage | Timeout | Exceeded Action |
|-------|---------|-----------------|
| 1 (Refinement) | 30 min | Escalate to user |
| 2 (Abstraction) | 30 min | Escalate to user |
| 3 (Classification) | 45 min | Proceed with HIGH confidence domains only |
| 4 (Perspectives) | 60 min total | Proceed with available (minimum 2 required) |
| 4 (per agent) | 20 min | Retry once, then skip |
| 5 (Verification) | 30 min | Skip verification (flag risk) |
| 6 (Convergence) | 40 min | Proceed with portfolio framing |
| 7 (Transfer) | 40 min | Lower confidence, proceed |
| 8 (Synthesis) | 90 min | Escalate to user |
| 9 (Strategic) | 40 min | Skip (optional) |
| 10 (Review) | 40 min | Deliver without review (flag risk) |
| 11 (Editorial) | 30 min | Deliver as-is |

**Global Workflow Timeout**: 8 hours
**On Global Timeout**: Save state, escalate to user

## Quality Gate Specifications

### Gate 1: Problem Abstraction (Stage 2)

| Check | Threshold | Automated? |
|-------|-----------|------------|
| Rasmussen levels populated | 5/5 levels | Yes |
| Structural elements identified | >= 3 elements | Yes |
| Constraints defined | >= 1 constraint | Yes |
| Success metrics measurable | >= 1 metric | Manual |

**Pass Threshold**: 4/4 checks pass
**On Failure**: Return to Stage 1 with feedback (max 2 iterations)

### Gate 2: Domain Classification (Stage 3)

| Check | Threshold | Automated? |
|-------|-----------|------------|
| Domains proposed | 4 domains | Yes |
| Near field present | >= 1 | Yes |
| Far field present | >= 1 | Yes |
| User confirmation | Required | No |

**Pass Threshold**: All checks pass
**On Failure**: Re-run with user guidance

### Gate 3: Perspective Generation (Stage 4)

| Check | Threshold | Automated? |
|-------|-----------|------------|
| Perspectives complete | >= 3 of 4 | Yes |
| Structural analogies per perspective | >= 2 | Yes |
| Domain-specific examples | >= 2 per perspective | Yes |
| No placeholder text | 0 instances | Yes |

**Pass Threshold**: All checks pass
**Marginal Pass**: 3/4 checks (flag in synthesis)
**On Failure**: Retry failed perspective or proceed with reduced coverage

### Gate 4: Transfer Feasibility (Stage 7)

| Check | Threshold | Automated? |
|-------|-----------|------------|
| All perspectives evaluated | 100% | Yes |
| Structural similarity scored | All have scores | Yes |
| At least one high-transfer analogy | >= 1 with score > 0.6 | Yes |
| Far field warnings added | All Far field have warnings | Yes |

**Pass Threshold**: 4/4 checks pass
**Marginal Pass**: 3/4 (proceed with caution flag)
**On Failure**: Return to Stage 3 with broadened domains

### Gate 5: Synthesis Coherence (Stage 8)

| Check | Threshold | Automated? |
|-------|-----------|------------|
| All perspectives referenced | 100% | Yes |
| Convergent insights highlighted | >= 1 | Yes |
| Contradictions resolved | 0 unresolved | Manual |
| Executive summary present | Yes | Yes |

**Pass Threshold**: All checks pass
**On Failure**: Return to synthesis with feedback

### Gate 6: Strategic Actionability (Stage 9)

| Check | Threshold | Automated? |
|-------|-----------|------------|
| Actionable recommendations | >= 3 | Yes |
| Implementation paths | Each has path | Yes |
| Timeline/effort estimates | Present | Yes |

**Pass Threshold**: All checks pass
**On Failure**: Add missing elements or skip stage

## Session Directory Structure

Location: `/tmp/pov-session-{workflow-id}/`

```
/tmp/pov-session-{workflow-id}/
├── workflow-state.yaml          # Workflow state for resume
├── progress.md                  # Progress tracking
├── stage-1-refinement.md        # Stage 1 output
├── stage-2-abstraction.yaml     # Stage 2 output (structured)
├── stage-3-domains.yaml         # Stage 3 output (structured)
├── stage-4-perspectives/        # Parallel outputs
│   ├── near-field-1.md
│   ├── near-field-2.md
│   ├── mid-field.md
│   └── far-field.md
├── stage-5-verification/        # Parallel outputs
│   ├── verification-1.md
│   └── verification-2.md
├── stage-6-convergence.md
├── stage-7-transfer.md
├── stage-8-synthesis.md
├── stage-9-strategic.md
├── stage-10-review.md
├── stage-11-final.md            # Final deliverable
└── failure-log.md               # Failure tracking
```

## Workflow State Schema

```yaml
workflow_state:
  workflow_id: string
  problem_statement: string
  started_at: ISO8601
  current_stage: integer  # 1-11
  stages_completed: [integer]

  artifacts:
    stage_1: path | null
    stage_2: path | null
    # ... etc

  status: enum[in_progress, paused, completed, failed]
  last_checkpoint: ISO8601
  error_log: [string]
```

## Error Handling Protocols

### Catastrophic Failure Protocol (Stage 4)

When ALL perspective agents fail:

1. **Detection**: 0 valid outputs from 4 parallel agents
2. **Check for partial outputs**: Any content > 500 chars is "partial"
3. **Present user options**:

```
PERSPECTIVE GENERATION FAILED: All 4 agents failed

Failures:
- Near Field 1: [status]
- Near Field 2: [status]
- Mid Field: [status]
- Far Field: [status]

Options:
(A) Retry SEQUENTIALLY (slower but more stable)
(B) Retry with 2 agents only (Near + Mid fields)
(C) Review partial outputs if available
(D) Abort workflow, preserve problem abstraction
```

### Convergence Floor Protocol

**Minimum Convergence Threshold**:
- 1+ themes appearing in 2+ reports, OR
- 2+ strategies appearing in 2+ reports, OR
- 1+ structural analogies across 2+ domains

**On Zero Convergence**:
1. Re-analyze with relaxed matching (related concepts, not exact)
2. Look for "anti-convergence" (all agree something doesn't work)
3. If still zero: Offer portfolio framing (present as independent options, not validated patterns)

### User Response Timeout Protocol

| Gate | Timeout | Default Action |
|------|---------|----------------|
| Domain Classification | 10 min | Proceed with HIGH confidence domains only |
| Transfer Feasibility | 5 min | Proceed with portfolio view |
| Strategic Priority | 5 min | Present all options without prioritization |

## Parallel Execution

### Stage 4: Perspective Generation (Fork-Join)

```yaml
parallel_execution:
  stage: 4
  pattern: fork_join
  max_agents: 4

  agent_assignment:
    agent_1: "Near Field 1"
    agent_2: "Near Field 2"
    agent_3: "Mid Field"
    agent_4: "Far Field"

  output_paths:
    agent_1: "/tmp/pov-session-{id}/stage-4-perspectives/near-field-1.md"
    agent_2: "/tmp/pov-session-{id}/stage-4-perspectives/near-field-2.md"
    agent_3: "/tmp/pov-session-{id}/stage-4-perspectives/mid-field.md"
    agent_4: "/tmp/pov-session-{id}/stage-4-perspectives/far-field.md"

  completion:
    minimum_required: 2
    proceed_on_timeout: true
    log_missing: true
```

### Stage 5: Verification (Parallel)

```yaml
parallel_execution:
  stage: 5
  pattern: parallel
  max_agents: 2

  assignment:
    agent_1: "Near + Mid perspectives"
    agent_2: "Far perspective + cross-check"
```

## Success Criteria

- [ ] **Diversity**: Perspectives span Near, Mid, and Far fields (4 domains minimum)
- [ ] **Depth**: Structural analogies via structure-mapping principles (>=2 per perspective)
- [ ] **Convergence**: Cross-perspective insights tracked with measurable algorithm
- [ ] **Transferability**: Rigorous transfer evaluation with Far field warnings
- [ ] **Actionability**: Strategic assessment with >=3 recommendations and implementation paths
