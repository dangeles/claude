---
name: lit-pm
last_updated: 2026-06-09
description: Use when coordinating COMPREHENSIVE adaptive literature reviews (4-24h, 9 stages with parallel review discovery, outline synthesis, section writing, fact-checking, editorial polish) with checkpoints calibrated to complexity. Orchestrates literature-researcher, lit-synthesizer, fact-checker, and editor. NOT for quick automated chains (use research-pipeline for fixed 5-stage 2-8h pipeline) or single-paper lookups (use researcher directly).

handoff:
  accepts_from:
    - "*"
  provides_to:
    - "*"
  schema_version: "3.0"
  schema_type: universal
categories:
  - research
  - analysis
---

# lit-pm: Literature Pipeline Manager

## Overview

lit-pm is a Tier 1 orchestrator that coordinates a 9-stage literature review pipeline. It manages parallel review discovery, adaptive checkpoints, and handoffs between specialist skills to produce comprehensive, decision-useful literature reviews.

## Delegation and scope

You coordinate specialists. Delegate specialist work when the subtask is substantial and independent — a search strategy for review discovery, researching and drafting a section, fact-checking a completed section, editorial polish of the full document. Handle it directly when you could finish it in a handful of tool calls, when the work is sequential, or when you need the context in your own loop. See `../references/delegation-and-scope.md`.

Work that stays with you either way:

- Session setup, directory creation, workflow state management
- Outline structure (section headers, thesis statements) — coordination, not prose
- Quality gate evaluation against specialist output
- User communication: summaries, approvals, status
- Complexity detection and checkpoint plan creation

If a required specialist is unavailable, stop and tell the user rather than doing the specialist work yourself.

### Tool selection

| Situation | Tool | Reason |
|-----------|------|--------|
| Specialist doing independent work | Task tool | Separate context, parallel execution |
| 2+ specialists working simultaneously | Task tool (multiple) | Only way to parallelize |
| Loading domain knowledge for your own decisions | Skill tool | Shared context needed |

### Progress tracking

Workflow state carries `stage_current` and `stage_completed`. Lead status updates with the stage you are in — "[Stage 5/8 - Parallel Section Research]" — and after a user interaction, say which stage you are returning to. While parallel agents run, keep a status board (agent, task, running/complete/failed).

### Three-tier architecture

- **Tier 1 — orchestrator (this skill)**: 9-stage pipeline, complexity detection, parallel execution with convergence tracking, workflow state, handoffs, quality gates, session storage.
- **Tier 2 — literature specialists**: `literature-researcher` (review discovery, section research); `lit-synthesizer` (narrative synthesis, introduction, conclusion).
- **Tier 3 — supporting**: `fact-checker` (quick validation + comprehensive review), `editor` (final polish), `requirements-analyst` (scope refinement), `devils-advocate` (adversarial review).

## When to Use This Skill

- **Internal research synthesis**: Decision-focused ("Should we pursue technology X?")
- **Literature surveys**: Landscape mapping ("What are current approaches to Y?")
- **Comprehensive reviews**: Grant proposals, research plans, scientific documents
- **Cross-domain literature**: Topics spanning multiple fields
- **High-stakes deliverables**: When thoroughness and accuracy matter

## When NOT to Use This Skill

- **Quick literature lookups**: Use researcher directly for single-topic searches
- **Single-paper analysis**: Use researcher for individual paper deep-dives
- **Non-scientific literature**: This skill is optimized for scientific/technical literature
- **Time-critical requests**: Pipeline takes 4-24 hours; use researcher for faster turnaround

---

## The 9-stage pipeline

Each stage is a distinct deliverable with a named owner. How a specialist reaches its deliverable is up to that specialist. Per-stage detail lives in `references/stage-specifications.md`.

### Stage 0: Session setup and archival guidelines

Owner: lit-pm. Runs automatically.

Create the session directory `/tmp/lit-pm-session-$(date +%Y%m%d-%H%M%S)-$$/` and store its path in workflow state. Then extract archival guidelines, preferring `.archive-metadata.yaml` at the repo root over `CLAUDE.md`:

1. Read `~/.claude/skills/archive-workflow/references/archival-compliance-check.md` and apply its 5-step pattern to file creation. If that file is missing, log a warning and proceed without the archival check.
2. Read `.archive-metadata.yaml` for naming conventions, directory structure, and project type.
3. If `.archive-metadata.yaml` is absent, fall back to `CLAUDE.md` prose extraction (repository organization, naming conventions, git rules, document structure, PDF paths) and warn: "Archival guidelines read from CLAUDE.md (fallback). Run archive-workflow to generate .archive-metadata.yaml." If `docs/organization/final-organization-report.md` exists, also warn that archival metadata was previously present and is now missing.
4. Write `archival-guidelines-summary.md` to the session directory recording source, project type, naming conventions, directory structure, and enforcement mode (default "advisory").

```yaml
session_setup:
  session_dir: "/tmp/lit-pm-session-{timestamp}-{pid}/"
  archival_summary_path: "{session_dir}/archival-guidelines-summary.md"
  guidelines_found: boolean
  guidelines_source: string  # ".archive-metadata.yaml" or "CLAUDE.md" or "defaults"
  enforcement_mode: string   # "advisory" | "soft-mandatory" | "hard-mandatory"
```

Include this block in every specialist dispatch (the standard archival context block from archival-compliance-check.md):

```yaml
archival_context:
  guidelines_present: true/false
  source: ".archive-metadata.yaml"  # or "CLAUDE.md" or "defaults"
  naming_convention: "kebab-case"
  output_directory: "docs/literature/{topic}/"
  enforcement_mode: "advisory"
  user_override: null
```

Malformed `.archive-metadata.yaml`: treat as missing. No `CLAUDE.md`: use defaults and log. Session directory creation failure is fatal — the pipeline needs session isolation, so stop there.

Session cleanup: delete the session directory once Stage 8 completes; on failure or abort, keep it and log the path for debugging.

### Stage 1: Scope refinement

Owner: requirements-analyst. Checkpoint: always.

Clarify the research question, define success criteria, set boundaries. Then classify the complexity tier and generate the matching depth_profile using the Depth Profile System table in `references/adaptive-orchestration.md`. Store it in `workflow_state.depth_profile`.

Present the checkpoint plan and depth profile together for approval using the combined format in `references/adaptive-orchestration.md` (User Presentation Format), and apply any user overrides.

Gate: specific research question, measurable criteria, clear boundaries, user approval.

### Stage 2: Parallel review discovery

Owner: 2-3 literature-researcher agents via Task tool. Checkpoint: high-stakes only.

Launch agents with diverse search strategies, then analyze convergence — reviews found independently by multiple agents are the strongest signal. Target 6-9 reviews.

Gate: 6-9 reviews, at least 2 showing convergence, major themes covered.

### Stage 3: Layered outline synthesis

Owner: lit-pm (structure) + literature-researcher (section detail). Checkpoint: medium/complex/high-stakes.

You write the outline structure — `depth_profile.sections.target_count` sections, each with a specific thesis. Delegate subsection proposals and review assignments per section.

Gate: balanced sections at the target count, specific theses, user approval where the plan calls for a checkpoint.

### Stage 4: Introduction writing

Owner: lit-synthesizer, then editor for a quick polish. Runs automatically.

The introduction frames the research question and previews the structure. Include the full depth_profile YAML block under a "DEPTH PROFILE:" header in the lit-synthesizer prompt.

Gate: clear framing, structure preview consistent with the outline.

### Stage 5: Parallel section research and writing

Owner: literature-researcher agents in parallel (max 3 concurrent). Gated by Stage 6a rather than a checkpoint.

Section writers are researchers, not review summarizers: targeted primary-literature search against their section thesis, per `depth_profile.research.papers_per_section`, plus a recency survey covering the last 6-12 months. Writers may add subsections; the assigned thesis stays fixed.

Send each section to Stage 6a as it completes rather than waiting for the whole set. If a section stalls out, proceed without it and flag the gap.

Gate: papers cited within the depth profile range, recency survey present, thesis addressed, no placeholder text.

### Stage 6a: Per-section quick validation (blocking)

Owner: fact-checker.

Mechanical checks: paper count, recency survey presence, no placeholders, thesis addressed. A section reaches synthesis only after it passes. Cap at 3 revision cycles; on the third failure, present the user with options (accept as-is, adjust requirements, reassign, drop the section) and record the decision in `workflow_state.quality_overrides`.

### Stage 6b: Comprehensive fact-check (non-blocking)

Owner: fact-checker.

Cross-section consistency, citation spot-checks against real sources, quantitative verification. Produces a P0/P1/P2 revision list that Stage 8 consumes.

### Stage 6c: Devil's advocate section review

Owner: devils-advocate. Always runs as a quality gate, not a user approval. Triggers once sections pass 6a and 6b.

Adversarial review per section: argument quality, assumptions, logical gaps. Up to 2 exchanges per section; if unresolved after that, pass with an uncertainty note.

Scope boundary against fact-checker — devils-advocate challenges argument strength, assumption validity, logical coherence, thesis appropriateness, and methodology context for claims. Citation accuracy, whether papers exist, and whether values match sources belong to fact-checker.

### Stage 7: Active synthesis and augmentation

Owner: lit-synthesizer as senior author. Checkpoint: high-stakes only.

Read all sections, identify cross-cutting themes, restructure for narrative flow, write the conclusion. The synthesizer has authority to add subsections and rewrite transitions, and flags additions over 20%. Include the full depth_profile block under a "DEPTH PROFILE:" header; it governs `synthesis.augmentation_budget` and `synthesis.conclusion_scope`.

Gate: logical flow, themes identified, gaps filled, conclusion synthesizes findings.

After Stage 7: if synthesis added >=20% content or complexity is high-stakes, run Stage 7.5; otherwise go to Stage 8.

### Stage 7.5: Devil's advocate synthesis review (conditional)

Owner: devils-advocate. Triggers on `addition_percentage >= 20%` or high-stakes complexity.

Strategic review of the synthesized document: thesis coherence across sections, cross-cutting theme validity, argument flow. Up to 2 exchanges, then proceed with documented uncertainty.

### Stage 8: Editorial polish

Owner: editor. Runs automatically.

Incorporate P0/P1 revisions from Stage 6b, polish for clarity, keep voice consistent, finalize formatting. Include the full depth_profile block under a "DEPTH PROFILE:" header; it governs `writing.density` and `writing.density_guidance`.

Gate: revisions incorporated, consistent voice, formatted.

### Post-pipeline: git strategy advisory (optional)

Once the final document is ready you can invoke `git-strategy-advisor` via Task tool in post-work mode:

```
Use git-strategy-advisor to determine git strategy for completed work.

mode: post-work
```

Read the advisor's `summary` field and include it in the completion summary. Skip the section silently if confidence is "none" or "low", if the advisor is unavailable or errors, or if the output files live outside the current git repository (the advisor only sees changes inside it). lit-pm has no built-in git logic of its own.

---

## Adaptive orchestration

Full complexity detection logic: `references/adaptive-orchestration.md`.

Complexity comes from three dimensions: **scope** (paper count — <10 simple, 10-30 medium, 30+ complex — plus topic breadth and literature maturity), **stakes** (from "quick survey" through "grant proposal"), and **user hints** (`--review-outline`, `--full-auto`, time constraints).

| Complexity | Stage 1 | Stage 2 | Stage 3 | Stage 6c | Stage 7 | Stage 7.5 | Rationale |
|------------|---------|---------|---------|----------|---------|-----------|-----------|
| Simple | CHECKPOINT | Auto | Auto | ACTIVE | Auto | Conditional* | Scope approval sufficient |
| Medium | CHECKPOINT | Auto | CHECKPOINT | ACTIVE | Auto | Conditional* | Direction check before heavy lifting |
| Complex | CHECKPOINT | Auto | CHECKPOINT | ACTIVE | CHECKPOINT | Conditional* | Multiple approval points |
| High-Stakes | CHECKPOINT | CHECKPOINT | CHECKPOINT | ACTIVE | CHECKPOINT | ACTIVE | Maximum oversight |

- **CHECKPOINT**: user approval before proceeding
- **Auto**: runs without user interaction
- **ACTIVE**: always runs as a quality gate, not a user approval
- **Conditional\***: triggers when synthesis adds >=20% content

After proposing the plan the user can accept it, skip specific checkpoints, add checkpoints, or pass `--full-auto` to skip all optional ones.

---

## Resource limits

```yaml
resource_limits:
  max_concurrent_agents: 4       # Hard ceiling
  max_parallel_researchers: 3    # Stage 2: Leave slot for orchestrator
  max_parallel_sections: 3       # Stage 5: Leave slot for fact-checker
  queue_behavior: FIFO           # When limits reached
```

## Parallel execution

Fan-out guidance: `../references/delegation-and-scope.md`.

**Stage 2 — review discovery.** Launch 2-3 literature-researcher agents in one message so they run concurrently, each with a different angle: broad keywords, specific technical terms, application focus. Each finds 3-5 review papers and writes to `{session_dir}/reviews/agent-{n}-{strategy}.yaml`. When all return, score convergence: found by 3/3 agents is must-read, 2/3 is medium priority, 1 agent contributes coverage breadth.

**Stage 5 — section writing.** Launch literature-researcher per section, max 3 concurrent. Each prompt carries the section thesis, assigned reviews, subsection structure, and the full depth_profile block under a "DEPTH PROFILE:" header, and writes to `{session_dir}/sections/{section_id}.md`. Route each finished section straight to Stage 6a.

---

## Handoffs

Full YAML schema: `references/handoff-schema.md`.

```yaml
handoff:
  version: "1.1"
  stage: integer           # 0-8
  status: enum             # pending | in_progress | complete | failed
  producer: string         # skill that produced handoff
  consumer: string         # skill that receives handoff
  workflow_id: string      # unique identifier
  timestamp: ISO8601       # when created
  session:                 # Added in v1.1
    session_dir: string    # Path to /tmp/lit-pm-session-{...}/
    archival_guidelines_path: string  # Path to archival-guidelines-summary.md
```

Each stage adds its own fields on top of this base. Every downstream agent also receives session context:

```yaml
session_context:
  archival_guidelines_path: "{session_dir}/archival-guidelines-summary.md"
  output_directory: string   # Where final document should be written
  pdf_storage_path: string   # Where to store downloaded PDFs
  naming_convention: object  # File prefix rules
```

`literature-researcher` uses the naming conventions for paper notes and PDF paths, `lit-synthesizer` the document structure and citation format, `editor` the writing style and formatting rules, `fact-checker` the citation format.

---

## Quality gates

Per-stage criteria: `references/quality-gates.md`.

Automated gates are programmatic: paper count threshold (>=15 per section by default), recency survey presence, placeholder detection ("TODO", "[CITE]"), word count range. Judgment gates need agent assessment: thesis specificity, narrative flow, cross-cutting theme identification.

Quality floor — these hold even under `--full-auto`:

- Stage 2: minimum 4 reviews, minimum 1 convergence
- Stage 3: minimum 2 sections, no section above 50% of total
- Stage 6c: devils-advocate runs and identifies a thesis for each section

Quality floor violations escalate to the user regardless of the checkpoint plan.

---

## Error handling and resumption

Compensation matrix and recovery detail: `references/error-handling.md`.

Key mechanisms: per-stage compensation actions, proceeding with partial results after repeated failures, atomic state writes (`.tmp` -> validate -> rename -> `.bak`), and graceful interrupt handling that preserves partial state.

```yaml
workflow_state:
  workflow_id: "lit-review-{topic}-{date}"
  stage_current: integer          # 0-8
  stage_completed: [list]
  checkpoints_remaining: [list]
  session:                        # Added for session management
    session_dir: string           # /tmp/lit-pm-session-{...}/
    archival_guidelines_path: string
    cleanup_on_complete: boolean  # Default true
  depth_profile: {}  # Full depth_profile block from adaptive-orchestration.md Depth Profile System
  artifacts:
    scope: "/path/to/scope.yaml"
    reviews: "/path/to/reviews.yaml"
    outline: "/path/to/outline.yaml"
    sections_complete: [list]
    sections_in_progress: [list]
```

Resume with `lit-pm --resume workflow-id`. Reuse the session directory if it still exists; if it is gone, re-run Stage 0 to recreate it (non-destructive).

---

## References

- `references/stage-specifications.md`: per-stage inputs, outputs, and gate criteria
- `references/handoff-schema.md`: YAML schema for stage-to-stage communication
- `references/error-handling.md`: compensation logic, recovery
- `references/adaptive-orchestration.md`: complexity detection, checkpoint plans, depth profiles
- `references/quality-gates.md`: per-stage quality validation criteria
- `../references/delegation-and-scope.md`: shared orchestrator conventions

## Examples

- `examples/hepatocyte-review-example.md`: complete walkthrough of a hepatocyte oxygenation review

---

## Invocation

```bash
# Standard invocation
lit-pm "Comprehensive review of [topic] for [purpose]"

# With explicit checkpoint control
lit-pm --review-outline "Survey of current approaches to [topic]"
lit-pm --full-auto "Quick survey of [topic]"

# Resume interrupted workflow
lit-pm --resume workflow-id
```

## Dependencies

- technical-pm parallel execution capabilities
- literature-researcher, lit-synthesizer, fact-checker, editor, requirements-analyst, devils-advocate
