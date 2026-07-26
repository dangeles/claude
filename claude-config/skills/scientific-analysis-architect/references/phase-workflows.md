# Phase Workflows

Per-phase detail for scientific-analysis-architect: dispatch briefs, document formats, and gate criteria. The orchestration sequence lives in SKILL.md.

## Contents

- Phase 0: Initialization
- Phase 1: Birds-Eye Planning
- Phase 2: Subsection Planning
- Phase 3: Structure Review
- Phase 4: Plan Review
- Phase 5: Document Generation
- Phase 6: Statistical Fact-Checking
- Phase 7: Audience Document Generation

---

## Phase 0: Initialization

Owner: orchestrator.

Ask for the output directory (default: current working directory) and validate it — exists, writable via a write test; offer alternatives or ask for a new path if either fails. Create the session directory (`{output_dir}/.scientific-analysis-session/`, falling back to `/tmp/scientific-analysis-architect-session-{timestamp}/`) and initialize `session-state.json` with `status: "initialized"`, `current_phase: 0`, and the output directory in config. Check for existing sessions and offer to resume any less than 72 hours old.

**Gate 0**: session directory created and writable, output directory validated, `session-state.json` initialized.

---

## Phase 1: Birds-Eye Planning

Owner: research-architect (Sonnet 4.6).

Ask the user: "Please describe your dataset and research goals. Include: data type (RNA-seq, proteomics, imaging, etc.), sample size and conditions, main research questions."

If the answer carries uncertainty signals ("not sure", "maybe", "could be", "uncertain"), spawn analysis-brainstormer ("Suggest analysis types for: {dataset}") and method-brainstormer ("Suggest methods for: {research_area}") in parallel, aggregate their suggestions, and ask the user which approaches interest them before proceeding.

Then write `research-structure.md` with 3-7 chapters, each with a title, goal, list of analyses, and dependencies on other chapters, and update session state.

```markdown
# Research Structure: {Project Title}

## Overview
{1-2 paragraph summary of research goals and approach}

## Dataset Description
{User-provided description, preserved verbatim}

## Chapters

### Chapter 1: {Chapter Title}
**Goal**: {atlas | hypothesis | mechanism}
**Analyses**:
1. {Analysis 1 name} - {brief description}
2. {Analysis 2 name} - {brief description}
**Dependencies**: None

### Chapter 2: {Chapter Title}
**Goal**: {atlas | hypothesis | mechanism}
**Analyses**:
1. {Analysis 1 name}
2. {Analysis 2 name}
**Dependencies**: Chapter 1 (requires normalized data)

...
```

**Gate 1**: file exists and is valid markdown, 3-7 chapters, each with a goal (atlas/hypothesis/mechanism) and at least one analysis, dependencies reference valid earlier chapters with no cycles.

---

## Phase 2: Subsection Planning

Owner: analysis-planner (Sonnet 4.6).

Read `research-structure.md`, then per chapter fan out to the expert panel in one message so they run concurrently, and aggregate their returns into `chapter{N}-notebook-plans.md`.

**Statistician consultant** — Description: "Statistician consultant: Review statistical approaches for chapter analyses". Prompt: include the chapter analysis requirements; ask for statistical method recommendations, power analysis considerations, multiple comparison handling, and validation approaches. Output to `{session_dir}/consultations/statistician-review.md`.

**Mathematician consultant** — Description: "Mathematician consultant: Review algorithms for chapter analyses". Prompt: include the chapter analysis requirements; ask for algorithm recommendations, complexity analysis, convergence properties, and numerical stability considerations. Output to `{session_dir}/consultations/mathematician-review.md`.

**Programmer consultant** — Description: "Programmer consultant: Review data requirements for chapter analyses". Prompt: include the chapter analysis requirements; ask for data format requirements, library recommendations, performance considerations, and implementation approach. Output to `{session_dir}/consultations/programmer-review.md`.

When consultants disagree, put the conflict to the user with each position rather than picking silently. On a consultant failure, retry once; if the retry fails, escalate for the statistician (critical) and proceed with a logged warning for the mathematician and programmer (optional).

**Gate 2**: all chapters have analysis plans, no unresolved critical conflicts, every analysis has a statistical approach, data flow is consistent (each analysis's inputs come from user data or an upstream output).

---

## Phase 3: Structure Review

Owner: structure-reviewer (Haiku).

Read `research-structure.md` and all chapter plans and check for missing dependencies, redundant analyses, logical flow problems, incomplete specifications, and missing quality controls. Categorize by severity — critical (blocks execution), major (affects quality), minor (improvement opportunity) — and write `structure-review-report.md`.

Present user approval gate 1 with the chapter count, analysis count, issue counts by severity, and the critical issues listed. On "A" record the approval and continue; on "c" ask what should change, route to the phase that owns it, re-run, and return to this gate; on "r" ask whether to return to Phase 1 or Phase 2.

**Gate 3** (user): user reviewed the summary and either approved or requested specific changes, with critical issues addressed.

---

## Phase 4: Plan Review

Owner: notebook-reviewer (Sonnet 4.6), one per chapter in parallel.

Each reviewer gets the chapter plan content and checks pseudocode completeness, statistical correctness, data flow consistency, and edge case coverage. Aggregate into `notebook-review-report.md`, proceeding with available reviews if some fail.

Present user approval gate 2 with per-chapter analysis and issue counts plus the critical issues; handle responses as in Phase 3.

---

## Phase 5: Document Generation

Owner: orchestrator (Step 1) and notebook-generator (Step 2).

Step 1: read `research-structure.md` and every approved `chapter{N}-notebook-plans.md`, synthesize `analysis-strategy-overview.md`, and write it to both `{output_dir}/` and `{session_dir}/`.

Step 2: fan out one generator per chapter in parallel. Prompt: "Generate markdown analysis documents for Chapter {N}: {chapter_plan_content}. Output to: {output_dir}/chapter{N}_{slug}/. Use hybrid prose + fenced pseudocode format. Required sections: Goal, Statistical Approach, Prerequisites, Analysis Steps, Expected Outputs, Notes and Caveats."

Then validate each generated document (required sections present, at least one fenced code block, balanced fences), copy to `{session_dir}/analyses/` as backup, and record `outputs.analyses` and `outputs.strategy_overview` in session state. If some chapters failed, ask the user whether to proceed with available documents, retry the failed chapters, or abort.

### Analysis document format

````markdown
# Analysis Title

<!-- Generated by: scientific-analysis-architect v2.0.0 -->
<!-- Session: {session_id} -->
<!-- Chapter: {N}, Analysis: {M} -->

## Goal
Prose description of what this analysis achieves.

## Statistical Approach
Method, justification, assumptions, and corrections.

## Prerequisites
- Input data and format
- Required libraries
- Upstream dependencies

## Analysis Steps

### Step 1: [Name]
Prose explanation of what this step does and why.

```python
# Pseudocode for step 1
```

### Step 2: [Name]
...

## Expected Outputs
- Output files/objects, format, characteristics

## Notes and Caveats
- Assumptions, limitations, alternatives
````

Code block rules: triple backticks with the `python` language identifier, never nested fences. If pseudocode needs triple-quoted strings, use single-quoted triple quotes in comments; for multi-line string literals use comment notation.

**Gate 5**: all analysis documents created (or partial with user approval), each with the required sections, at least one fenced code block, and balanced fences; master strategy overview present with its required sections; backups in `{session_dir}/analyses/`; file naming follows `analysis{N}_{M}_{slug}.md`.

---

## Phase 6: Statistical Fact-Checking

Owner: statistical-fact-checker (Sonnet 4.6).

Read all generated analysis documents and look for test mismatches, multiple testing issues, interpretation errors, assumption violations, and effect size gaps. Categorize by severity — critical (incorrect conclusions), standard (best practice violation), minor (improvement opportunity).

With zero concerns, show the justification and proceed. With 5 or fewer, go straight to interview mode one concern at a time. With more, present the summary and batch options first. The interview format and batch handling are specified in the interview protocol reference.

Accumulate accepted, rejected, and skipped decisions, then write `statistical-review-report.md` and `corrections-manifest.json`. Ask whether to apply the accepted corrections now or save the manifest for later.

If corrections were applied, regenerate `analysis-strategy-overview.md` from the corrected documents — updating the consolidated methods table and affected chapter summaries — and re-validate it against the Gate 5 overview criteria.

**Gate 6** (user): all concerns reviewed or explicitly skipped, user approved the correction application decision, review report generated, manifest saved, and if corrections were applied the overview refreshed and validated.

---

## Phase 7: Audience Document Generation

Owner: orchestrator.

Confirm the Tier 1 artifacts exist (`analysis-strategy-overview.md`, `research-structure.md`, `session-state.json`), and if `corrections-manifest.json` exists with accepted corrections, that the overview was modified after it. Missing critical artifacts means stopping with a clear error.

Create `{output_dir}/.research-architecture/` if needed, then generate one document at a time:

- **Researcher plan** (Template A) from Tier 1 plus Tier 2 (`chapter{N}-notebook-plans.md`) sources — prose-only, no code blocks, domain language. Write to `{output_dir}/researcher-plan.md`.
- **Architect handoff** (Template B) from Tier 1 plus `session-state.json`, and Tier 4 sources if they exist (review reports, `corrections-manifest.json`) — design rationale, method log, current state, open questions. Write to `{output_dir}/.research-architecture/architect-handoff.md`.
- **Engineering translation** (Template C) from Tier 1 and Tier 2, reading Tier 3 analysis documents per-chapter as needed — system overview, pipeline architecture, data specs, processing stages. Write to `{output_dir}/.research-architecture/engineering-translation.md`.

Copy all three to `{session_dir}/audience-documents/`, then set `current_phase: 7`, add `7` to `completed_phases`, record paths in `outputs.audience_documents`, and set `status: "completed"`.

On resume, check which audience documents already exist, validate them against the Gate 7 section requirements, and regenerate only the missing or invalid ones.

**Gate 7**:

- `{output_dir}/researcher-plan.md` exists and is non-empty (> 500 bytes)
- `{output_dir}/.research-architecture/architect-handoff.md` exists and is non-empty
- `{output_dir}/.research-architecture/engineering-translation.md` exists and is non-empty
- Researcher plan has required sections (case-insensitive): Research Overview, Research Questions, Expected Outcomes, Decision Points
- Architect handoff has required sections: Design Rationale, Current State, Open Questions, Continuation Guidance
- Engineering translation has required sections: System Overview, Pipeline Architecture, Data Specifications, Processing Stages, Resource Requirements, Dependencies
- No fenced code blocks (```) in `researcher-plan.md`
- Backup copies exist in `{session_dir}/audience-documents/`
- All three documents include provenance metadata (HTML comments with session ID)
