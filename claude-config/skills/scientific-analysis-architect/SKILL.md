---
name: scientific-analysis-architect
description: Use when planning multi-chapter scientific research analyses with expert consultation — produces markdown analysis documents with pseudocode for RNA-seq, proteomics, or other data analysis workflows. NOT for software architecture (use systems-architect) or single-component design within an existing analysis (use senior-developer).
version: 2.0.0
last_updated: 2026-05-24
tags: [scientific-analysis, multi-agent, markdown, research-planning, pseudocode]
---

# scientific-analysis-architect

Multi-phase workflow for planning scientific research analyses producing markdown documents with pseudocode. Biology-agnostic design: agents request context via user prompts and never inject biological interpretation.

## Delegation and scope

You are the architect: you plan how specialists work together. Delegate specialist work when the subtask is substantial and independent — statistical design for a chapter, algorithm requirements, generating a chapter's analysis documents, statistical fact-checking of the finished set. Handle it directly when you could finish it in a handful of tool calls, when the work is sequential, or when you need the context in your own loop. See `../references/delegation-and-scope.md`.

Work that stays with you either way:

- Session setup, directory creation, state file management
- Quality gate evaluation and validation commands (markdown structure checks, dependency verification)
- User communication: summaries, approvals, status
- Workflow coordination: reading state, tracking progress, managing handoffs

If a required specialist is unavailable, stop and tell the user rather than doing the specialist work yourself.

### Tool selection

| Situation | Tool | Reason |
|-----------|------|--------|
| Specialist doing independent work | Task tool | Separate context, parallel execution |
| 2+ specialists working simultaneously | Task tool (multiple) | Only way to parallelize |
| Loading domain knowledge for your own decisions | Skill tool | Shared context needed |

### Progress tracking

`{session_dir}/session-state.json` carries `current_phase` and `completed_phases`. Lead status updates with the phase you are in — "[Phase 4/7 - Plan Review]" — and after a user interaction, say which phase you are returning to.

## When to Use

- Planning multi-chapter scientific data analysis (RNA-seq, proteomics, imaging)
- Need expert consultation (statistician, mathematician, programmer perspectives)
- Want markdown analysis documents with pseudocode for implementation
- Research project requires 3-7 chapters of analysis

## When NOT to Use

- Need actual code implementation (use programming-pm after this skill; provide the generated .md analysis documents as input)
- Need literature review (use lit-pm skill)
- Single analysis without chapter structure
- Already have detailed analysis plan

## Workflow Overview

| Phase | Owner | Deliverable |
|-------|-------|-------------|
| 0: Initialization | orchestrator | Session directory, validated output directory, `session-state.json` |
| 1: Birds-eye planning | research-architect (Sonnet 4.6) | `research-structure.md` (3-7 chapters) |
| 2: Subsection planning | analysis-planner (Sonnet 4.6) + 3 consultants | `chapter{N}-notebook-plans.md` |
| 3: Structure review | structure-reviewer (Haiku) | `structure-review-report.md` -> user approval gate 1 |
| 4: Plan review | notebook-reviewer (Sonnet 4.6), parallel | `notebook-review-report.md` -> user approval gate 2 |
| 5: Document generation | orchestrator + notebook-generator | `analysis-strategy-overview.md` + `.md` analysis documents with pseudocode |
| 6: Statistical fact-checking | statistical-fact-checker (Sonnet 4.6) | Corrected analysis documents, refreshed overview |
| 7: Audience documents | orchestrator | `researcher-plan.md` + `.research-architecture/{architect-handoff,engineering-translation}.md` |

**Estimated Runtime**: 61-81 minutes for 3 chapters

## Phase 0: Initialization

Create the session directory — primary `{output_directory}/.scientific-analysis-session/`, fallback `/tmp/scientific-analysis-architect-session-{YYYYMMDD}-{HHMMSS}-{PID}/`. Check the output directory exists and is writable with a write test, and offer alternatives if it is not. Initialize `session-state.json` with status "initialized".

Then read `~/.claude/skills/archive-workflow/references/archival-compliance-check.md` and apply its 5-step pattern to file creation; if that file is missing, log a warning and proceed without the archival check. Store the guidelines in `session-state.json`, validate proposed analysis document and directory paths against them, and pass `archival_context` to all downstream agent dispatches.

**Quality Gate 0**: Session directory created, output directory validated.

## Phase 1: Birds-Eye Planning

Owner: research-architect (Sonnet 4.6).

Ask the user to describe their dataset and research goals. If the answer signals uncertainty ("not sure", "maybe"), fan out to analysis-brainstormer and method-brainstormer (Haiku 4.5) and present their suggestions before continuing. Then generate `{session_dir}/research-structure.md` with 3-7 chapters, each carrying a goal and its analyses.

Biology-agnostic behavior: agents ask "What biological questions are you trying to answer?" rather than asserting "You should look at cell types".

**Quality Gate 1**: Structure has 3-7 chapters, each with goal and analyses.

## Phase 2: Subsection Planning

Owner: analysis-planner (Sonnet 4.6).

For each chapter, fan out to the expert panel in parallel (all Haiku): statistician-consultant for statistical approach validation, mathematician-consultant for algorithm requirements, programmer-consultant for data requirements. Aggregate their recommendations into `{session_dir}/chapter{N}-notebook-plans.md`, and put any consultant disagreement to the user rather than picking silently.

On parallel failures, wait for all retries and escalate once with all failures in a single prompt. The statistician is critical; the other two are optional.

**Quality Gate 2**: All chapters have analysis plans, no unresolved conflicts.

## Phase 3: Structure Review

Owner: structure-reviewer (Haiku).

Review `research-structure.md` and all chapter plans for missing dependencies, redundancies, and logical issues, and write `{session_dir}/structure-review-report.md`.

**USER APPROVAL GATE 1**:
```
Structure Review Complete

Summary:
- {N} chapters planned
- {M} analyses total
- {K} issues identified

Approve / Request changes / Reject? [A/c/r]
```

## Phase 4: Plan Review

Owner: notebook-reviewer (Sonnet 4.6), one per chapter in parallel.

Reviewers check pseudocode completeness, statistical correctness, and data flow. Aggregate into `{session_dir}/notebook-review-report.md`.

**USER APPROVAL GATE 2**:
```
Plan Review Complete

Per-Chapter Summary:
- Chapter 1: {N} analyses, {K} issues
...

Approve / Request changes / Reject? [A/c/r]
```

## Phase 5: Document Generation

### Step 1: Master strategy overview (orchestrator)

Read `research-structure.md` and every `chapter{N}-notebook-plans.md`, then write `analysis-strategy-overview.md` to both `{output_dir}/` and `{session_dir}/`. It contains: project objective, dataset summary, Strategy at a Glance table (chapter, title, goal, analyses, key method), chapter summaries, data flow between chapters (text diagram), consolidated methods table with justification, required libraries, execution order with dependency notes, and assumptions and limitations.

This document synthesizes already-approved content; it does not introduce new analyses or methods. Keep it to 1-3 pages with 2-4 sentence chapter summaries (what and why, not how). For <= 4 chapters, omit Execution Order if the flow is linear; for >= 6 chapters, include a dependency graph.

### Step 2: Analysis documents (notebook-generator)

Fan out one generator per chapter in parallel. Each produces `.md` files mixing prose with fenced pseudocode, written to both the output directory and the session directory as backup. Every analysis document uses these section headings, which the Gate 5 validator matches on:

- `## Goal`: what this analysis achieves
- `## Statistical Approach`: method, justification, assumptions, corrections
- `## Prerequisites`: input data, required libraries, upstream dependencies
- `## Analysis Steps`: numbered steps with prose plus fenced Python pseudocode blocks
- `## Expected Outputs`: output files/objects, format, characteristics
- `## Notes and Caveats`: assumptions, limitations, alternatives

Code block rules: triple backticks with the `python` language identifier, never nested fences. If pseudocode needs triple-quoted strings, use single-quoted triple quotes inside comments; for multi-line string literals use comment notation instead.

If some chapters fail, offer to proceed with what is available and allow per-chapter regeneration later.

**Output**:
- `{output_dir}/chapter{N}_{slug}/analysis{N}_{M}_{slug}.md`
- `{output_dir}/analysis-strategy-overview.md`
- `{session_dir}/analyses/` and `{session_dir}/analysis-strategy-overview.md` (backups)

**Quality Gate 5**: Every analysis document has the required sections (Goal, Statistical Approach, Analysis Steps, Expected Outputs), at least one fenced code block, and balanced fences. Master strategy overview exists with required sections. Run the validator in `references/quality-gates.md`.

## Phase 6: Statistical Fact-Checking

Owner: statistical-fact-checker (Sonnet 4.6). Full protocol: `references/interview-protocol.md`.

**INTERVIEW MODE**: with <= 5 concerns, present them one at a time; with more, present a summary first and offer batch options.

```
Statistical Concern {N} of {total}

Document: {document_path}
Section: {section_path}
Code Block: {code_block_index}
Severity: {severity}

Issue: {description}

Current: {current_content}
Recommendation: {recommended_fix}

Accept? [yes/no/skip/explain]
```

Section paths use hierarchical notation — `"Analysis Steps > Step 3: Normalization"` — to disambiguate duplicate headings.

Batch options after 5 concerns: continue one-by-one, accept all remaining, reject all remaining, or accept critical/standard and skip minor. Then:

```
Summary:
- {X} accepted, {Y} rejected, {Z} skipped

Apply corrections? [yes/no]
```

If corrections were applied, regenerate `analysis-strategy-overview.md` so the methods and approaches match the corrected documents, and re-validate it against Gate 5.

**Output**:
- `{session_dir}/statistical-review-report.md`
- `{session_dir}/corrections-manifest.json`
- Updated `.md` analysis documents and refreshed `analysis-strategy-overview.md` (if corrections applied)

### Post-workflow: git strategy advisory (optional)

Once analysis documents are finalized you can invoke `git-strategy-advisor` via Task tool in post-work mode:

```
Use git-strategy-advisor to determine git strategy for completed work.

mode: post-work
```

Read the advisor's `summary` field into the completion summary. Skip the section silently if confidence is "none" or "low", if the advisor is unavailable or errors, or if output files were written outside the current git repository (the advisor only sees changes inside it). This skill has no built-in git logic of its own.

## Phase 7: Audience Document Generation

Owner: orchestrator.

Three audience-targeted documents, all synthesized from already-approved and fact-checked material — no new analyses or methods.

**Input artifacts**, tiered:
- Tier 1 (always read): `analysis-strategy-overview.md`, `research-structure.md`, `session-state.json`
- Tier 2 (per-chapter): `chapter{N}-notebook-plans.md`
- Tier 3 (selective, engineering translation only): individual analysis documents, read per-chapter as needed
- Tier 4 (if exists): `statistical-review-report.md`, `corrections-manifest.json`, review reports

Generate documents one at a time. The researcher plan and architect handoff draw on Tier 1 and Tier 2; the engineering translation reads Tier 3 per-chapter rather than loading everything at once.

Before generating, confirm the Tier 1 artifacts exist and are non-empty, and if `corrections-manifest.json` exists with accepted corrections, that `analysis-strategy-overview.md` was modified after it. Missing critical artifacts means stopping with a clear error.

Then create `{output_dir}/.research-architecture/` if absent and write:

- `{output_dir}/researcher-plan.md` — see [audience-document-templates.md](references/audience-document-templates.md) Template A
- `{output_dir}/.research-architecture/architect-handoff.md` — Template B
- `{output_dir}/.research-architecture/engineering-translation.md` — Template C

Copy all three to `{session_dir}/audience-documents/`, then update session state: `current_phase: 7`, add `7` to `completed_phases`, record paths in `outputs.audience_documents`, set `status: "completed"`.

On resume, skip regenerating any audience document that already exists and passes section validation.

**Quality Gate 7**: All 3 audience documents exist with their required sections, and backups exist. Validator: [quality-gates.md](references/quality-gates.md).

**Completion Announcement**:
```
Audience Documents Generated

Three audience-targeted documents have been created:

1. Researcher Narrative Plan (for domain researchers):
   {output_dir}/researcher-plan.md

2. Architect Handoff (for analysis architects):
   {output_dir}/.research-architecture/architect-handoff.md

3. Engineering Translation (for systems engineers):
   {output_dir}/.research-architecture/engineering-translation.md

Backup copies saved to: {session_dir}/audience-documents/

Workflow complete.
```

## Session Management

### Session Directory Structure

```
{session_dir}/
+-- session-state.json          # Resumable state
+-- research-structure.md       # Phase 1 output
+-- chapter1-notebook-plans.md  # Phase 2 output
+-- chapter2-notebook-plans.md
+-- structure-review-report.md  # Phase 3 output
+-- notebook-review-report.md   # Phase 4 output
+-- analysis-strategy-overview.md # Phase 5 output (master overview)
+-- statistical-review-report.md # Phase 6 output
+-- corrections-manifest.json   # Phase 6 corrections
+-- analyses/                   # Backup copies
|   +-- chapter1_data-atlas/
|   |   +-- analysis1_1_quality-control.md
|   |   +-- analysis1_2_normalization.md
|   +-- chapter2_hypothesis-testing/
|       +-- analysis2_1_differential-expression.md
+-- audience-documents/            # Phase 7 output (backup)
|   +-- researcher-plan.md
|   +-- architect-handoff.md
|   +-- engineering-translation.md
+-- logs/
    +-- workflow.log
```

### Resume Protocol

On invocation, check for existing sessions (output_dir first, then /tmp). If one is found and is less than 72 hours old:

```
Found incomplete session from {timestamp}
Project: {research_goals}
Status: Phase {N}

Resume? [yes/no]
```

On yes, load state and continue from the current phase; on no, archive the old session and start fresh. On Ctrl+C, save state with status "interrupted" and print "Session saved. Resume with: /scientific-analysis-architect".

## Error Handling

Full specification: [error-handling.md](references/error-handling.md).

Retry a failed agent once. If it fails again, ask the user whether to proceed without it or abort — with the caveat that Phase 2's statistician is critical while the mathematician and programmer consultants are optional. When a phase stalls, prefer proceeding with the partial results you have and telling the user what is missing over blocking: Phase 2 can proceed with available consultants, Phase 4 with available reviews, Phase 5 with partial generation plus a retry offer, Phase 7 with available documents plus a warning. A Phase 0 failure is fatal.

## Quality Gates Summary

| Gate | Phase | Owner | Pass Criteria |
|------|-------|-------|---------------|
| 0 | 0 | Orchestrator | Session created, output directory validated |
| 1 | 1 | research-architect | 3-7 chapters with goals |
| 2 | 2 | analysis-planner | All chapter plans, no critical conflicts |
| 3 | 3 | User | Approve structure |
| 4 | 4 | User | Approve analysis plans |
| 5 | 5 | Orchestrator + notebook-generator | Valid .md analysis documents (structure validated), master overview present |
| 6 | 6 | User | Interview complete, corrections applied |
| 7 | 7 | Orchestrator | All 3 audience documents exist, required sections present, backups exist |

## Dependencies

- **Tools**: Task, AskUserQuestion, Read, Write, Bash
- **Python Packages**: None required
- **Complements**: lit-pm (literature), programming-pm (implementation -- accepts markdown analysis documents as input)
- **Output**: Pseudocode analysis documents (.md) for manual or programming-pm implementation

## References

- [agent-definitions.md](references/agent-definitions.md)
- [phase-workflows.md](references/phase-workflows.md)
- [interview-protocol.md](references/interview-protocol.md)
- [notebook-templates.md](references/notebook-templates.md)
- [audience-document-templates.md](references/audience-document-templates.md)
- [session-schema.md](references/session-schema.md)
- [error-handling.md](references/error-handling.md)
- [quality-gates.md](references/quality-gates.md)
- [../references/delegation-and-scope.md](../references/delegation-and-scope.md)

## Examples

- [rnaseq-analysis-plan.md](examples/rnaseq-analysis-plan.md)
- [statistical-interview-session.md](examples/statistical-interview-session.md)
