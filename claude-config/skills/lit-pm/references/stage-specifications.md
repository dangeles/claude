# Stage Specifications

Inputs, deliverables, output schemas, and gate criteria for the 9-stage pipeline (Stage 0-8). The orchestration logic lives in SKILL.md; this file is the per-stage detail.

## Contents

- Stage 0: Archival Guidelines Review
- Stage 1: Scope Refinement
- Stage 2: Parallel Review Discovery
- Stage 3: Layered Outline Synthesis
- Stage 4: Introduction Writing
- Stage 5: Parallel Section Research & Writing
- Stage 6a: Per-Section Quick Validation
- Stage 6b: Comprehensive Final Fact-Check
- Stage 7: Active Synthesis & Augmentation
- Stage 8: Editorial Polish
- Stage Summary Table
- Session Management Summary

---

## Stage 0: Archival Guidelines Review

Owner: lit-pm, automatic.

Initialize the workflow session and extract project-specific archival guidelines so downstream agents follow consistent naming, directory structure, and formatting.

```bash
SESSION_DIR="/tmp/lit-pm-session-$(date +%Y%m%d-%H%M%S)-$$"
mkdir -p "$SESSION_DIR"
```

Locate the project `CLAUDE.md` (working directory, then up to 3 parent levels; defaults if absent) and extract repository organization, naming conventions, version control rules, document structure requirements, citation format, and PDF acquisition paths. Write the summary to `{SESSION_DIR}/archival-guidelines-summary.md` and record the session in workflow state.

```yaml
stage_0_output:
  session:
    session_dir: string       # Full path to session directory
    created_at: ISO8601       # When session was created
    archival_summary_path: string  # Path to guidelines summary
  guidelines:
    source: string            # Path to CLAUDE.md or "defaults"
    found: boolean            # Whether CLAUDE.md was found
    directory_structure:
      literature: string      # e.g., "docs/literature/<topic>/"
      pdf_storage: string     # e.g., "docs/literature/<topic>/pdfs/"
      reports: string         # e.g., "docs/reports/"
    naming_conventions:
      review: string          # e.g., "review-"
      analysis: string        # e.g., "analysis-"
      reference: string       # e.g., "reference-"
      paper_notes: string     # e.g., "<author>-<year>-"
    document_structure: list  # Required sections
    citation_format: string   # e.g., "Nature-style inline"
    git_rules:
      commit_after_edit: boolean
      no_version_numbers: boolean
```

### Archival summary contents

`archival-guidelines-summary.md` records, for the session:

- **Directory structure**: literature reviews in `docs/literature/<topic>/`, PDFs in `docs/literature/<topic>/pdfs/`, reports and analyses in `docs/reports/`, paper notes alongside the review.
- **Naming conventions**: `review-` for literature reviews, `analysis-` for analyses, `reference-` for reference documents, `<author>-<year>-` for paper notes.
- **Document structure**: title; metadata (version, date, sources); executive summary (1-3 paragraphs); table of contents for 3+ sections; body numbered hierarchically; key parameters table with sources; gaps and limitations; references with DOIs; revision history.
- **Citation format**: Nature-style inline. Superscript numbers in text ("as shown previously¹,²"), bracketed numbers in tables ("[1]"), bibliography as Author(s). Title. *Journal* **volume**, pages (year). DOI.
- **Git rules**: commit after every edit to `docs/` or `modules/`, no version-numbered files (use git history), edit in place.
- **Restriction**: never download patent documents (PDFs, applications, or patent-related files).

### Gate and failure handling

Gate: session directory created and writable, archival summary written, session path in workflow state.

| Condition | Action |
|-----------|--------|
| Session directory creation fails | Stop — the pipeline needs session isolation |
| CLAUDE.md not found | Use defaults, log warning, continue |
| Parsing error in CLAUDE.md | Partial extraction plus defaults, log warning |

Cleanup: `rm -rf "$SESSION_DIR"` once Stage 8 completes. On failure or abort, keep the directory and log "Session preserved at: {SESSION_DIR}".

---

## Stage 1: Scope Refinement

Owner: requirements-analyst. Checkpoint: always.

Input: the user's initial prompt. Clarify the research question, audience, and intended use; define success criteria; set in-scope and out-of-scope boundaries; classify complexity; propose the checkpoint plan and depth profile for approval.

```yaml
scope_document:
  research_question: string  # Specific, testable
  success_criteria: list     # Measurable outcomes
  in_scope: list            # Topics to cover
  out_of_scope: list        # Explicit exclusions
  complexity_tier: enum     # Simple | Medium | Complex | High-Stakes
  checkpoint_plan: object   # Which stages have checkpoints
```

Gate: research question is specific (not "what is X?" but "what governs X?"), success criteria measurable, boundaries clear, user approves scope and checkpoint plan.

If the user rejects the scope, requirements-analyst revises from their feedback.

---

## Stage 2: Parallel Review Discovery

Owner: lit-pm orchestrating 2-3 literature-researcher agents in parallel. Checkpoint: high-stakes only.

Input: the refined scope document. Generate 2-3 diverse search strategies — broad keywords ("CAR-T manufacturing review"), specific technical terms ("transduction efficiency optimization"), application focus ("clinical CAR-T production challenges") — and give one to each agent. Each finds 2-3 high-quality reviews with a brief annotation covering key findings, coverage, and quality.

### Convergence tracking

Collect reviews from all agents with source attribution, then match duplicates in this order: DOI exact match (highest confidence), PubMed ID exact match, title fuzzy match (Levenshtein distance < 0.2), first author + year + journal. Score by how many agents found each review — 3/3 is high priority (must-read), 2/3 medium, 1/3 a unique perspective contributing coverage breadth.

```yaml
review_collection:
  reviews:
    - title: string
      doi: string | null
      authors: string
      year: integer
      source_agent: string
      convergence_count: integer
      annotation: string
  convergence_analysis:
    high_priority: list
    medium_priority: list
    unique_perspectives: list
    themes_covered: list
    gaps_identified: list
```

Gate: 6-9 reviews total, at least 2 showing convergence, published in the last 10 years unless older work is definitive, all major themes in scope covered.

Recovery: convergence below 2 — adjust search strategies and re-run. Fewer than 4 reviews total — broaden terms, extend to 15 years. Low quality — prioritize high-citation reviews. Stalled — proceed with what came back and flag incomplete coverage.

---

## Stage 3: Layered Outline Synthesis

Owner: lit-pm for the high-level structure, literature-researcher agents for section detail. Checkpoint: medium complexity and above.

Input: the review collection with convergence data. You read the reviews and create `depth_profile.sections.target_count` sections, each with a title, a core thesis (a testable claim or question), and the key question it addresses. Each section's detail work goes to a literature-researcher, which proposes 2-4 subsections with specific questions, identifies which reviews and papers to cite per subsection, and notes gaps needing targeted research. Present the complete outline to the user where the checkpoint plan calls for it.

```yaml
outline:
  introduction:
    thesis: string
    assigned_to: lit-synthesizer
  sections:
    - title: string
      thesis: string
      subsections: list
      assigned_to: string  # literature-researcher agent
      key_reviews: list
      research_mandate: string
  balance_check:
    largest_section_pct: float
    smallest_section_pct: float
```

Gate: target section count, each with a specific thesis, all major review themes covered, roughly 15-30% of total per section, user approval where applicable.

Recovery: revise from user feedback; after a second rejection, return to Stage 2 with adjusted scope. Merge small sections and split large ones when unbalanced.

---

## Stage 4: Introduction Writing

Owner: lit-synthesizer, then editor for a quick polish. Automatic.

Input: the approved outline. The introduction frames the research question and its importance, previews the structure and each section's thesis, sets depth and scope expectations, and supplies context that section writers reference. The editor's pass is a clarity check and an outline-consistency check, not fact-checking.

Deliverable: a polished 2-4 paragraph introduction, available to Stage 5 section writers as context.

Gate: research question framed clearly, structure preview matches the approved outline, 2-4 paragraphs, editor polish applied.

---

## Stage 5: Parallel Section Research & Writing

Owner: literature-researcher agents, continuing from their Stage 3 sections. Gated by per-section fact-checking rather than a checkpoint.

Input: the introduction as context plus the section assignment.

Section writers are active researchers, not review summarizers. The Stage 2 reviews give them the high-level landscape (what is known, what is controversial), key research groups and seminal papers, entry points for deeper investigation, and identified gaps. From there, writers conduct targeted primary literature search against their section thesis, forward and backward citation tracking from the reviews, parameter-specific searches where quantitative data is needed, and a sweep for recent papers the reviews may have missed.

Papers per section, totalling `depth_profile.research.papers_per_section`: 10-15 foundational (established findings), 5-10 recent (last 2-3 years), 3-5 very recent (last 6-12 months) for the recency survey.

### Recency survey subsection

Each section includes a brief (1-2 paragraph) subsection on the most recent literature:

```markdown
### X.X.X Recent Developments (2025-2026)

Three recent studies have begun addressing [specific challenge]. [Author] et al.
(2025) demonstrated that [finding] [citation]. Concurrently, [another group]
reported [finding] [citation]. Most recently, [third study] identified [finding],
suggesting [implication] [citation].

These developments suggest [emerging trend], though [remaining uncertainty].
```

`depth_profile.research.recency_survey` controls the format (BRIEF/STANDARD/COMPREHENSIVE) and `depth_profile.research.quantitative_table` controls whether a quantitative table is included (OPTIONAL/RECOMMENDED/REQUIRED). Both come from the depth profile in the handoff.

### Writer autonomy

The assigned thesis is fixed and the section boundaries hold, but writers may add subsections when research reveals the need and reorder subsections where logical flow improves.

Gate: papers cited within the depth profile range and primary rather than review-only, recency survey present with 3-5 papers from the last 6-12 months, citations carry publication dates, section addresses its thesis, no placeholder text ("TODO", "[CITE]", "[INSERT]"), length per `depth_profile.sections.depth_per_section` (FOCUSED 1000-2000, STANDARD 1500-2500, COMPREHENSIVE 2000-3500; default 2000-3000 words).

Recovery: a stalled section is skipped with the gap flagged for synthesis. Low paper count goes back to the writer with specific research directions.

---

## Stage 6a: Per-Section Quick Validation (BLOCKING)

Owner: fact-checker. A section reaches synthesis only after it passes.

Checks against the completed draft:

- Primary papers cited, per `depth_profile.research.papers_per_section` lower bound
- Recency survey present with 3-5 papers from the last 6-12 months
- Citations include publication dates
- Quantitative data carries units and measurement context
- Section addresses its assigned thesis
- No contradictions with the introduction
- No placeholder text
- Length within the depth profile range

Result: PASS or REVISION-NEEDED with specific issues.

Escape hatch — cap at 3 revision cycles. On the third failure, present the user with the specific issues and these options: accept the section as-is (waiving the requirement), adjust requirements for this section, assign a different researcher, or remove the section from the outline. Record the decision in `workflow_state.quality_overrides[]`.

Minor issues go back to the writer as a revision list; major issues warrant substantial revision.

---

## Stage 6b: Comprehensive Final Fact-Check (NON-BLOCKING)

Owner: fact-checker.

Input: all sections plus the introduction. Deep checks: cross-section consistency (no contradictory claims), citation accuracy by spot-checking 10 random citations against the retrieved papers, quantitative claim verification (values match sources, units correct, context preserved), gap analysis against the stated scope, and whether measurement methods and conditions are noted for key data.

Deliverable: a revision list with priorities — P0 critical (factual errors, contradictions; fix before delivery), P1 important (missing citations, unclear claims; should fix), P2 nice-to-have (formatting, style; editor handles). The list goes to Stage 8. Major P0 issues may warrant returning to Stage 7.

---

## Stage 7: Active Synthesis & Augmentation

Owner: lit-synthesizer as senior author. Checkpoint: high-stakes only.

Input: the introduction plus all validated section drafts.

This is active curation, not assembly. Read the sections as a scientist: do they collectively answer the research question, what narrative arc do they tell, and where do they contradict, overlap, or leave gaps? Identify cross-cutting themes — patterns recurring across sections, implicit connections writers hinted at, insights visible only across the whole set.

The senior author has authority to reorder sections for logical flow, merge similar subsections across sections, add subsections to fill gaps with targeted research, and rewrite transitions and section boundaries for narrative continuity. The conclusion synthesizes key takeaways across sections, states remaining uncertainties and contradictions, draws implications back to the original research question, and identifies future directions.

Additions over 20% of content are flagged to the user with rationale, and trigger Stage 7.5.

Gate: logical narrative flow, cross-cutting themes identified and woven through, gaps filled, conclusion synthesizes across all sections, major additions flagged.

Recovery: revise from user feedback at the checkpoint; cross-section inconsistency goes back to the sections with specific feedback.

---

## Stage 8: Editorial Polish

Owner: editor. Automatic.

Input: the synthesized document plus the Stage 6b revision list. Incorporate P0 and P1 items, remove redundancy, clarify ambiguous phrasing, keep terminology consistent, normalize citation format (Nature-style inline), section numbering, and figure/table references, and smooth stylistic differences between sections.

Gate: P0/P1 revisions incorporated, consistent voice, clear and properly formatted.

If the polish pass surfaces major content problems, return to Stage 7. If it cannot finish, deliver as-is with a note about incomplete polish.

---

## Stage Summary Table

| Stage | Owner | Checkpoint |
|-------|-------|------------|
| 0: Archival | lit-pm | Never |
| 1: Scope | requirements-analyst | Always |
| 2: Reviews | lit-pm + researchers | High-Stakes |
| 3: Outline | lit-pm + researchers | Medium+ |
| 4: Intro | lit-synthesizer + editor | Never |
| 5: Sections | researchers (parallel) | Never |
| 6a: Quick FC | fact-checker | Blocking |
| 6b: Deep FC | fact-checker | Never |
| 6c: DA sections | devils-advocate | Active |
| 7: Synthesis | lit-synthesizer | High-Stakes |
| 7.5: DA synthesis | devils-advocate | Conditional |
| 8: Polish | editor | Never |

## Session Management Summary

| Event | Action |
|-------|--------|
| Workflow start | Create `/tmp/lit-pm-session-{timestamp}-{pid}/` |
| Stage 0 complete | Write `archival-guidelines-summary.md` |
| All stages complete | Delete session directory |
| Workflow abort/failure | Preserve session directory, log path |
| Resume workflow | Reuse session if exists, recreate if missing |
