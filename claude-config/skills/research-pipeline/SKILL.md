---
name: research-pipeline
version: 1.0
last_updated: 2026-02-03
description: "Use when you need a QUICK automated research-to-polish chain (2-8h, fixed 5 stages, no checkpoints). NOT for adaptive multi-stage reviews with parallel agents and complexity-based checkpoints (use lit-pm for 4-24h comprehensive) or custom sequences (use technical-pm)."
prerequisites:
  - Research topic or question clearly defined
  - Access to literature databases (PubMed, bioRxiv, OpenAlex)
  - Target output format (review, analysis, summary)
  - Estimated scope (comprehensive vs focused)
success_criteria:
  - Literature reviewed and synthesized into cohesive document
  - Arguments challenged and strengthened via adversarial review
  - All citations verified and properly formatted
  - Prose polished and ready for publication/archival
  - Pipeline completion with no manual handoffs required
estimated_duration: 4-8 hours (comprehensive review), 2-4 hours (focused summary)
---

# Research Pipeline Skill

## Purpose

research-pipeline automates the research workflow by chaining five specialized skills in sequence. Instead of manually invoking researcher, synthesizer, devils-advocate, fact-checker, and editor with handoffs between each, this skill runs the whole chain with structured context passing.

This is the "pipeline pattern" for skill orchestration: a predefined sequence of skills that together accomplish a goal requiring several specialized capabilities.

## When to Use This Skill

1. **Complete research workflow needed**: from research question to polished document
2. **Quality matters**: output must withstand adversarial review and fact-checking
3. **Hands-off execution desired**: set the goal, receive the final product
4. **Standard research pipeline fits**: researcher -> synthesizer -> review -> edit matches the need

### Clear Indicators for Use

- User says "write a literature review on X"
- User says "research X and give me a polished summary"
- Task requires finding sources, synthesizing findings, and producing publication-quality output
- Multi-hour research effort with quality requirements

### When NOT to Use

- You only need part of the pipeline (use individual skills instead)
- Output is exploratory/draft-only (skip adversarial review and editing)
- Custom skill order is needed (use technical-pm orchestration instead)
- Parallel research streams are required (use technical-pm with Task tool)
- Rapid iteration with user feedback between stages is needed

## Orchestration stance

You keep the user's original goal intact across the chain and make sure each stage receives the context the previous one produced. The pipeline is a convenience, not a constraint: when a stage fails validation or genuinely needs user input, pause and escalate rather than passing weak output downstream.

Each stage runs as a separate specialist via Task tool. Fan-out conventions: `../references/delegation-and-scope.md`.

## Pipeline Architecture

```
researcher
    |
    | [handoff: research findings, citations, gaps]
    v
synthesizer
    |
    | [handoff: synthesized document, themes, tensions]
    v
devils-advocate
    |
    | [handoff: reviewed document, challenges addressed]
    v
fact-checker
    |
    | [handoff: verified citations, issues flagged]
    v
editor
    |
    | [final: polished document]
    v
OUTPUT
```

| Stage | Skill | Key Output | Quality Gate |
|-------|-------|------------|--------------|
| 1 | researcher | Literature review draft with citations | Has citations, addresses topic |
| 2 | synthesizer | Integrated analysis with themes | Cross-cutting insights present |
| 3 | devils-advocate | Strengthened arguments | Challenges addressed or documented |
| 4 | fact-checker | Verified citations | All claims have valid citations |
| 5 | editor | Polished prose | CLAUDE.md style compliant |

## Initializing the pipeline

Parse the user's goal into topic, scope (comprehensive vs focused), constraints (time, page count, focus areas), and a `workflow_id`.

### Archival compliance check

Before creating the pipeline context, read `~/.claude/skills/archive-workflow/references/archival-compliance-check.md` and apply its 5-step pattern to file creation. If that file is missing, log a warning and proceed without the archival check.

```yaml
pipeline:
  archival:
    guidelines_present: true/false
    naming_convention: "{from YAML}"
    output_directory_override: "{if archival says docs go elsewhere}"
    enforcement_mode: "advisory"
```

When setting the output location (default `docs/literature/{topic}/`), validate it against the archival structure guidelines. On a violation, present batch advisory options, record the user's choice in the pipeline context, and pass `archival_context` to every downstream stage.

### Initial context

```yaml
pipeline:
  workflow_id: "research-{uuid}"
  goal: "{user's original goal}"
  topic: "{extracted topic}"
  scope: comprehensive | focused
  constraints:
    max_pages: {N}
    focus_areas: [list]
    time_limit: {hours}
  current_stage: 1
  started_at: "{timestamp}"
```

## The five stages

Each stage is dispatched via Task tool with the previous stage's handoff file path in the prompt. Judge each output against the one-line quality gate in the table above before writing the handoff; the specialists own how they get there.

**Stage 1 — researcher.** Literature search (PubMed, bioRxiv, OpenAlex), paper notes, draft review with Nature-style citations. If the draft has no citations or misses the topic, pause and offer the user a retry with narrowed scope, acceptance of partial output, or abort.

**Stage 2 — synthesizer.** Reads the researcher output, identifies cross-cutting themes, surfaces tensions and contradictions, draws project-specific implications. Should add value beyond a summary and carry citations forward.

**Stage 3 — devils-advocate.** Identifies the thesis, evaluates strategic coherence, challenges thesis-critical claims, proposes counter-arguments, exchanges with the synthesizer for up to 2 rounds. The handoff records which challenges were addressed, which uncertainties remain, and the approval status.

**Stage 4 — fact-checker.** Inventories quantitative claims and verifies each citation against a real source, checks superscript format, flags missing or wrong citations, and produces a verification report. Minor issues get corrected inline; major ones go back to the synthesizer and are then re-checked.

**Stage 5 — editor.** Applies CLAUDE.md style guidelines: bullets to prose where appropriate, bridging transitions, glossary placement, final polish for publication.

### Completion report

```markdown
# Research Pipeline Complete

**Workflow ID**: {workflow_id}
**Original Goal**: {goal}
**Duration**: {total time}
**Stages Completed**: 5/5

## Final Output
**Location**: {path to final document}

## Pipeline Summary
| Stage | Status | Notes |
|-------|--------|-------|
| Researcher | Complete | 8 papers reviewed |
| Synthesizer | Complete | 3 themes identified |
| Devil's Advocate | Complete | 2 challenges, 1 uncertainty |
| Fact-Checker | Complete | 12 citations verified |
| Editor | Complete | CLAUDE.md compliant |

## Quality Indicators
- Citations verified: 12/12
- Challenges addressed: 2/2
- Uncertainties documented: 1
- Style compliance: PASS

## Output Files
- Final document: {path}
- Researcher notes: {path}
- Synthesis draft: {path}
- Review report: {path}
- Fact-check report: {path}
```

### Optional: git strategy advisory

After the completion report you can invoke `git-strategy-advisor` via Task tool in post-work mode:

```
Use git-strategy-advisor to determine git strategy for completed work.

mode: post-work
```

Read the advisor's `summary` field into the completion report. Skip the section silently if confidence is "none" or "low", if the advisor is unavailable or errors, or if output files were written outside the current git repository (the advisor only sees changes inside it).

## Handoff Format

Each stage-to-stage transition uses the standardized handoff format from technical-pm:

```yaml
handoff:
  version: "1.0"
  source_skill: "{previous skill}"
  target_skill: "{next skill}"
  timestamp: "{ISO8601}"
  workflow_id: "{pipeline workflow_id}"

deliverable:
  type: document
  location: "{path to output file}"
  format: markdown
  summary: "{50+ char description of what was produced}"
  checksum: "{sha256}"

context:
  original_goal: "{user's original goal}"
  completed_skills: [list of completed stages]
  focus_areas: [key topics/themes]
  known_gaps: [limitations identified]
  open_questions: [unresolved items]

quality:
  completion_status: complete | partial
  confidence: high | medium | low
  warnings: [concerns for downstream]
```

Validate before each handoff: required fields present, summary >= 50 chars, file exists, checksum recomputes and matches. On validation failure, stop the pipeline and report to the user.

## Configuration Options

### Scope Settings

**Comprehensive** (default):
- Researcher: Full literature review (8+ papers)
- Synthesizer: Multi-theme analysis
- Devil's Advocate: 2 exchange rounds
- Fact-Checker: All citations verified
- Editor: Full CLAUDE.md polish

**Focused**:
- Researcher: Targeted review (3-5 papers)
- Synthesizer: Single-theme summary
- Devil's Advocate: 1 exchange round
- Fact-Checker: Critical citations only
- Editor: Essential polish only

### Skip Options

```
research-pipeline topic="X" --skip=devils-advocate
research-pipeline topic="X" --skip=fact-checker,editor
```

Skipping stages reduces quality guarantees; note skipped stages in the completion report.

### Output Location

Default: `docs/literature/{topic}/`

Override: `research-pipeline topic="X" --output="{custom path}"`

## Error Handling

If a stage fails: keep all previous stage outputs, save the current stage's work-in-progress, stop before the next stage, and report with options.

```
Pipeline paused: Stage 3 (devils-advocate) failed

Error: Could not identify thesis in synthesizer output

Completed outputs preserved:
- docs/literature/topic/researcher-draft.md
- docs/literature/topic/synthesis.md (partial)

Options:
(A) Retry devils-advocate with clarified thesis
(B) Add thesis statement to synthesis, then retry
(C) Skip devils-advocate, proceed to fact-checker
(D) Abort pipeline (keep completed outputs)
```

### Interruption Recovery

Pipeline state is saved after each stage completes, at `/tmp/pipeline-state-{workflow_id}.yaml`.

```
Found interrupted pipeline: research-abc123
Topic: "hepatocyte oxygenation"
Progress: 3/5 stages complete

Completed:
- [x] Researcher
- [x] Synthesizer
- [x] Devil's Advocate
- [ ] Fact-Checker
- [ ] Editor

Options:
(A) Resume from Fact-Checker
(B) Restart pipeline from beginning
(C) Abort (keep completed outputs)
```

## Outputs

### Primary Output
- Final polished document: `docs/literature/{topic}/review-{topic}.md`

### Intermediate Outputs (preserved for reference)
- Researcher draft: `docs/literature/{topic}/researcher-draft.md`
- Paper notes: `docs/literature/{topic}/notes/`
- Synthesis document: `docs/literature/{topic}/synthesis-{topic}.md`
- Devil's advocate review: `docs/literature/{topic}/review-report.md`
- Fact-check report: `docs/literature/{topic}/fact-check-report.md`

### Pipeline Artifacts
- Handoff documents: `/tmp/handoff-{workflow_id}-*.yaml`
- State file: `/tmp/pipeline-state-{workflow_id}.yaml`
- Completion report: Displayed to user, optionally saved

## Integration with Other Skills

Use technical-pm instead of research-pipeline when multiple independent research streams need parallel execution, when custom skill ordering is required, when non-research skills belong in the workflow, or when dependency management gets complex.

After the pipeline completes, the user may invoke archive-workflow separately for project organization:

```
Skill(archive-workflow, project="{project root}")
```

The pipeline does not invoke archive-workflow automatically — organization decisions stay with the user. It may invoke git-strategy-advisor (via Task tool, not Skill tool) for advisory branch/push/PR recommendations.

## Example Invocations

### Basic Research Review

**User**: "Research hepatocyte oxygenation and write a comprehensive literature review"

1. Researcher: reviews 8-10 papers on hepatocyte oxygen consumption
2. Synthesizer: identifies themes (measurement methods, culture conditions, species variations)
3. Devil's Advocate: challenges the assumption that in vitro values apply to bioreactor design
4. Fact-Checker: verifies all K_oA values trace to primary sources
5. Editor: polishes into a CLAUDE.md-compliant document

**Output**: `docs/literature/hepatocyte-oxygenation/review-hepatocyte-oxygenation.md`

### Focused Summary

**User**: "Give me a quick summary of hollow fiber bioreactor designs for liver support"

Focused mode: 3-4 key papers, single-theme synthesis, 1 devil's advocate round, critical citations only, essential polish.

**Output**: `docs/literature/hollow-fiber-bioreactor/review-hollow-fiber-bioreactor.md`

### With Constraints

**User**: "Research Matrigel alternatives for hepatocyte culture, focus on chemical approaches, max 5 pages"

Topic is Matrigel alternatives for hepatocyte culture, focus is chemical rather than biological approaches, and the 5-page limit propagates: the researcher filters for chemical/synthetic matrix papers, the synthesizer respects the page limit, and every stage honors the focus constraint.

## Common Pitfalls

1. **Scope creep in researcher stage**
   - **Symptom**: Researcher spends 6+ hours on a "quick summary" request
   - **Fix**: Set explicit scope at pipeline init (comprehensive vs focused)

2. **Thesis drift between stages**
   - **Symptom**: Final document doesn't answer the original question
   - **Fix**: Handoff includes original_goal; each stage checks alignment

3. **Citation format inconsistency**
   - **Symptom**: Mixed citation formats (superscripts and brackets)
   - **Fix**: Fact-checker enforces format; editor normalizes

4. **Skipping stages reduces quality**
   - **Symptom**: Unchallenged arguments, unverified citations in final output
   - **Fix**: Warn the user when stages are skipped; note it in the completion report

## Handoffs

| Condition | Hand off to |
|-----------|-------------|
| Pipeline complete | **User** (with completion report) |
| Pipeline complete, organization needed | **archive-workflow** (manual invocation) |
| Stage failure, needs diagnosis | **User** (with error context) |
| Custom workflow needed | **technical-pm** (for flexible orchestration) |
| Parallel research streams | **technical-pm** with Task tool |

---

## Supporting Resources

See `examples/` directory for:
- `pipeline-invocation-example.md` - Sample pipeline run with all stages
- `completion-report-example.md` - Example completion report format
