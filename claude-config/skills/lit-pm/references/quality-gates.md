# Quality Gates

Per-stage quality validation criteria for the lit-pm pipeline.

## Contents

- Gate Types
- Per-Stage Quality Gates (1, 2, 3, 4, 5, 6a, 6b, 6c, 7, 7.5, 8)
- Quality Floor
- Gate Result Schema
- Escape Hatches
- Threshold Provenance

---

## Gate Types

Automated gates are checkable without judgment: count thresholds (papers, words, sections), presence checks (recency survey exists, all sections present), pattern matches (no placeholders), and range checks (word count, section balance).

Judgment gates need an agent's read: thesis specificity and claim accuracy (fact-checker), narrative flow and cross-cutting themes (lit-synthesizer).

---

## Per-Stage Quality Gates

### Stage 1: Scope Refinement

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Research question present | Presence | Not empty | Yes | Yes |
| Research question specific | Semantic | Not "what is X?" | No | Yes |
| Success criteria defined | Count | >= 1 | Yes | Yes |
| In-scope defined | Count | >= 1 | Yes | Yes |
| Out-of-scope defined | Count | >= 0 | Yes | No |
| User approval | User action | Approved | Yes | Yes |

### Stage 2: Review Discovery

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Total reviews collected | Count | >= 6 | Yes | No |
| Minimum reviews | Count | >= 4 | Yes | Yes |
| Convergence achieved | Count | >= 2 | Yes | Yes |
| Reviews with annotations | Count | = total | Yes | Yes |
| Coverage of themes | Semantic | All major themes | No | Yes |

### Stage 3: Outline Synthesis

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Section count | Range | [depth_profile.sections.target_count, default 3-5] (Simple: 2-3, Medium: 3-4, Complex: 3-5, High-Stakes: 4-6) | Yes | Yes |
| Minimum sections | Count | >= 2 | Yes | Yes |
| Each section has thesis | Presence | All sections | Yes | Yes |
| Theses are specific | Semantic | Testable claims | No | Yes |
| Section balance | Range | 15-40% each | Yes | No |
| Max section size | Range | < 50% total | Yes | Yes |
| User approval (if checkpoint) | User action | Approved | Yes | Yes |

### Stage 4: Introduction

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Introduction present | Presence | Not empty | Yes | Yes |
| Word count | Range | Per depth_profile.synthesis.introduction_scope: BRIEF 150-400, STANDARD 300-800 (default), COMPREHENSIVE 500-1200 | Yes | No |
| Structure preview present | Semantic | Sections mentioned | No | Yes |
| Matches outline | Semantic | Consistent | No | Yes |
| Editor polish applied | Presence | True | Yes | Yes |

### Stage 5: Section Writing

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Primary papers cited | Count | >= depth_profile.research.papers_per_section lower bound (Simple: >= 10, default: >= 15) | Yes | Yes |
| Maximum papers | Count | <= depth_profile.research.papers_per_section upper bound (default <= 30) | Yes | No |
| Recency survey present | Presence | Subsection exists (format per depth_profile.research.recency_survey) | Yes | Yes |
| Recency papers (6-12 mo) | Count | >= 3 | Yes | Yes |
| Section addresses thesis | Semantic | Agent judgment | No | Yes |
| No contradictions with intro | Semantic | Agent judgment | No | Yes |
| No placeholder text | Pattern | 0 matches | Yes | Yes |
| Word count | Range | Per depth_profile.sections.depth_per_section: FOCUSED 1000-2000, STANDARD 1500-2500 (default), COMPREHENSIVE 2000-3500 | Yes | No |

**Placeholder Patterns to Detect**:
```regex
TODO|FIXME|\[CITE\]|\[INSERT\]|\[TBD\]|\[PLACEHOLDER\]|XXX
```

### Stage 6a: Quick Validation (Per-Section)

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Paper count | Count | >= depth_profile.research.papers_per_section lower bound (Simple: >= 10, default: >= 15) | Yes | Yes |
| Recency survey | Presence | True | Yes | Yes |
| Recency paper count | Count | >= 3 | Yes | Yes |
| Citations have dates | Pattern | All citations | Yes | No |
| Thesis addressed | Semantic | Agent judgment | No | Yes |
| No placeholders | Pattern | 0 matches | Yes | Yes |
| Word count range | Range | Per depth_profile.sections.depth_per_section: FOCUSED 1000-2000, STANDARD 1500-2500, COMPREHENSIVE 2000-3500 (default 1500-3500) | Yes | No |

### Stage 6b: Comprehensive Fact-Check

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Cross-section consistency | Semantic | No contradictions | No | P0 |
| Citation accuracy (spot-check) | Manual | 10 random | No | P0 if major |
| Quantitative verification | Manual | Values match sources | No | P0 if wrong |
| Gap analysis | Semantic | No obvious gaps | No | P1 |
| Methodological context | Semantic | Present for key data | No | P2 |

Priority levels: **P0** critical (fix before delivery), **P1** important (should fix), **P2** nice-to-have (editor handles).

### Stage 6c: Devil's Advocate Section Review

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| DA review executed | Presence | True | Yes | Yes |
| Thesis identified | Presence | Not empty | No | Yes |
| Strategic challenges addressed | Semantic | All resolved OR documented | No | Yes |
| Exchange count | Range | 1-2 | Yes | Yes |
| Uncertainty documented (if 2 exchanges) | Presence | Required if unresolved | No | Yes |

Passes when all strategic challenges are resolved, or when 2 exchanges are complete with `uncertainty_notes` recording what stayed unresolved.

```yaml
stage_6c_gate_result:
  section_id: string
  status: enum  # PASS | PASS_WITH_UNCERTAINTY
  thesis_identified: string
  exchanges_completed: integer
  challenges:
    strategic_count: integer
    tactical_count: integer
    resolved_count: integer
  uncertainty_notes: list | null
```

### Stage 7: Synthesis

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Narrative flow logical | Semantic | Sections build | No | Yes |
| Cross-cutting themes identified | Presence | >= 1 | No | Yes |
| Conclusion present | Presence | Not empty | Yes | Yes |
| Conclusion synthesizes findings | Semantic | All sections referenced | No | Yes |
| Major additions flagged | Presence | If >20% new | Yes | Yes |
| User approval (if checkpoint) | User action | Approved | Yes | Yes |

### Stage 7.5: Devil's Advocate Synthesis Review

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| DA synthesis review executed | Presence | True (if triggered) | Yes | Yes |
| Document thesis identified | Presence | Not empty | No | Yes |
| Cross-section coherence evaluated | Semantic | All sections covered | No | Yes |
| Strategic challenges addressed | Semantic | All resolved OR documented | No | Yes |
| Exchange count | Range | 1-2 | Yes | Yes |
| Uncertainty documented (if 2 exchanges) | Presence | Required if unresolved | No | Yes |

Passes when thesis-coherence and cross-cutting-theme challenges are resolved, or when 2 exchanges are complete with strategic uncertainty notes.

```yaml
stage_7_5_gate_result:
  status: enum  # PASS | PASS_WITH_UNCERTAINTY | SKIPPED
  trigger_reason: string | null  # Why 7.5 ran (or why skipped)
  document_thesis: string
  exchanges_completed: integer
  strategic_assessment:
    thesis_coherence: enum  # STRONG | ADEQUATE | WEAK_DOCUMENTED
    cross_cutting_themes: enum  # VALID | QUESTIONABLE_DOCUMENTED
    argument_flow: enum  # LOGICAL | NEEDS_WORK_NOTED
  uncertainty_notes: list | null
```

### Stage 8: Editorial Polish

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| P0 revisions incorporated | Count | All P0s | Yes | Yes |
| P1 revisions incorporated | Count | All P1s | Yes | No |
| Voice consistent | Semantic | Agent judgment | No | Yes |
| Formatting consistent | Pattern | Uniform citations | Yes | No |

---

## Quality Floor

These hold even under `--full-auto`, and violations escalate to the user regardless of checkpoint plan or automation settings:

```yaml
quality_floor:
  stage_2:
    minimum_reviews: 4
    minimum_convergence: 1  # At least 1 review found by multiple agents
    on_failure: "HALT: Review discovery below minimum threshold. Found {count} reviews, need 4."

  stage_3:
    minimum_sections: 2
    max_section_imbalance: 50%  # No section > 50% of total
    on_failure: "HALT: Outline severely unbalanced. Largest section is {pct}% of total."

  stage_5:
    minimum_papers_per_section: 10
    recency_survey_required: true
    on_failure: "HALT: Section research insufficient. {section} has only {count} papers."

  stage_6a:
    no_placeholders: true
    on_failure: "HALT: Section contains placeholder text. Found: {placeholders}"

  stage_6c:
    minimum_da_engagement: true  # DA must actually run
    thesis_must_be_identified: true  # Cannot pass without identifying thesis
    on_failure: "HALT: Devil's advocate could not identify thesis for section {id}. Manual intervention required."
```

---

## Gate Result Schema

```yaml
gate_result:
  stage: integer
  timestamp: ISO8601
  overall_status: enum  # PASS | FAIL | PASS_WITH_WARNINGS
  checks:
    - name: string
      type: enum  # automated | judgment
      blocking: boolean
      passed: boolean
      value: any  # actual value
      threshold: any  # expected value
      details: string | null
  warnings: list
  action_taken: string | null
```

A failed blocking check goes back to the producer with the specific issues. A failed non-blocking check is recorded as a warning and the pipeline continues.

---

## Escape Hatches

**Stage 6a** caps revision cycles at 3 per section. On the third failure, present the accumulated issues and offer: accept the section as-is (waiving the requirement), adjust requirements for this section, assign a different researcher, or remove the section. Record the outcome in `workflow_state.quality_overrides[]` with the section, decision, timestamp, and waived issues.

**Devil's advocate stages** escalate when Stage 6c leaves more than half of sections with unresolved strategic challenges (or identifies a thesis as fundamentally weak), or when Stage 7.5 marks document thesis coherence WEAK or cross-cutting themes INVALID after 2 exchanges. Offer: accept with uncertainty notes and proceed to the editor, return to section writing with specific guidance, return to the outline stage to restructure, or abort if the issues are fundamental. Record the outcome in `workflow_state.da_escape_decisions[]`.

---

## Threshold Provenance

Thresholds come from literature review practice (15-30 papers per section is standard), prior technical-pm experience (word count ranges), and observed error rates (placeholder detection patterns). If they prove too strict or too loose in practice, log the violations with context and update this document with the new values and rationale.
