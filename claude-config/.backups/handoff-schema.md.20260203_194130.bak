# Handoff Schema

YAML schema for stage-to-stage communication in the lit-pm pipeline.

---

## Base Schema

All handoffs include these common fields:

```yaml
handoff:
  version: "1.0"
  stage: integer           # 1-8
  status: enum             # pending | in_progress | complete | failed
  producer: string         # skill name that produced this handoff
  consumer: string         # skill name that will receive this handoff
  workflow_id: string      # unique workflow identifier
  timestamp: ISO8601       # when handoff was created
```

---

## Stage-Specific Schemas

### Stage 1 -> Stage 2: Scope to Review Discovery

```yaml
stage_1_to_2:
  scope:
    research_question: string
    success_criteria:
      - string
    in_scope:
      - string
    out_of_scope:
      - string
  complexity:
    tier: enum  # Simple | Medium | Complex | High-Stakes
    indicators:
      scope_breadth: string
      stakes_level: string
      user_hints: string
  checkpoint_plan:
    stage_1: true  # Always
    stage_2: boolean
    stage_3: boolean
    stage_7: boolean
  user_approval:
    approved: boolean
    timestamp: ISO8601
    notes: string | null
```

### Stage 2 -> Stage 3: Reviews to Outline

```yaml
stage_2_to_3:
  reviews:
    - title: string
      doi: string | null
      pmid: string | null
      authors: string
      year: integer
      journal: string
      source_agent: string    # which agent found this review
      convergence_count: integer  # how many agents found it
      annotation: string      # key findings, coverage, quality
  convergence_analysis:
    high_priority: list       # reviews found by 3/3 agents
    medium_priority: list     # reviews found by 2/3 agents
    unique_perspectives: list # reviews found by 1 agent
    themes_covered:
      - theme: string
        reviews: list
    gaps_identified:
      - gap: string
        severity: enum  # minor | significant
```

### Stage 3 -> Stage 4: Outline to Introduction

```yaml
stage_3_to_4:
  outline:
    introduction:
      thesis: string
      assigned_to: lit-synthesizer
    sections:
      - id: string
        title: string
        thesis: string
        subsections:
          - title: string
            questions: list
        assigned_to: string
        key_reviews: list
        research_mandate: string
    balance:
      section_weights:
        section_1: float  # percentage of total
        section_2: float
      balanced: boolean
  user_approval:
    approved: boolean
    timestamp: ISO8601
    modifications: list | null
```

### Stage 4 -> Stage 5: Introduction to Section Writing

```yaml
stage_4_to_5:
  introduction:
    content: string  # markdown
    word_count: integer
    editor_polished: boolean
  section_assignments:
    - section_id: string
      assigned_to: string
      thesis: string
      research_mandate: string
      key_reviews: list
      context_from_intro: string  # relevant framing
```

### Stage 5 -> Stage 6a: Section to Quick Validation

```yaml
stage_5_to_6a:
  section:
    id: string
    title: string
    content: string  # markdown
    word_count: integer
    paper_count: integer
    recency_survey_present: boolean
    recency_paper_count: integer
    writer_notes: string | null
```

### Stage 6a -> Stage 5 (Revision) or Stage 7 (Pass)

```yaml
stage_6a_result:
  section_id: string
  status: enum  # PASS | REVISION_NEEDED
  checks:
    paper_count:
      passed: boolean
      value: integer
      threshold: 15
    recency_survey:
      passed: boolean
      paper_count: integer
    thesis_addressed:
      passed: boolean
      notes: string
    no_placeholders:
      passed: boolean
      found: list | null
    word_count:
      passed: boolean
      value: integer
  revision_needed:
    issues:
      - type: string
        severity: enum  # minor | major
        description: string
    cycle_count: integer
    max_cycles: 3
```

### Stage 5 (All Sections) -> Stage 6b: Comprehensive Fact-Check

```yaml
stage_5_to_6b:
  sections:
    - id: string
      title: string
      content: string
      validation_status: PASS
  introduction:
    content: string
```

### Stage 6b -> Stage 7: Fact-Check Results to Synthesis

```yaml
stage_6b_to_7:
  revision_list:
    p0_critical:
      - section: string
        issue: string
        citation: string | null
        recommendation: string
    p1_important:
      - section: string
        issue: string
        recommendation: string
    p2_nice_to_have:
      - section: string
        issue: string
        recommendation: string
  cross_section_issues:
    contradictions:
      - sections: list
        claim_1: string
        claim_2: string
        recommendation: string
    gaps:
      - description: string
        severity: string
```

### Stage 7 -> Stage 8: Synthesis to Editorial Polish

```yaml
stage_7_to_8:
  document:
    introduction: string
    sections:
      - id: string
        title: string
        content: string
    conclusion: string
  synthesis_notes:
    cross_cutting_themes:
      - theme: string
        sections_involved: list
        synthesis_added: string
    structural_changes:
      - type: enum  # reorder | merge | add | rewrite
        description: string
    major_additions:
      added: boolean
      percentage: float
      rationale: string | null
  revision_list_from_6b: object  # Pass through for editor
```

### Stage 8 -> Delivery: Final Document

```yaml
stage_8_final:
  document:
    title: string
    content: string  # Full markdown document
    word_count: integer
    sections: integer
    papers_cited: integer
  quality_summary:
    p0_resolved: boolean
    p1_resolved: boolean
    voice_consistent: boolean
    final_read_complete: boolean
  workflow_summary:
    workflow_id: string
    duration_total: string  # ISO 8601 duration
    complexity_tier: string
    checkpoints_used: list
    quality_overrides: list | null
```

---

## Validation Rules

### Required Fields per Transition

| Transition | Required Fields |
|------------|-----------------|
| 1 -> 2 | research_question, complexity.tier, checkpoint_plan |
| 2 -> 3 | reviews (min 4), convergence_analysis |
| 3 -> 4 | outline.sections (min 2), user_approval |
| 4 -> 5 | introduction.content, section_assignments |
| 5 -> 6a | section.content, section.paper_count |
| 6a -> 5/7 | status, checks |
| 5 -> 6b | all sections with PASS status |
| 6b -> 7 | revision_list |
| 7 -> 8 | document, synthesis_notes |
| 8 -> delivery | document.content, quality_summary |

### Validation Protocol

Before accepting a handoff:
1. Parse YAML
2. Check all required fields present
3. Validate enum values
4. Check cross-references (e.g., section_ids match)
5. Log validation result

On validation failure:
- Log specific missing/invalid fields
- Return to producer with error
- Do NOT proceed to next stage

---

## Example Handoff: Stage 2 -> Stage 3

```yaml
handoff:
  version: "1.0"
  stage: 2
  status: complete
  producer: literature-researcher-agents
  consumer: lit-pm
  workflow_id: "lit-review-cart-manufacturing-20260203"
  timestamp: "2026-02-03T14:30:00Z"

stage_2_to_3:
  reviews:
    - title: "CAR-T Manufacturing: Challenges and Solutions"
      doi: "10.1016/j.tibtech.2023.05.001"
      pmid: "37234567"
      authors: "Smith et al."
      year: 2023
      journal: "Trends in Biotechnology"
      source_agent: "literature-researcher-agent-1"
      convergence_count: 3
      annotation: "Comprehensive review of manufacturing challenges. Strong coverage of vector production and expansion. Key reference for Section 1."

    - title: "Scaling CAR-T Cell Production"
      doi: "10.1038/s41587-022-01234-5"
      authors: "Zhang et al."
      year: 2022
      journal: "Nature Biotechnology"
      source_agent: "literature-researcher-agent-2"
      convergence_count: 2
      annotation: "Focus on bioreactor scale-up. Good quantitative data on cell densities."

  convergence_analysis:
    high_priority:
      - "CAR-T Manufacturing: Challenges and Solutions"
    medium_priority:
      - "Scaling CAR-T Cell Production"
      - "Lentiviral Vector Production Review"
    unique_perspectives:
      - "Clinical-Scale T-Cell Expansion"
      - "Quality Control in CAR-T Manufacturing"
    themes_covered:
      - theme: "Vector production"
        reviews: ["review_1", "review_3"]
      - theme: "Cell expansion"
        reviews: ["review_2", "review_4"]
    gaps_identified:
      - gap: "Regulatory considerations"
        severity: minor
```
