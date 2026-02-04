# Inter-Agent Data Contracts

This document defines the data contracts between POV-expansion pipeline stages.

## Overview

The POV-expansion workflow uses filesystem-based state management with structured YAML outputs. Each stage produces a handoff file that the next stage consumes.

## Stage 2 → Stage 3: Abstract Representation

**Producer**: pov-abstractor-classifier
**Consumer**: pov-abstractor-classifier (Stage 3), pov-perspective-analyst (Stage 4)
**File**: `/tmp/pov-session-{id}/stage-2-abstraction.yaml`

```yaml
abstract_representation:
  version: "1.0"
  original_problem: string
  rasmussen_hierarchy:
    purpose: string               # Why the system exists
    abstract_function: string     # What flows and constraints
    generalized_function: string  # General process type
    physical_function: string     # Specific mechanism
    physical_form: string         # Concrete implementation
  structural_elements:
    - element: string
      type: enum[entity, relation, attribute, function]
      abstraction_level: enum[surface, structural, deep]
  constraints:
    - constraint: string
      type: enum[hard, soft, preference]
  success_metrics:
    - metric: string
      measurable: boolean
```

## Stage 3 → Stage 4: Domain Classification

**Producer**: pov-abstractor-classifier
**Consumer**: pov-expansion-pm (for agent dispatch), pov-perspective-analyst (Stage 4)
**File**: `/tmp/pov-session-{id}/stage-3-domains.yaml`

```yaml
domain_classification:
  version: "1.0"
  near_field:
    - domain: string
      confidence: enum[HIGH, MEDIUM, LOW]
      rationale: string
      structural_similarity: float  # 0.0-1.0
    - domain: string
      confidence: enum[HIGH, MEDIUM, LOW]
      rationale: string
      structural_similarity: float
  mid_field:
    - domain: string
      confidence: enum[HIGH, MEDIUM, LOW]
      rationale: string
      structural_similarity: float
  far_field:
    - domain: string
      confidence: enum[HIGH, MEDIUM, LOW]
      rationale: string
      structural_similarity: float
      warning: "Far transfer is empirically rare. This analogy requires rigorous evaluation."
  user_confirmation:
    timestamp: ISO8601
    action: enum[approved, replaced, skipped]
    notes: string
```

## Stage 4 → Stage 6: Perspective Reports

**Producer**: pov-perspective-analyst (4 parallel instances)
**Consumer**: pov-synthesizer (Stage 6), pov-transfer-evaluator (Stage 7)
**Files**:
- `/tmp/pov-session-{id}/stage-4-perspectives/near-field-1.md`
- `/tmp/pov-session-{id}/stage-4-perspectives/near-field-2.md`
- `/tmp/pov-session-{id}/stage-4-perspectives/mid-field.md`
- `/tmp/pov-session-{id}/stage-4-perspectives/far-field.md`

```yaml
perspective_report:
  version: "1.0"
  domain: string
  domain_distance: enum[near, mid, far]

  analogous_problems:
    - problem: string
      source_domain: string
      structural_similarity: float  # 0.0-1.0
      mapping_depth: enum[surface, structural, deep]
      preserved_relations:
        - relation: string
          confidence: float

  solutions:
    - solution: string
      mechanism: string  # How it works
      evidence: string  # Citations or examples
      transfer_barriers: [string]
      adaptation_cost: enum[low, medium, high]
      prerequisites: [string]

  structural_analogies:
    - source_structure: string
      target_structure: string
      relationship_preserved: string
      confidence: float

  convergence_markers:
    - marker: string  # Theme, strategy, or pattern
      description: string
      relevance: float
```

## Stage 6 → Stage 7: Convergence Results

**Producer**: pov-synthesizer (Stage 6)
**Consumer**: pov-synthesizer (Stage 8), pov-transfer-evaluator (Stage 7)
**File**: `/tmp/pov-session-{id}/stage-6-convergence.md`

```yaml
convergence_results:
  version: "1.0"
  high_confidence:  # 3-4 reports
    - concept: string
      reports: [string]  # e.g., ["P1", "P3", "P4"]
      confidence: float  # 0.9-1.0
      convergence_type: enum[theme, strategy, structure, mechanism]
      relevance: float

  medium_confidence:  # 2 reports
    - concept: string
      reports: [string]
      confidence: float  # 0.6-0.8
      convergence_type: enum[theme, strategy, structure, mechanism]
      relevance: float

  novel_insights:  # 1 report but high relevance
    - concept: string
      reports: [string]
      confidence: float  # 0.3-0.5
      note: string

  convergence_floor_met: boolean
  convergence_floor_reason: string
```

## Stage 7 → Stage 8: Transfer Evaluation

**Producer**: pov-transfer-evaluator
**Consumer**: pov-synthesizer (Stage 8)
**File**: `/tmp/pov-session-{id}/stage-7-transfer.md`

```yaml
transfer_evaluation:
  version: "1.0"
  evaluations:
    - analogy_id: string
      analogy: string
      source_domain: string
      domain_distance: enum[near, mid, far]

      scores:
        structural_depth: float  # 0.0-1.0
        mechanism_clarity: float
        barrier_awareness: float
        adaptation_requirements: float
        evidence_quality: float
        overall: float  # weighted average

      classification: enum[quick_win, strategic_bet, easy_win, avoid]

      impact_assessment:
        level: enum[high, low]
        rationale: string

      barriers:
        - barrier: string
          severity: enum[low, medium, high]
          mitigation: string

      adaptation:
        cost: enum[low, medium, high]
        requirements: [string]
        estimated_effort: string

      far_field_warning: string | null  # REQUIRED if domain_distance = far

      recommendation: string

  summary:
    total_analogies_evaluated: integer
    quick_wins: integer
    strategic_bets: integer
    easy_wins: integer
    avoid: integer
    average_feasibility: float
```

## Handoff Protocol

Each stage produces a handoff file to enable resumability and debugging:

```yaml
handoff:
  stage: integer              # Stage number (1-11)
  producer: string            # Agent name
  consumer: string            # Next agent name
  status: enum[complete, partial, failed]
  output_path: string         # Path to output file
  timestamp: ISO8601
  duration_seconds: integer
  key_findings:               # Brief summary for progress tracking
    - string
  open_questions:             # Issues for next stage
    - string
  quality_score: integer      # 1-5 self-assessment
```

**Location**: `/tmp/pov-session-{id}/handoffs/stage-{N}-handoff.yaml`

## Workflow State Schema

**Producer**: pov-expansion-pm (maintained throughout)
**File**: `/tmp/pov-session-{id}/workflow-state.yaml`

```yaml
workflow_state:
  workflow_id: string         # e.g., "pov-20260203-154530"
  problem_statement: string
  started_at: ISO8601
  current_stage: integer      # 1-11
  stages_completed: [integer]

  artifacts:
    stage_1: string | null    # Path to output file
    stage_2: string | null
    stage_3: string | null
    stage_4: string | null
    stage_5: string | null
    stage_6: string | null
    stage_7: string | null
    stage_8: string | null
    stage_9: string | null
    stage_10: string | null
    stage_11: string | null

  status: enum[in_progress, paused, completed, failed]
  last_checkpoint: ISO8601
  error_log:
    - timestamp: ISO8601
      stage: integer
      error: string
      resolution: string

  user_decisions:             # Track user input for traceability
    domain_classification:
      timestamp: ISO8601
      action: enum[approved, replaced, skipped]
      notes: string
```

## Progress Tracking Schema

**Producer**: pov-expansion-pm
**File**: `/tmp/pov-session-{id}/progress.md`

Format (Markdown):
```markdown
# POV-Expansion Progress

**Workflow ID**: pov-session-{id}
**Problem**: [Brief problem statement]
**Started**: {timestamp}
**Current Stage**: {N} of 11

## Stage Progress

| Stage | Status | Started | Completed | Duration |
|-------|--------|---------|-----------|----------|
| 1 | COMPLETE | {time} | {time} | {duration} |
| 2 | IN_PROGRESS | {time} | - | {elapsed} |
...

## Recent Events

- {timestamp}: Stage 2 started
- {timestamp}: Stage 1 completed (Quality Gate 1: PASS)
...
```

## Verification Output Schema

**Producer**: fact-checker (Stage 5, 2 parallel instances)
**Consumer**: pov-synthesizer (Stage 6)
**Files**:
- `/tmp/pov-session-{id}/stage-5-verification/verification-1.md`
- `/tmp/pov-session-{id}/stage-5-verification/verification-2.md`

Format (Markdown with embedded YAML):
```yaml
verification_results:
  version: "1.0"
  perspectives_checked: [string]  # Which reports were verified
  claims_verified: integer
  errors_found: integer

  findings:
    - claim: string
      location: string  # e.g., "near-field-1.md line 45"
      status: enum[verified, incorrect, unverifiable]
      correction: string | null
      confidence: float

  overall_assessment: enum[high_confidence, medium_confidence, low_confidence]
```

## Data Flow Diagram

```
Stage 1 (requirements-analyst)
  ↓ [stage-1-refinement.md]
Stage 2 (pov-abstractor-classifier)
  ↓ [stage-2-abstraction.yaml]
Stage 3 (pov-abstractor-classifier)
  ↓ [stage-3-domains.yaml]
Stage 4 (pov-perspective-analyst x4) [PARALLEL]
  ↓ [stage-4-perspectives/*.md]
Stage 5 (fact-checker x2) [PARALLEL]
  ↓ [stage-5-verification/*.md]
Stage 6 (pov-synthesizer)
  ↓ [stage-6-convergence.md]
Stage 7 (pov-transfer-evaluator)
  ↓ [stage-7-transfer.md]
Stage 8 (pov-synthesizer)
  ↓ [stage-8-synthesis.md]
Stage 9 (strategist)
  ↓ [stage-9-strategic.md]
Stage 10 (devils-advocate)
  ↓ [stage-10-review.md]
Stage 11 (editor)
  ↓ [stage-11-final.md] [DELIVERABLE]
```

## Validation Rules

### Required Fields

Each data contract specifies required vs optional fields:
- `version`: REQUIRED (for schema evolution)
- Core data fields: REQUIRED (varies by schema)
- Metadata fields: OPTIONAL

### Type Safety

Agents MUST validate:
- Enums match allowed values
- Floats are in specified ranges (e.g., 0.0-1.0)
- Arrays are non-empty where specified
- Timestamps are valid ISO8601

### Backward Compatibility

Schema version changes:
- **Patch** (1.0 → 1.1): Add optional fields only
- **Minor** (1.0 → 2.0): Add required fields with defaults
- **Major** (1.0 → 2.0): Breaking changes (avoid)

## Error Handling

### Malformed Data

If consumer receives malformed data:
1. Log error to `workflow-state.yaml` error_log
2. Attempt graceful degradation (use partial data if possible)
3. If critical: Escalate to user via pov-expansion-pm

### Missing Files

If expected file missing:
1. Check handoff status (was stage completed?)
2. If stage marked complete but file missing: CRITICAL ERROR
3. If stage in progress: Wait with timeout
4. After timeout: Escalate

### Schema Version Mismatch

If consumer expects v1.0 but receives v2.0:
1. Check if backward compatible (added fields only)
2. If compatible: Proceed with warning
3. If incompatible: Fail with clear error message

## Usage Examples

### Reading Abstract Representation (Stage 4 Agent)

```python
import yaml

# Read abstract problem
with open(f'/tmp/pov-session-{workflow_id}/stage-2-abstraction.yaml') as f:
    abstract = yaml.safe_load(f)

purpose = abstract['abstract_representation']['rasmussen_hierarchy']['purpose']
constraints = abstract['abstract_representation']['constraints']

# Use for domain research
print(f"Looking for solutions to: {purpose}")
```

### Writing Perspective Report (Stage 4 Agent)

```python
import yaml

report = {
    'perspective_report': {
        'version': '1.0',
        'domain': 'Healthcare Patient Retention',
        'domain_distance': 'mid',
        'analogous_problems': [...],
        'solutions': [...],
        'structural_analogies': [...],
        'convergence_markers': [...]
    }
}

output_path = f'/tmp/pov-session-{workflow_id}/stage-4-perspectives/mid-field.md'
with open(output_path, 'w') as f:
    f.write('# Perspective Report: Healthcare Patient Retention\n\n')
    f.write('```yaml\n')
    yaml.dump(report, f)
    f.write('```\n')
```

### Reading Convergence Results (Stage 7 Agent)

```python
import yaml

# Read convergence analysis
with open(f'/tmp/pov-session-{workflow_id}/stage-6-convergence.md') as f:
    content = f.read()
    yaml_block = content.split('```yaml')[1].split('```')[0]
    convergence = yaml.safe_load(yaml_block)

high_conf_themes = convergence['convergence_results']['high_confidence']

# Prioritize high-confidence insights for transfer evaluation
for theme in high_conf_themes:
    print(f"HIGH CONFIDENCE: {theme['concept']} (from {len(theme['reports'])} reports)")
```

## Best Practices

1. **Always include version**: Enables schema evolution
2. **Use descriptive IDs**: Not just "1", but "healthcare-continuity"
3. **Preserve units**: "30%" not "0.3" for percentages
4. **Document assumptions**: In notes/rationale fields
5. **Validate before write**: Check schema compliance
6. **Log handoffs**: Create handoff file for every stage
7. **Enable resume**: Workflow state supports pause/resume
