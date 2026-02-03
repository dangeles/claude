# Handoff Document Format

## Overview

Every skill-to-skill handoff in an orchestrated workflow MUST use this standardized format. This ensures reliable context passing and prevents silent failures.

## Schema (YAML)

```yaml
# Required fields
handoff:
  version: "1.0"                    # Schema version for compatibility
  source_skill: string              # Skill that produced this handoff
  target_skill: string              # Skill that will receive this handoff
  timestamp: ISO8601                # When handoff was created
  workflow_id: string               # Unique ID for the workflow session

# Deliverable information
deliverable:
  type: document | data | analysis  # What was produced
  location: string                  # File path to output
  format: markdown | json | yaml    # Output format
  summary: string                   # Brief description (min 50 chars)
  checksum: sha256                  # For corruption detection

# Context for receiving skill
context:
  original_goal: string             # User's original request
  completed_skills: list[string]    # Skills already executed
  focus_areas: list[string]         # What to pay attention to
  known_gaps: list[string]          # Identified gaps or limitations
  open_questions: list[string]      # Unresolved items

# Quality indicators
quality:
  completion_status: complete | partial | failed
  confidence: high | medium | low
  warnings: list[string]            # Any concerns for downstream
```

## Validation Rules

Before passing handoff to next skill, validate:

1. **Schema validation**: All required fields present
2. **Content validation**:
   - `deliverable.summary` >= 50 characters
   - `deliverable.location` file exists
   - `context.completed_skills` is non-empty
3. **Checksum validation**: Recompute and compare
4. **Cross-reference**: Deliverable file actually exists

## On Validation Failure

STOP workflow immediately. Do NOT proceed with invalid handoff.

Display to user:
```
Handoff validation failed: [specific error]
Source skill: [name]
Target skill: [name]

Options:
(A) Regenerate handoff from source skill
(B) Manually fix handoff document
(C) Abort workflow (preserve completed outputs)
```

## Example Handoff

```yaml
handoff:
  version: "1.0"
  source_skill: researcher
  target_skill: synthesizer
  timestamp: "2026-02-03T14:30:00Z"
  workflow_id: "workflow-abc123"

deliverable:
  type: document
  location: "docs/literature/hepatocyte-oxygenation/review-draft.md"
  format: markdown
  summary: "Literature review covering 8 key papers on hepatocyte oxygen consumption rates in various culture formats"
  checksum: "sha256:abc123..."

context:
  original_goal: "Write comprehensive literature review on hepatocyte oxygenation"
  completed_skills: ["researcher"]
  focus_areas: ["oxygen consumption rates", "culture format comparison", "measurement methods"]
  known_gaps: ["Limited data on 3D spheroid formats", "No papers after 2024"]
  open_questions: ["Should we include non-mammalian hepatocytes?"]

quality:
  completion_status: complete
  confidence: high
  warnings: ["Some older papers (pre-2020) may have outdated methods"]
```

## Handoff Generation Protocol

When generating a handoff document:

1. **Create workflow_id** if not exists: `workflow-{uuid4()[:8]}`
2. **Compute checksum**: `sha256sum deliverable.location`
3. **Write summary**: Must be descriptive, minimum 50 characters
4. **Include all context**: Don't assume receiving skill knows background
5. **List known gaps**: Be explicit about limitations
6. **Set quality indicators**: Be honest about confidence level

## Receiving a Handoff

When receiving a handoff as the target skill:

1. **Validate immediately**: Run all validation rules before starting work
2. **Read context carefully**: Understand original goal and focus areas
3. **Note known gaps**: Plan to address or acknowledge in output
4. **Check open questions**: Decide whether to resolve or pass through
5. **Verify deliverable**: Confirm you can read the referenced file

## Versioning

Schema version "1.0" is current. Future versions will be backward compatible.

If you receive a handoff with higher version number:
- Process known fields normally
- Ignore unknown fields
- Log warning for downstream awareness
