# Handoff Schema v2.0

This document defines the generic YAML schema for handing off from perspective-swarm to any handoff-eligible workflow.

## Schema Version

**Current Version**: 2.0.0
**Breaking Changes from v1.0**: `target.skill` is now dynamic (not hardcoded to lit-pm)

## When Handoff Occurs

Handoff is triggered when:
1. User selects a handoff option (C/D/E...) in Stage 4
2. Synthesis is complete but user wants deeper exploration
3. Convergence is low and user wants more evidence/perspectives

## Skill Eligibility Metadata

For a skill to appear as a handoff option, its SKILL.md must include `handoff:` frontmatter:

```yaml
---
name: skill-name
description: ...

handoff:
  accepts_handoff: true                    # REQUIRED: Must be true
  handoff_categories: [research, analysis] # REQUIRED: At least one category
  handoff_description: "Brief description (shown in menu)"  # REQUIRED
  handoff_trigger: "{payload_path}"        # OPTIONAL: Default is just payload path
  protocol_version: "2.0"                  # OPTIONAL: Default "2.0"

  # OPTIONAL: Health check command
  health_check: "skill-name --health"

  # OPTIONAL: Capability declarations
  requires:                                # Fields this skill needs in payload
    - context.original_prompt
    - context.problem_type
  optional_consumes:                       # Nice-to-have fields
    - insights.uncertainties
    - research_seeds.open_questions
---
```

### Core Categories (Recommended)

| Category | Description | Example Skills |
|----------|-------------|----------------|
| `research` | Deep investigation, literature review | lit-pm, researcher |
| `implementation` | Code/build artifacts | programming-pm, software-developer |
| `analysis` | Quantitative/rigorous analysis | statistician, mathematician |
| `architecture` | System/technical design | systems-architect |
| `verification` | Fact-checking, validation | fact-checker |
| `creative` | Ideation, brainstorming, cross-domain | pov-expansion, perspective-swarm |

Additional free-form categories are allowed and will be matched using exact string comparison.

## Handoff Payload Schema (v2.0)

```yaml
handoff:
  version: "2.0"                    # Schema version
  timestamp: ISO8601                # When handoff was created
  expires_at: ISO8601               # Payload expiration (default: +1 hour)

  source:
    skill: perspective-swarm        # Source skill name
    workflow_id: string             # Unique identifier
    session_path: string            # /tmp/swarm-session-{id}/

  target:
    skill: string                   # Dynamic - populated at handoff time
    invocation: string              # How to invoke (from skill metadata)
    category: string                # Primary category of selected skill

  # Core context any workflow can use
  context:
    original_prompt: string         # User's original question
    reframed_challenge: string      # Neutrally framed challenge
    problem_type: string            # decision | creative | analytical | strategic
    synthesis_summary: string       # Executive summary from Stage 3

  # Structured insights (optional - target may ignore)
  insights:
    convergent:
      - theme: string
        confidence_score: float
        contributing_archetypes: [string]
        key_evidence: [string]
    divergent:
      - archetype: string
        insight: string
        confidence: integer
    uncertainties: [string]         # Areas needing more research
    blind_spots: [string]           # Aggregated from all perspectives

  # Research seeds (useful for research-type handoffs)
  research_seeds:
    suggested_terms:
      - term: string
        rationale: string
    open_questions: [string]

  # Metadata for handoff tracking
  meta:
    perspectives_completed: integer
    convergence_level: string       # high | medium | low | none
    user_feedback: string           # Any refinement feedback from user
    handoff_reason: string          # Why user chose this handoff
    handoff_chain: [string]         # Track handoff history (prevent loops)
    payload_hash: string            # sha256 for integrity verification
    payload_size_bytes: integer
```

## Validation Rules

### Required Fields (Handoff MUST fail if missing)

| Field | Type | Description |
|-------|------|-------------|
| `handoff.version` | string | Must be "2.0" |
| `handoff.timestamp` | ISO8601 | Must be valid ISO8601 |
| `handoff.source.skill` | string | Must be "perspective-swarm" |
| `handoff.source.session_path` | string | Must exist and be readable |
| `handoff.target.skill` | string | Must be non-empty |
| `handoff.context.original_prompt` | string | Must be non-empty |
| `handoff.context.problem_type` | enum | Must be one of: decision, creative, analytical, strategic |

### Optional Fields (Default values)

| Field | Default |
|-------|---------|
| `handoff.expires_at` | `timestamp + 1 hour` |
| `handoff.context.synthesis_summary` | Empty string |
| `handoff.insights.*` | Empty arrays |
| `handoff.research_seeds.*` | Empty arrays |
| `handoff.meta.handoff_chain` | `["perspective-swarm"]` |

### Cross-Field Validation

1. `handoff.meta.perspectives_completed` must match number of perspective files in session
2. `handoff.meta.convergence_level` must be one of: high, medium, low, none
3. If `handoff.target.skill` appears in `handoff.meta.handoff_chain`, warn about potential loop

## Handoff Lifecycle

### Phase 1: Discovery (Stage 1, cached for Stage 4)

```
1. DISCOVERY - Scan for eligible skills (see workflow-discovery.md)
   - Cache results in available-workflows.yaml
```

### Phase 2: Selection (Stage 4)

```
2. SELECTION - User chooses target from relevance-ranked list
```

### Phase 3: Validation (Before payload write)

```
3. VALIDATION
   - Verify target skill still exists
   - If target has health_check, run it (optional)
   - Verify session_path is accessible
```

### Phase 4: Transfer (Payload creation)

```
4. TRANSFER
   - Generate handoff-payload.yaml
   - Calculate payload_hash (sha256)
   - Write to session directory
```

### Phase 5: Invocation

```
5. INVOCATION
   - Invoke target: /skill-name {payload_path}
```

### Phase 6: Acknowledgment (Future enhancement)

```
6. ACKNOWLEDGMENT (optional)
   - Target skill confirms receipt (not required for v2.0)
```

## Example Handoff Payload

```yaml
handoff:
  version: "2.0"
  timestamp: "2026-02-04T19:30:00Z"
  expires_at: "2026-02-04T20:30:00Z"

  source:
    skill: perspective-swarm
    workflow_id: swarm-session-20260204-183000-a1b2c3d4
    session_path: /tmp/swarm-session-20260204-183000-a1b2c3d4/

  target:
    skill: lit-pm
    invocation: "/lit-pm --handoff {payload_path}"
    category: research

  context:
    original_prompt: "Should we expand into the European market?"
    reframed_challenge: "Evaluate market expansion opportunity: European market entry for B2B SaaS product"
    problem_type: strategic
    synthesis_summary: |
      European market expansion presents a significant growth opportunity with
      regulatory complexity as the primary risk factor. Market data supports
      demand, but implementation feasibility requires careful planning around
      GDPR compliance and local market adaptation.

  insights:
    convergent:
      - theme: "Growth opportunity in EU market"
        confidence_score: 7.2
        contributing_archetypes: [optimist, pragmatist]
        key_evidence:
          - "EU B2B SaaS market growing at 12% CAGR"
          - "Low competition in specific vertical"

    divergent:
      - archetype: critic
        insight: "GDPR compliance costs underestimated in most expansion plans"
        confidence: 5
      - archetype: innovator
        insight: "Partnership model may be better than direct entry"
        confidence: 3

    uncertainties:
      - "True cost of GDPR compliance for our data model"
      - "Local competitor response timeline"
      - "Currency fluctuation impact on pricing"

    blind_spots:
      - "May underweight execution complexity (Optimist)"
      - "May miss partnership opportunities (Pragmatist)"

  research_seeds:
    suggested_terms:
      - term: "B2B SaaS GDPR compliance cost studies"
        rationale: "Address critic's concern about underestimated compliance"
      - term: "European market entry partnership vs direct"
        rationale: "Explore innovator's partnership suggestion"
    open_questions:
      - "What is the average GDPR compliance cost for similar SaaS products?"
      - "What partnership models have succeeded for US B2B SaaS in EU?"

  meta:
    perspectives_completed: 5
    convergence_level: medium
    user_feedback: "Particularly interested in GDPR compliance depth"
    handoff_reason: "User selected deep literature research for regulatory clarity"
    handoff_chain: ["perspective-swarm"]
    payload_hash: "sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    payload_size_bytes: 2847
```

## Error Response Format

If handoff fails, return structured error:

```yaml
error:
  code: "INVALID_PAYLOAD" | "TARGET_NOT_FOUND" | "VALIDATION_FAILED" | "HEALTH_CHECK_FAILED"
  message: "Human-readable error description"
  details:
    missing_fields: [string]        # For INVALID_PAYLOAD
    validation_errors: [string]     # For VALIDATION_FAILED
    target_skill: string            # For TARGET_NOT_FOUND
  recoverable: boolean              # Can user retry with fixes?
  payload_preserved: string         # Path to saved payload for manual retry
```

## Session State Preservation

**Policy**: Source skill (perspective-swarm) SHOULD preserve session state until:
1. Target skill acknowledges handoff (future), OR
2. User explicitly closes session, OR
3. 24-hour expiry reached

This enables recovery if handoff fails.

## Backward Compatibility

For existing `lit-pm` integrations:
- v1.0 payloads are still valid (v2.0 is a superset)
- Target skills MUST ignore unknown fields (graceful handling)
- `target.skill: lit-pm` works exactly as before
