# Handoff Schema

This document defines the YAML schema for handing off from perspective-swarm to lit-pm for comprehensive literature review.

## When Handoff Occurs

Handoff to lit-pm is triggered when:
1. User selects "Deep dive" option in Stage 4
2. Synthesis is complete but user wants deeper research
3. Convergence is low and user wants more evidence

## Handoff Payload Schema

```yaml
# handoff-payload.yaml
handoff:
  version: "1.0"
  timestamp: ISO8601

  source:
    skill: perspective-swarm
    workflow_id: string
    session_path: string          # /tmp/swarm-session-{id}/

  target:
    skill: lit-pm

  payload:
    original_prompt: string       # User's original question
    reframed_challenge: string    # Neutrally framed challenge
    problem_type: string          # decision | creative | analytical | strategic

    synthesis_summary: string     # Executive summary from Stage 3

    convergent_insights:          # Themes with agreement
      - theme: string
        confidence_score: float
        contributing_archetypes: [string]
        key_evidence: [string]

    divergent_insights:           # Unique perspectives
      - archetype: string
        insight: string
        confidence: integer

    key_uncertainties: [string]   # Areas needing more research

    suggested_search_terms:       # For lit-pm literature search
      - term: string
        rationale: string

    blind_spots_identified: [string]  # Aggregated from all perspectives

  context:
    perspectives_completed: integer
    convergence_level: string     # high | medium | low | none
    user_feedback: string         # Any refinement feedback from user
```

## Example Handoff

```yaml
handoff:
  version: "1.0"
  timestamp: 2026-02-04T19:00:00Z

  source:
    skill: perspective-swarm
    workflow_id: swarm-20260204-183000-a1b2c3d4
    session_path: /tmp/swarm-session-20260204-183000-a1b2c3d4/

  target:
    skill: lit-pm

  payload:
    original_prompt: "Should we expand into the European market?"
    reframed_challenge: "Evaluate market expansion opportunity: European market entry for B2B SaaS product"
    problem_type: strategic

    synthesis_summary: |
      European market expansion presents a significant growth opportunity with
      regulatory complexity as the primary risk factor. Market data supports
      demand, but implementation feasibility requires careful planning around
      GDPR compliance and local market adaptation.

    convergent_insights:
      - theme: "Growth opportunity in EU market"
        confidence_score: 7.2
        contributing_archetypes: [optimist, pragmatist]
        key_evidence:
          - "EU B2B SaaS market growing at 12% CAGR"
          - "Low competition in specific vertical"

    divergent_insights:
      - archetype: critic
        insight: "GDPR compliance costs underestimated in most expansion plans"
        confidence: 5
      - archetype: innovator
        insight: "Partnership model may be better than direct entry"
        confidence: 3

    key_uncertainties:
      - "True cost of GDPR compliance for our data model"
      - "Local competitor response timeline"
      - "Currency fluctuation impact on pricing"

    suggested_search_terms:
      - term: "B2B SaaS GDPR compliance cost studies"
        rationale: "Address critic's concern about underestimated compliance"
      - term: "European market entry partnership vs direct"
        rationale: "Explore innovator's partnership suggestion"
      - term: "SaaS expansion case studies Europe 2024-2025"
        rationale: "Find comparable examples for analyst validation"

    blind_spots_identified:
      - "May underweight execution complexity (Optimist)"
      - "May miss partnership opportunities (Pragmatist)"
      - "Limited quantitative data on our specific vertical (Analyst)"

  context:
    perspectives_completed: 5
    convergence_level: medium
    user_feedback: "Particularly interested in GDPR compliance depth"
```

## Handoff Protocol

### Pre-Handoff Validation

Before generating handoff:
1. Verify lit-pm skill is available
2. Ensure all required payload fields are populated
3. Copy critical artifacts (synthesis, perspectives) to handoff payload directly (don't rely on session path alone)

### Handoff Execution

1. Generate `handoff-payload.yaml` in session directory
2. Update workflow-state.yaml with handoff status
3. Invoke lit-pm with:
   ```
   /lit-pm --handoff {session_path}/handoff-payload.yaml
   ```
4. Wait for acknowledgment
5. On success: Mark workflow as COMPLETED
6. On failure: Present user options (retry, save standalone, manual research directions)

### Post-Handoff

The perspective-swarm session remains available for reference. lit-pm may:
- Access original perspective files
- Reference the synthesis document
- Build upon the suggested search terms

Session cleanup:
- Session path remains valid for 24 hours
- lit-pm should copy needed data to its own session
- perspective-swarm archives session after expiry
