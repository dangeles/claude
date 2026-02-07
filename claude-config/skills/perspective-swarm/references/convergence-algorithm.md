# Convergence Algorithm

This document specifies the confidence-weighted convergence algorithm used in Stage 3 synthesis.

## Overview

The algorithm:
1. Groups insights by semantic similarity
2. Applies convergence multipliers based on how many archetypes agree
3. Weights by self-reported confidence
4. Ranks insights by final weighted score
5. Preserves divergent (unique) insights with attribution

## Algorithm Specification

### Step 1: Insight Extraction

From each perspective output, extract:
- `archetype`: string (optimist | critic | analyst | innovator | pragmatist)
- `key_insight`: string (1-2 sentences)
- `confidence`: integer (1-5)
- `evidence`: list of strings
- `research_backed`: boolean (true if WebSearch succeeded)

### Step 2: Semantic Grouping

Group insights by thematic similarity using LLM-based theme extraction:

**Primary Method: LLM Theme Extraction**
```
Prompt: "Analyze these {N} insights and group them by thematic similarity.
Return a JSON structure with:
- theme_id: unique identifier
- theme_description: 1-sentence description of the theme
- insight_ids: list of insight IDs belonging to this theme

Insights:
{formatted_insights}
"
```

**Fallback Method: Keyword Overlap (Jaccard Similarity)**
```python
def keyword_overlap_jaccard(insight_a, insight_b):
    keywords_a = extract_keywords(insight_a)
    keywords_b = extract_keywords(insight_b)
    intersection = len(keywords_a & keywords_b)
    union = len(keywords_a | keywords_b)
    return intersection / union if union > 0 else 0

# Group if Jaccard similarity > 0.3
if keyword_overlap_jaccard(i, j) > 0.3:
    group_together(i, j)
```

### Step 2a: LLM Grouping Implementation

The LLM grouping call is performed inline by the brainstorming-pm orchestrator during Stage 3 convergence analysis. This is not delegated to a separate agent because it requires cross-perspective access to all insight outputs simultaneously.

See `../brainstorming-pm/references/model-selection.md` for full model selection rationale and fallback chain details.

**JSON Schema for Grouping Output**:

```json
{
  "type": "object",
  "properties": {
    "themes": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "theme_id": { "type": "string" },
          "theme_description": { "type": "string" },
          "insight_ids": { "type": "array", "items": { "type": "string" } }
        },
        "required": ["theme_id", "theme_description", "insight_ids"]
      }
    }
  },
  "required": ["themes"]
}
```

Each `insight_id` maps to `{archetype}_insight` (e.g., `optimist_insight`, `critic_insight`). After grouping, verify that all returned `insight_ids` match the expected `{archetype}_insight` format and that every active archetype's insight appears exactly once across all themes. On mismatch, retry with explicit ID format instruction before falling back to Jaccard.

**Refined Prompt Template**:

```
Analyze these {N} insights and group them by thematic similarity.
Return ONLY valid JSON matching the schema below. No preamble or explanation.

Schema: {json_schema}

Rules:
- Each insight must appear in exactly one theme
- Insights with no thematic overlap should be in a single-member theme
- Use the exact insight_ids provided

Insights:
{formatted_insights}
```

**Fallback Trigger Conditions**:

The following conditions trigger a fallback from LLM-based grouping to the Jaccard keyword overlap method:

1. **Invalid JSON**: LLM returns invalid JSON after 2 retries -> Jaccard fallback
2. **Degenerate grouping (collapsed)**: 0 themes returned, or all insights placed in a single theme -> Jaccard fallback
3. **Latency exceeded**: LLM response takes > 30 seconds -> Jaccard fallback
4. **Model unavailable**: API error or service unavailability -> Jaccard fallback
5. **Degenerate grouping (fragmented)**: All themes are singletons (each containing exactly 1 insight, no grouping occurred) -> Jaccard fallback

On any fallback: log `grouping_method: jaccard_fallback` in `workflow-state.yaml` under the `stage_3` entry. This enables post-hoc analysis of grouping reliability.

**Fallback Escalation Chain**:

Current (Task tool inherits model):
```
Orchestrator model inline (primary)
  -> retry with orchestrator model (1x, 5-second delay)
  -> Jaccard keyword overlap (deterministic)
```

Target (when Task tool supports model selection):
```
Haiku 4.5 (primary)
  -> retry Haiku 4.5 (1x, 5-second delay)
  -> Sonnet 4.5 (1x escalation)
  -> Jaccard keyword overlap (deterministic)
```

### Step 3: Convergence Scoring

For each theme group:

```python
def calculate_convergence_score(group):
    convergence_count = len(group)

    # Convergence multipliers (research-backed)
    if convergence_count >= 4:
        multiplier = 2.5
    elif convergence_count == 3:
        multiplier = 2.0
    elif convergence_count == 2:
        multiplier = 1.5
    else:
        multiplier = 1.0  # Unique insight, no bonus (but not penalized)

    # Average confidence across contributing agents
    avg_confidence = mean([i.confidence for i in group])

    # Research backing bonus
    research_bonus = 1.0 + (0.1 * sum(1 for i in group if i.research_backed))

    # Final score
    return avg_confidence * multiplier * research_bonus
```

### Step 4: Adjusted Thresholds for Partial Completion

When only 4 agents complete (1 missing):

```python
# Reduced convergence multipliers
if convergence_count >= 3:
    multiplier = 2.0  # (was 2.5 for 4+)
elif convergence_count == 2:
    multiplier = 1.3  # (was 1.5)
else:
    multiplier = 1.0
```

### Step 5: Output Generation

**Convergent Insights** (ranked by score):
```yaml
convergent_insights:
  - theme: "Theme description"
    score: 8.5
    convergence_count: 3
    contributing_archetypes: [optimist, analyst, pragmatist]
    synthesized_insight: "Combined insight statement"
    evidence:
      - "Evidence from optimist"
      - "Evidence from analyst"
      - "Evidence from pragmatist"
```

**Divergent Insights** (attributed to source):
```yaml
divergent_insights:
  - archetype: innovator
    insight: "Unique insight from innovator"
    confidence: 4
    evidence: ["Supporting point"]
    note: "This perspective was not echoed by other archetypes"
```

## Confidence Validation

```python
def validate_confidence(raw_confidence, archetype):
    validation_warnings = []

    # Clamp to valid range
    if raw_confidence < 1:
        confidence = 1
        validation_warnings.append(f"{archetype}: Confidence below minimum, set to 1")
    elif raw_confidence > 5:
        confidence = 5
        validation_warnings.append(f"{archetype}: Confidence above maximum, set to 5")
    else:
        confidence = raw_confidence

    return confidence, validation_warnings
```

**Warning Propagation**: Validation warnings should be:
1. Logged in `workflow-state.yaml` under the agent's entry
2. Included in the synthesis output if any clamping occurred
3. Surfaced to user in Stage 4 summary if significant

## Research Backing Impact

| Search Result | Confidence Adjustment |
|---------------|----------------------|
| 2 successful searches | No adjustment |
| 1 successful search | No adjustment |
| 0 successful searches | -1 confidence (floor of 1) |
| Service unavailable | -1 confidence, note "limited research" |

## No Convergence Handling

When zero themes appear in 2+ perspectives:

1. **Diagnose the situation**:
   - All high confidence + no overlap = Genuinely complex problem
   - All low confidence + no overlap = Prompt may need refinement
   - Mixed = Normal divergence

2. **Reframe output**:
   ```markdown
   ## Portfolio of Options

   No strong consensus emerged across perspectives. This suggests:
   - The problem may be genuinely complex with no clear answer
   - Different framings lead to different conclusions
   - Multiple valid approaches exist

   Below are 5 distinct perspectives, each offering a different path forward.
   ```

3. **Offer user choice**:
   - Proceed with portfolio synthesis
   - Refine the original prompt for another attempt
   - Select most relevant perspective for deeper exploration

## Example Calculation

Given 5 perspectives with insights:

| Archetype | Insight Theme | Confidence | Research |
|-----------|---------------|------------|----------|
| Optimist | Growth opportunity | 4 | Yes |
| Critic | Regulatory risk | 5 | Yes |
| Analyst | Market data supports | 4 | Yes |
| Innovator | Platform approach | 3 | No |
| Pragmatist | Growth opportunity | 4 | Yes |

**Grouping**:
- Group A: "Growth opportunity" (Optimist, Pragmatist) - count: 2
- Group B: "Regulatory risk" (Critic) - count: 1
- Group C: "Market data" (Analyst) - count: 1
- Group D: "Platform approach" (Innovator) - count: 1

**Scoring**:
- Group A: avg_conf=4, multiplier=1.5, research=1.2 -> score = 4 * 1.5 * 1.2 = **7.2**
- Group B: avg_conf=5, multiplier=1.0, research=1.1 -> score = 5 * 1.0 * 1.1 = **5.5**
- Group C: avg_conf=4, multiplier=1.0, research=1.1 -> score = 4 * 1.0 * 1.1 = **4.4**
- Group D: avg_conf=3, multiplier=1.0, research=1.0 -> score = 3 * 1.0 * 1.0 = **3.0**

**Output**:
1. Convergent: "Growth opportunity" (score 7.2, 2 archetypes)
2. Divergent: "Regulatory risk" (Critic, confidence 5)
3. Divergent: "Market data supports" (Analyst, confidence 4)
4. Divergent: "Platform approach" (Innovator, confidence 3)
