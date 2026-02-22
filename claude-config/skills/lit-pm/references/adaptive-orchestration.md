# Adaptive Orchestration

Complexity detection logic and checkpoint plan generation for lit-pm.

---

## Complexity Detection

When lit-pm receives a refined scope from requirements-analyst, it analyzes three dimensions to determine the appropriate checkpoint plan.

### Dimension 1: Scope Indicators

| Indicator | Simple | Medium | Complex |
|-----------|--------|--------|---------|
| Paper count estimate | <10 | 10-30 | 30+ |
| Topic breadth | Single topic | 2-3 themes | Cross-domain |
| Literature maturity | Established field | Emerging field | Contradictory |

**Scoring**:
- Each indicator scores 0 (Simple), 1 (Medium), or 2 (Complex)
- Sum scores: 0-1 = Simple, 2-3 = Medium, 4-5 = Complex, 6 = Complex

### Dimension 2: Stakes Indicators

**Keyword Detection**:

| Stakes Level | Keywords |
|--------------|----------|
| Low | "quick survey", "what's known about", "explore", "overview" |
| Medium | "inform decision", "compare approaches", "assess feasibility", "evaluate" |
| High | "grant background", "comprehensive review", "publication", "proposal" |

**Context Signals**:

| Stakes Level | Signals |
|--------------|---------|
| Low | Time constraint mentioned, exploratory language |
| Medium | Decision mentioned, options being compared |
| High | External audience, deliverable mentioned, dollar amounts |

### Dimension 3: User Hints

**Explicit Flags**:
- `--review-outline`: Add checkpoint at Stage 3
- `--review-drafts`: Add checkpoint after Stage 5
- `--full-auto`: Skip all optional checkpoints

**Time Constraints**:
- "quick" or "fast": Reduce checkpoints
- "thorough" or "detailed": Add checkpoints
- "comprehensive": Maximum checkpoints

**Prior Context**:
- User has iterated on similar documents before: Trust built, can reduce checkpoints
- First time user: Default to more checkpoints

---

## Complexity Classification Algorithm

```yaml
classify_complexity:
  input:
    scope_document: object
    user_prompt: string
    flags: list

  process:
    1_score_scope:
      paper_count: estimate_from_topic(scope)
      breadth: count_themes(scope.in_scope)
      maturity: assess_literature_state(scope.research_question)

    2_detect_stakes:
      keywords: scan_for_keywords(user_prompt)
      context: analyze_context_signals(user_prompt)

    3_check_hints:
      flags: parse_flags(flags)
      time_words: scan_time_constraints(user_prompt)

    4_combine:
      scope_tier: map_score_to_tier(scope_score)
      stakes_boost: high_stakes ? +1 : 0
      hint_adjustment: apply_hint_adjustments()

  output:
    tier: enum  # Simple | Medium | Complex | High-Stakes
    confidence: float
    rationale: string
```

### Classification Matrix

| Scope Score | Low Stakes | Medium Stakes | High Stakes |
|-------------|------------|---------------|-------------|
| Simple (0-1) | Simple | Medium | High-Stakes |
| Medium (2-3) | Medium | Medium | High-Stakes |
| Complex (4+) | Complex | Complex | High-Stakes |

---

## Checkpoint Plan Table

Based on complexity tier, lit-pm proposes this checkpoint plan:

| Complexity | Stage 0 | Stage 1 | Stage 2 | Stage 3 | Stage 6c | Stage 7 | Stage 7.5 | Rationale |
|------------|---------|---------|---------|---------|----------|---------|-----------|-----------|
| Simple | Auto | CHECKPOINT | Auto | Auto | ACTIVE | Auto | Conditional | Scope approval sufficient |
| Medium | Auto | CHECKPOINT | Auto | CHECKPOINT | ACTIVE | Auto | Conditional | Direction check before heavy lifting |
| Complex | Auto | CHECKPOINT | Auto | CHECKPOINT | ACTIVE | CHECKPOINT | Conditional | Multiple approval points |
| High-Stakes | Auto | CHECKPOINT | CHECKPOINT | CHECKPOINT | ACTIVE | CHECKPOINT | ACTIVE | Maximum oversight |

**Checkpoint Types**:
- **CHECKPOINT**: User approval required before proceeding
- **Auto**: Runs automatically, no user interaction
- **ACTIVE**: Always runs as quality gate (not user approval gate)
- **Conditional**: Triggers only when condition met (>=20% additions for 7.5)

**Note**: Stage 0 (Archival Guidelines) always runs automatically and is never a checkpoint.

### Checkpoint Purposes

| Stage | Checkpoint Purpose |
|-------|-------------------|
| 0 | (Auto) Initialize session, extract archival guidelines |
| 1 | Validate research question and boundaries |
| 2 | Confirm foundational reviews are correct |
| 3 | Approve structure before section writing (expensive) |
| 6c | Challenge argument quality, test assumptions per section |
| 7 | Review synthesis before final polish |
| 7.5 | Strategic coherence review when synthesis adds significant content |

---

## Depth Profile System

The depth profile extends the complexity tier classification to calibrate output depth, length, and detail for all downstream specialists. It is generated alongside the checkpoint plan and propagates through every specialist handoff.

### Schema

```yaml
depth_profile:
  tier: string          # mirrors complexity tier: Simple|Medium|Complex|High-Stakes

  # Qualitative directive — the "soul" of the profile, prepended to every specialist prompt
  directive: string

  sections:
    target_count: string      # range: "2-3" | "3-4" | "3-5" | "4-6"
    depth_per_section: string # FOCUSED | STANDARD | COMPREHENSIVE

  research:
    papers_per_section: string  # range: "10-15" | "12-20" | "15-25" | "15-30"
    recency_survey: string      # BRIEF | STANDARD | COMPREHENSIVE
    quantitative_table: string  # OPTIONAL | RECOMMENDED | REQUIRED

  writing:
    density: string             # DENSE | STANDARD | EXPANSIVE
    density_guidance: string    # explicit behavioral translation (see table below)

  synthesis:
    augmentation_budget: string   # minimal | moderate | generous
    introduction_scope: string    # BRIEF | STANDARD | COMPREHENSIVE
    conclusion_scope: string      # BRIEF | STANDARD | COMPREHENSIVE
```

### Tier-to-Depth Mapping

| Dimension | Simple | Medium | Complex | High-Stakes |
|-----------|--------|--------|---------|-------------|
| **directive** | "Prioritize succinctness. Every sentence must earn its place. Prefer one precise sentence over two hedged ones." | "Balance thoroughness with readability. Prefer concise coverage; expand only when complexity demands it." | "Cover thoroughly but avoid redundancy. Add depth where it genuinely matters; trim where it does not." | "Be comprehensive and rigorous. Depth and detail are expected; cover all angles." |
| **sections.target_count** | 2-3 | 3-4 | 3-5 | 4-6 |
| **sections.depth_per_section** | FOCUSED | STANDARD | STANDARD | COMPREHENSIVE |
| **research.papers_per_section** | 10-15 | 12-20 | 15-25 | 15-30 |
| **research.recency_survey** | BRIEF | STANDARD | STANDARD | COMPREHENSIVE |
| **research.quantitative_table** | OPTIONAL | RECOMMENDED | RECOMMENDED | REQUIRED |
| **writing.density** | DENSE | DENSE | STANDARD | EXPANSIVE |
| **synthesis.augmentation_budget** | minimal | moderate | moderate | generous |
| **synthesis.introduction_scope** | BRIEF | STANDARD | STANDARD | COMPREHENSIVE |
| **synthesis.conclusion_scope** | BRIEF | STANDARD | STANDARD | COMPREHENSIVE |

**Default for absent profile**: If no depth_profile is provided to a specialist (e.g., standalone invocation), all specialists fall back to High-Stakes/COMPREHENSIVE defaults — equivalent to current behavior.

### Density Guidance (Behavioral Translations)

These values for `writing.density_guidance` give specialists explicit behavioral instructions to prevent inconsistent interpretation of density enum values:

**DENSE** (Simple and Medium tiers):
> State conclusions directly with 1-2 supporting citations. Do NOT narrate each paper's methodology individually. Merge findings by theme, not by paper. Eliminate hedging phrases ("it may be argued that", "one could speculate"). Prefer active voice. No paragraph should restate a point made in the preceding paragraph. The goal is integration of information, not enumeration.

**STANDARD** (Complex tier):
> Balance coverage and readability. Each major finding warrants 2-3 sentences. Methodological context is appropriate when it affects interpretation. Transitions are encouraged. Let content determine length within quality gate ranges.

**EXPANSIVE** (High-Stakes tier):
> Cover methodological nuances, contradictions, and subtleties. Extended discussion of limitations is expected. Transitional analysis connecting sections is encouraged. This is a comprehensive reference document.

### User Presentation Format

After generating the checkpoint plan, also generate the depth_profile and present both together:

```
ORCHESTRATION PLAN

**Complexity Tier**: [TIER] — [one-sentence rationale]

**Checkpoint Plan**: [checkpoint plan table]

**Depth Profile**: [TIER] — "[directive first sentence]"
- Sections: [target_count] sections (depth: [depth_per_section])
- Research: [papers_per_section] papers/section | recency: [recency_survey] | table: [quantitative_table]
- Writing: [density] density
- Synthesis: [augmentation_budget] augmentation | intro: [introduction_scope] | conclusion: [conclusion_scope]

Override options:
1. **Accept** — proceed with this profile
2. **More depth** — upgrade to next tier's profile
3. **Less depth** — use minimal depth (all FOCUSED/BRIEF dimensions)
4. **Custom** — specify what to change

[Existing checkpoint override options continue here]
```

### Mid-Pipeline Depth Profile Adaptation

When the complexity tier is upgraded or downgraded mid-pipeline (see Mid-Pipeline Adaptation section below), the depth_profile MUST also be updated to match the new tier. The updated profile must be explicitly communicated to any downstream specialists not yet invoked. Already-completed sections are not retroactively changed, but a note should be added to the session state flagging the inconsistency for the synthesizer.

---

## User Override

After proposing checkpoint plan, present to user for approval:

### Presentation Format

```markdown
**Proposed Checkpoint Plan**

Based on your request, I've classified this as a **{tier}** complexity review.

Rationale: {rationale}

Proposed checkpoints:
- [x] Stage 1 (Scope) - ALWAYS
- [{x/?}] Stage 2 (Reviews) - {Yes/No}
- [{x/?}] Stage 3 (Outline) - {Yes/No}
- [{x/?}] Stage 7 (Synthesis) - {Yes/No}

Options:
1. **Accept** - Proceed with this plan
2. **Add checkpoints** - Tell me which stages to add
3. **Remove checkpoints** - Tell me which stages to skip
4. **Full-auto** - Only pause at Stage 1
```

### Override Options

| User Response | Action |
|---------------|--------|
| "Accept" / "Looks good" / "Proceed" | Use proposed plan |
| "Add Stage X" | Add checkpoint at specified stage |
| "Skip Stage X" | Remove checkpoint (except Stage 1) |
| "Add all checkpoints" | Use High-Stakes plan |
| "--full-auto" | Only Stage 1 checkpoint |

### Recording Override

```yaml
checkpoint_state:
  tier_detected: string
  tier_reason: string
  plan_proposed:
    stage_1: true
    stage_2: boolean
    stage_3: boolean
    stage_6c: "ACTIVE"  # Always runs, distinct from checkpoint
    stage_7: boolean
    stage_7_5:          # Conditional, not boolean
      mode: enum  # AUTO | CONDITIONAL | ACTIVE
      trigger_condition: ">=20% OR HIGH-STAKES"
  user_override:
    action: enum  # accept | add | remove | full_auto
    modifications:
      - stage: integer | string  # Allow "6c" and "7.5"
        change: enum  # added | removed
    timestamp: ISO8601
  plan_final:
    stage_1: true
    stage_2: boolean
    stage_3: boolean
    stage_6c: "ACTIVE"
    stage_7: boolean
    stage_7_5:
      mode: enum
      trigger_condition: string
```

---

## Mid-Pipeline Adaptation

During execution, lit-pm monitors for conditions that suggest complexity was misjudged.

### Trigger Conditions

| Condition | Detection | Significance |
|-----------|-----------|--------------|
| Highly contradictory literature | >50% of reviews disagree on key claims | Higher than expected complexity |
| Large outline | >5 sections identified | More work than expected |
| Slow progress | Stages running >50% over estimate | Possible complexity underestimate |
| Fast progress | Stages completing <50% of estimate | Possible complexity overestimate |

### Adaptation Protocol

1. **Detect trigger condition**
2. **Log**: "Detected higher complexity than estimated. Reason: {reason}"
3. **Offer user**: "Would you like to add a checkpoint before Stage {N}? (yes/no)"
4. **Record decision** in checkpoint_state.modifications[]
5. **Show modification history** on resume

### Adaptation Message Format

```markdown
**Complexity Adjustment**

I've detected this review may be more complex than initially estimated.

Reason: {specific reason}

Current checkpoint plan:
- [x] Stage 1 (complete)
- [ ] Stage 3 (upcoming)
- [ ] Stage 7 (upcoming)

Recommendation: Add checkpoint at Stage {N}

Would you like to:
1. Add the recommended checkpoint
2. Add a different checkpoint
3. Continue with current plan
```

### Adaptation Constraints

- Cannot remove checkpoints mid-pipeline (only add)
- Cannot adapt after Stage 7 begins
- Maximum 2 adaptations per workflow (prevent thrashing)
- When tier changes, depth_profile MUST also be updated to match the new tier (see Depth Profile System section above)

### Mid-Pipeline Adaptation: Devil's Advocate Triggers

```yaml
adaptation_triggers_da:
  stage_6c:
    trigger: >40% of sections have PASS_WITH_UNCERTAINTY status
    significance: "Arguments weaker than expected, consider adding Stage 7 checkpoint"
    adaptation_offer: "Would you like to add a checkpoint at Stage 7 to review before synthesis?"

  stage_7_5:
    trigger: Strategic thesis coherence marked "WEAK"
    significance: "Document may need restructuring"
    adaptation_offer: "Synthesis review found thesis coherence issues. Options: proceed with warnings, return to outline, abort"
```

---

## Example Detection Scenarios

### Example 1: Simple

```yaml
user_prompt: "Quick survey of CAR-T manufacturing challenges"

detection:
  scope:
    paper_count: ~8 (simple)
    breadth: 1 domain (simple)
    maturity: established (simple)
    score: 0
  stakes:
    keywords: ["quick", "survey"]
    level: Low
  hints:
    flags: none
    time_words: ["quick"]

classification:
  tier: Simple
  confidence: 0.9
  rationale: "Single-domain, established literature, quick survey requested"

proposed_plan:
  stage_1: true
  stage_2: false
  stage_3: false
  stage_7: false

message: |
  I'll run the full pipeline automatically after scope approval.
  You can add --review-outline if you want to check structure first.
```

### Example 2: High-Stakes

```yaml
user_prompt: |
  Comprehensive review of mRNA delivery for tissue engineering -
  this will inform our $2M grant proposal

detection:
  scope:
    paper_count: ~35 (complex)
    breadth: cross-domain - mRNA + tissue eng (complex)
    maturity: emerging (medium)
    score: 5
  stakes:
    keywords: ["comprehensive", "grant proposal", "$2M"]
    level: High
  hints:
    flags: none
    time_words: ["comprehensive"]

classification:
  tier: High-Stakes
  confidence: 0.95
  rationale: "Cross-domain topic, substantial paper count, high-stakes grant context"

proposed_plan:
  stage_1: true
  stage_2: true
  stage_3: true
  stage_7: true

message: |
  Given the grant context, I'll pause for approval at:
  - Review discovery (ensure we have the right foundational papers)
  - Outline structure (get section framing right before writing)
  - Final synthesis (review before editorial polish)
```

### Example 3: Medium with Override

```yaml
user_prompt: "Literature review on bioreactor oxygenation strategies"

detection:
  scope:
    paper_count: ~20 (medium)
    breadth: 2 themes (medium)
    maturity: established (simple)
    score: 3
  stakes:
    keywords: none detected
    level: Medium (inferred from "literature review")
  hints:
    flags: ["--review-outline"]

classification:
  tier: Medium
  confidence: 0.8

proposed_plan:
  stage_1: true
  stage_2: false
  stage_3: true  # From flag
  stage_7: false

user_response: "Also add Stage 7 checkpoint"

final_plan:
  stage_1: true
  stage_2: false
  stage_3: true
  stage_7: true  # Added by user

checkpoint_state:
  user_override:
    action: add
    modifications:
      - stage: 7
        change: added
```

---

## Quality Floor Protection

Even with `--full-auto`, these checks cannot be skipped:

```yaml
quality_floor:
  stage_2:
    minimum_reviews: 4
    minimum_convergence: 1
    on_failure: "HALT: Review discovery below minimum threshold."

  stage_3:
    minimum_sections: 2
    max_section_imbalance: 50%  # No section > 50% of total
    on_failure: "HALT: Outline severely unbalanced."

  stage_5:
    minimum_papers_per_section: 10
    recency_survey_required: true
    on_failure: "HALT: Section research insufficient."
```

Quality floor violations always escalate to user, regardless of checkpoint plan.
