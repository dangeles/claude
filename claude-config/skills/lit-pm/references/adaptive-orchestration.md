# Adaptive Orchestration

Complexity detection, checkpoint plan generation, and the depth profile system for lit-pm.

## Contents

- Complexity Detection
- Classification Matrix
- Checkpoint Plan Table
- Depth Profile System
- User Presentation Format
- User Override
- Mid-Pipeline Adaptation

---

## Complexity Detection

Given a refined scope from requirements-analyst, three dimensions set the checkpoint plan.

### Dimension 1: Scope

| Indicator | Simple | Medium | Complex |
|-----------|--------|--------|---------|
| Paper count estimate | <10 | 10-30 | 30+ |
| Topic breadth | Single topic | 2-3 themes | Cross-domain |
| Literature maturity | Established field | Emerging field | Contradictory |

Score each indicator 0 (Simple), 1 (Medium), or 2 (Complex) and sum: 0-1 = Simple, 2-3 = Medium, 4+ = Complex.

### Dimension 2: Stakes

| Stakes Level | Keywords | Context signals |
|--------------|----------|-----------------|
| Low | "quick survey", "what's known about", "explore", "overview" | Time constraint mentioned, exploratory language |
| Medium | "inform decision", "compare approaches", "assess feasibility", "evaluate" | Decision mentioned, options being compared |
| High | "grant background", "comprehensive review", "publication", "proposal" | External audience, deliverable mentioned, dollar amounts |

### Dimension 3: User Hints

Explicit flags: `--review-outline` adds a checkpoint at Stage 3, `--review-drafts` adds one after Stage 5, `--full-auto` skips all optional checkpoints. Time words shift the plan too — "quick" or "fast" reduces checkpoints, "thorough" or "comprehensive" adds them. A user who has iterated on similar documents before needs fewer checkpoints than a first-time user.

---

## Classification Matrix

| Scope Score | Low Stakes | Medium Stakes | High Stakes |
|-------------|------------|---------------|-------------|
| Simple (0-1) | Simple | Medium | High-Stakes |
| Medium (2-3) | Medium | Medium | High-Stakes |
| Complex (4+) | Complex | Complex | High-Stakes |

Record the resulting tier with a confidence value and a one-sentence rationale.

Worked example — "Comprehensive review of mRNA delivery for tissue engineering; this will inform our $2M grant proposal" scores complex on paper count (~35) and breadth (cross-domain), medium on maturity (emerging), and high on stakes ("comprehensive", "grant proposal", "$2M"), classifying as High-Stakes: checkpoints at scope, review discovery, outline, and synthesis. By contrast "Quick survey of CAR-T manufacturing challenges" scores 0 on scope with low stakes, classifying as Simple: scope approval, then automatic to delivery.

---

## Checkpoint Plan Table

| Complexity | Stage 0 | Stage 1 | Stage 2 | Stage 3 | Stage 6c | Stage 7 | Stage 7.5 | Rationale |
|------------|---------|---------|---------|---------|----------|---------|-----------|-----------|
| Simple | Auto | CHECKPOINT | Auto | Auto | ACTIVE | Auto | Conditional | Scope approval sufficient |
| Medium | Auto | CHECKPOINT | Auto | CHECKPOINT | ACTIVE | Auto | Conditional | Direction check before heavy lifting |
| Complex | Auto | CHECKPOINT | Auto | CHECKPOINT | ACTIVE | CHECKPOINT | Conditional | Multiple approval points |
| High-Stakes | Auto | CHECKPOINT | CHECKPOINT | CHECKPOINT | ACTIVE | CHECKPOINT | ACTIVE | Maximum oversight |

- **CHECKPOINT**: user approval before proceeding
- **Auto**: runs without user interaction
- **ACTIVE**: always runs as a quality gate, not a user approval gate
- **Conditional**: triggers when the condition is met (>=20% additions for 7.5)

Stage 0 (archival guidelines) always runs automatically and is never a checkpoint.

| Stage | Checkpoint purpose |
|-------|-------------------|
| 1 | Validate research question and boundaries |
| 2 | Confirm foundational reviews are correct |
| 3 | Approve structure before section writing (expensive) |
| 6c | Challenge argument quality, test assumptions per section |
| 7 | Review synthesis before final polish |
| 7.5 | Strategic coherence review when synthesis adds significant content |

---

## Depth Profile System

The depth profile extends the tier classification to calibrate output depth, length, and detail for all downstream specialists. It is generated alongside the checkpoint plan and propagates through every specialist handoff.

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

**Default for absent profile**: If no depth_profile reaches a specialist (e.g. standalone invocation), it falls back to High-Stakes/COMPREHENSIVE defaults.

### Density Guidance (Behavioral Translations)

These values for `writing.density_guidance` give specialists explicit behavioral instructions so the density enum is interpreted consistently:

**DENSE** (Simple and Medium tiers):
> State conclusions directly with 1-2 supporting citations. Do not narrate each paper's methodology individually. Merge findings by theme, not by paper. Eliminate hedging phrases ("it may be argued that", "one could speculate"). Prefer active voice. No paragraph should restate a point made in the preceding paragraph. The goal is integration of information, not enumeration.

**STANDARD** (Complex tier):
> Balance coverage and readability. Each major finding warrants 2-3 sentences. Methodological context is appropriate when it affects interpretation. Transitions are encouraged. Let content determine length within quality gate ranges.

**EXPANSIVE** (High-Stakes tier):
> Cover methodological nuances, contradictions, and subtleties. Extended discussion of limitations is expected. Transitional analysis connecting sections is encouraged. This is a comprehensive reference document.

---

## User Presentation Format

Present the checkpoint plan and depth profile together:

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
1. **Accept** — proceed with this profile and plan
2. **More depth** — upgrade to next tier's profile
3. **Less depth** — use minimal depth (all FOCUSED/BRIEF dimensions)
4. **Add checkpoints** — name the stages to add
5. **Remove checkpoints** — name the stages to skip (Stage 1 stays)
6. **Full-auto** — only pause at Stage 1
```

---

## User Override

| User Response | Action |
|---------------|--------|
| "Accept" / "Looks good" / "Proceed" | Use proposed plan |
| "Add Stage X" | Add checkpoint at specified stage |
| "Skip Stage X" | Remove checkpoint (except Stage 1) |
| "Add all checkpoints" | Use High-Stakes plan |
| "--full-auto" | Only Stage 1 checkpoint |

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

The quality floor holds regardless of overrides — see the Quality gates section of SKILL.md.

---

## Mid-Pipeline Adaptation

Some conditions suggest the complexity was misjudged:

| Condition | Detection | Significance |
|-----------|-----------|--------------|
| Highly contradictory literature | >50% of reviews disagree on key claims | Higher complexity than estimated |
| Large outline | >5 sections identified | More work than expected |
| Stage 6c uncertainty | >40% of sections PASS_WITH_UNCERTAINTY | Arguments weaker than expected; consider a Stage 7 checkpoint |
| Stage 7.5 coherence | Thesis coherence marked WEAK | Document may need restructuring; offer proceed-with-warnings, return to outline, or abort |

On detecting one: log it ("Detected higher complexity than estimated. Reason: {reason}"), offer the user a checkpoint before the upcoming stage, record the decision in `checkpoint_state.modifications[]`, and show the modification history on resume.

Constraints: checkpoints can be added mid-pipeline but not removed; no adaptation after Stage 7 begins; at most 2 adaptations per workflow. When the tier changes, update `depth_profile` to match and communicate it to specialists not yet invoked. Already-completed sections are not redone, but note the inconsistency in session state for the synthesizer.
