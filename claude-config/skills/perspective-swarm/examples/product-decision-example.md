# Example: Product Decision Walkthrough

This example demonstrates a complete perspective-swarm workflow for a product decision.

## User Input

```
Should we add AI-powered features to our productivity app, or focus on
improving core performance and stability?
```

## Stage 1: Problem Framing

**Input Validation**: PASS (78 characters, clear question, no prohibited content)

**Problem Type Identified**: decision

**Reframed Challenge**:
> Evaluate product direction trade-off: AI feature investment versus core performance optimization for productivity application

**Agent Prompts Generated**: 5 prompts tailored to each archetype

**Session Created**: `/tmp/swarm-session-20260204-140000-x7y8z9ab/`

**Session Lock Created**: `.session.lock` with workflow_id and timestamp

**Quality Gate**: PASS

---

## Stage 2: Parallel Perspective Generation

### Optimist Perspective

**Key Insight**: AI features could differentiate us in a crowded market and attract premium customers willing to pay more for intelligent automation.

**Supporting Evidence**:
- Market research shows 67% of productivity app users interested in AI features (Source: Gartner 2025)
- Competitors with AI features command 40% price premium
- AI feature adoption in enterprise productivity tools grew 3x in 2025

**Confidence**: 4/5

**Blind Spots**: May underweight the development cost and timeline; assumes AI features will work well on first release.

---

### Critic Perspective

**Key Insight**: Rushing AI features risks introducing instability when our core product already has performance complaints, potentially accelerating churn.

**Supporting Evidence**:
- Our NPS dropped 12 points in Q4, with "slow performance" as top complaint
- 3 competitors failed after poorly-executed AI launches (ProductHunt analysis)
- AI feature bugs take 3x longer to fix than traditional bugs (industry benchmark)

**Confidence**: 5/5

**Blind Spots**: May be overly focused on current complaints rather than future market position.

---

### Analyst Perspective

**Key Insight**: Data suggests a sequenced approach: performance improvements take 3 months, AI features take 6 months, but doing AI first delays both by 9+ months.

**Supporting Evidence**:
- Internal metrics: 23% of users hit performance issues weekly
- Development estimate: AI features require 2 engineers for 6 months
- Historical data: Our last major feature delayed performance work by 4 months

**Confidence**: 4/5

**Blind Spots**: Limited external data on competitor timelines; internal estimates may be optimistic.

---

### Innovator Perspective

**Key Insight**: Consider a "performance-first AI" approach: use AI to automatically optimize app performance, solving both problems simultaneously.

**Supporting Evidence**:
- Microsoft's Edge browser used AI for memory optimization (shipped 2024)
- Emerging pattern: "AI for efficiency" rather than "AI as feature"
- Could be implemented with smaller team than full AI features

**Confidence**: 3/5

**Blind Spots**: This is a novel approach; may be harder to market than visible AI features.

---

### Pragmatist Perspective

**Key Insight**: We can do performance in Q1, basic AI in Q2, with current team - but only if we don't add other features and reduce technical debt work.

**Supporting Evidence**:
- Team capacity: 4 engineers, 2 allocated to maintenance
- Performance work is well-scoped; AI features are not
- Customers have waited 2 quarters for performance fixes

**Confidence**: 4/5

**Blind Spots**: May be underestimating stakeholder pressure for AI features; doesn't account for competitive moves.

---

**Stage 2 Summary**:
- Agents completed: 5/5
- Time elapsed: 12 minutes
- All quality gates passed

---

## Stage 3: Convergence Analysis

### Convergent Insights (2+ archetypes)

1. **Performance problems are urgent** (Score: 7.8)
   - Contributing: Critic (5), Analyst (4), Pragmatist (4)
   - Multiplier: 2.0x (3 agents)
   - Summary: Current performance issues are causing measurable user dissatisfaction and should be addressed before adding complexity.

2. **Sequencing matters** (Score: 6.0)
   - Contributing: Analyst (4), Pragmatist (4)
   - Multiplier: 1.5x (2 agents)
   - Summary: The order of work significantly impacts total timeline; performance-first enables faster AI delivery.

### Divergent Insights (unique)

1. **Market differentiation opportunity** (Optimist, confidence 4)
   - AI features could command premium pricing and differentiate in market

2. **AI execution risk** (Critic, confidence 5)
   - Poorly-executed AI features could accelerate churn

3. **Performance-first AI** (Innovator, confidence 3)
   - Novel approach using AI for performance optimization rather than visible features

### Conflicts Identified

| Tension | Position A | Position B | Resolution Approach |
|---------|-----------|-----------|---------------------|
| Timing | Ship AI fast for market (Optimist) | Fix performance first (Critic, Analyst, Pragmatist) | Convergence favors performance-first |
| AI approach | Visible AI features (Optimist) | AI for optimization (Innovator) | Present as options, not resolved |

### Key Uncertainties

- True timeline for AI feature development (estimates vary)
- Competitor AI launch timelines (unknown)
- Whether "performance-first AI" is technically feasible for our stack

### Aggregated Blind Spots

- May be underweighting competitive pressure for AI features
- Internal estimates may be optimistic
- Novel approaches (Innovator) are untested

---

## Stage 4: User Review

### Synthesis Presented

**Executive Summary**:
Three of five perspectives converge on addressing performance issues before AI features, with strong evidence of current user dissatisfaction. However, the market opportunity for AI features is real, suggesting a sequenced approach rather than choosing one over the other. An innovative "performance-first AI" option was identified but requires feasibility assessment.

**Recommended Next Steps**:
1. Commit to performance improvements in Q1
2. Scope minimal viable AI features for Q2
3. Investigate "AI for optimization" as potential quick win
4. Monitor competitor AI launches for timing adjustment

### User Decision

User selected: **(A) Accept**

Workflow state: COMPLETED

---

## Output Files

```
/tmp/swarm-session-20260204-140000-x7y8z9ab/
├── .session.lock               # Session lock file
├── workflow-state.yaml         # Final state: COMPLETED
├── stage-1-framing.yaml        # Problem framing
├── perspectives/
│   ├── optimist.md
│   ├── critic.md
│   ├── analyst.md
│   ├── innovator.md
│   └── pragmatist.md
└── stage-3-synthesis.md        # Final synthesis document
```

---

## Key Observations

1. **Convergence revealed priority**: 3/5 perspectives independently identified performance as urgent
2. **Divergent insight was valuable**: Innovator's "performance-first AI" was non-obvious and potentially best option
3. **Conflict was informative**: Tension between Optimist and Critic highlighted real strategic trade-off
4. **Blind spots aggregated well**: Combined view revealed overconfidence in estimates
5. **Total time**: 22 minutes (within 30-minute target)

---

## Alternate Stage 4: User Review with Multiple Handoff Options

This section demonstrates the generic handoff mechanism when the user chooses to continue with another workflow.

### Synthesis Presented with Handoff Options

```
## Synthesis Complete

European market expansion presents a significant growth opportunity with
regulatory complexity as the primary risk factor. Five perspectives were
analyzed, with strong convergence on market potential and divergent views
on entry strategy.

**Core Options:**
(A) Accept - Workflow complete
(B) Refine - Provide feedback, return to Stage 3

**Continue with another workflow:**
(C) [research] lit-pm - Comprehensive literature review (4-24 hours)
    Relevance: HIGH (strategic problem + 3 uncertainties)
(D) [creative] pov-expansion - Cross-domain perspective analysis (2-3 hours)
    Relevance: MEDIUM (strategic problem)
(E) [implementation] programming-pm - Software implementation (2-8 hours)
    Relevance: LOW (no implementation signals)

Select an option: C
```

### User Selects lit-pm (Option C)

**Handoff payload generated**: `/tmp/swarm-session-20260204-183000-a1b2c3d4/handoff-payload.yaml`

```yaml
handoff:
  version: "2.0"
  timestamp: "2026-02-04T19:00:00Z"
  expires_at: "2026-02-04T20:00:00Z"

  source:
    skill: perspective-swarm
    workflow_id: swarm-session-20260204-183000-a1b2c3d4
    session_path: /tmp/swarm-session-20260204-183000-a1b2c3d4/

  target:
    skill: lit-pm
    invocation: "--handoff {payload_path}"
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
      - "Limited quantitative data on our specific vertical (Analyst)"

  research_seeds:
    suggested_terms:
      - term: "B2B SaaS GDPR compliance cost studies"
        rationale: "Address critic's concern about underestimated compliance"
      - term: "European market entry partnership vs direct"
        rationale: "Explore innovator's partnership suggestion"
      - term: "SaaS expansion case studies Europe 2024-2025"
        rationale: "Find comparable examples for analyst validation"
    open_questions:
      - "What is the average GDPR compliance cost for similar SaaS products?"
      - "What partnership models have succeeded for US B2B SaaS in EU?"
      - "How long does typical EU market entry take for our category?"

  meta:
    perspectives_completed: 5
    convergence_level: medium
    user_feedback: "Particularly interested in GDPR compliance depth"
    handoff_reason: "User selected deep literature research for regulatory clarity"
    handoff_chain: ["perspective-swarm"]
    payload_hash: "sha256:a1b2c3d4e5f6..."
    payload_size_bytes: 3247
```

**Invocation**: `/lit-pm --handoff /tmp/swarm-session-20260204-183000-a1b2c3d4/handoff-payload.yaml`

**State transition**: AWAITING_USER -> COMPLETED
