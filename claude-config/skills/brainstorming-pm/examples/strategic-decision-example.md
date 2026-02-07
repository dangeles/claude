# Strategic Decision Example: European Market Expansion

This walkthrough demonstrates brainstorming-pm orchestrating a complete 4-stage pipeline for a strategic business decision.

## User Prompt

```
Should we pivot our B2B SaaS product to focus on the European market?
```

---

## Stage 1: Framing (Elapsed: 3 min)

**State**: `[Stage 1/4 - FRAMING]`

The orchestrator validates the prompt (38 characters, clear question, no prohibited content) and performs framing inline.

**Output** (`stage-1-framing.yaml`):

```yaml
problem_type: strategic
reframed_challenge: >
  Evaluate the strategic opportunity and risks of pivoting a B2B SaaS product
  to focus on the European market, considering regulatory, competitive,
  financial, and operational dimensions.

archetype_prompts:
  optimist: >
    Analyze the European market expansion opportunity for a B2B SaaS product.
    Focus on growth potential, market gaps, and competitive advantages.
    Conduct 1-2 web searches for supporting data.
  critic: >
    Critically evaluate the risks and challenges of pivoting to the European
    market. Consider regulatory barriers (GDPR), cost overruns, and
    competitive threats. Conduct 1-2 web searches for evidence.
  analyst: >
    Provide a data-driven assessment of the European B2B SaaS market.
    Analyze market size, growth rates, customer acquisition costs, and
    unit economics. Conduct 1-2 web searches for current data.
  innovator: >
    Propose unconventional approaches to European market entry that go
    beyond direct expansion. Consider partnerships, platform models,
    or vertical-specific strategies. Conduct 1-2 web searches.
  pragmatist: >
    Assess the operational feasibility of European expansion. Consider
    timeline, resource requirements, team structure, and incremental
    vs. full-pivot approaches. Conduct 1-2 web searches.
```

**Optional checkpoint**: Reframed challenge presented to user. User approves and continues.

---

## Stage 2: Parallel Perspective Generation (Elapsed: 15 min total, 12 min for this stage)

**State**: `[Stage 2/4 - DIVERGING]`

The orchestrator launches 5 Task tool invocations concurrently.

### Status Board

```
[Stage 2/4 - DIVERGING] 4 of 5 agents complete (1 timed out)

| Agent      | Status   | Elapsed | Confidence |
|------------|----------|---------|------------|
| Optimist   | complete | 6:30    | 4          |
| Critic     | complete | 8:45    | 5          |
| Analyst    | complete | 7:15    | 4          |
| Innovator  | timeout  | 10:00   | --         |
| Pragmatist | complete | 9:20    | 4          |
```

### Abbreviated Perspective Outputs

**Optimist** (perspectives/optimist.md):
> **Key Insight**: The European B2B SaaS market is projected to grow at 12% CAGR through 2028, with significant whitespace in vertical-specific solutions. Early movers in underserved verticals (construction tech, agricultural SaaS) command premium pricing.
>
> **Evidence**: (1) Gartner 2025 EU SaaS report shows 12.3% CAGR; (2) Low penetration in mid-market segment compared to US.
>
> **Confidence**: 4/5
> **Blind spots**: May underweight execution complexity and GDPR compliance costs.

**Critic** (perspectives/critic.md):
> **Key Insight**: GDPR compliance costs are systematically underestimated in expansion plans. A 2025 McKinsey survey found that 67% of US SaaS companies exceeded their GDPR compliance budget by 40-120%. Additionally, EU customers have longer sales cycles (avg 6.2 months vs 3.8 months in US).
>
> **Evidence**: (1) McKinsey 2025 survey of 200 SaaS companies; (2) Salesforce EU sales cycle benchmarks.
>
> **Confidence**: 5/5
> **Blind spots**: May undervalue market potential by overweighting short-term costs.

**Analyst** (perspectives/analyst.md):
> **Key Insight**: Unit economics analysis suggests breakeven at 18-24 months in EU expansion, assuming 40% higher CAC (customer acquisition cost) and 15% lower initial ARPU compared to US market. The break-even horizon extends to 30+ months if GDPR compliance costs exceed baseline by >60%.
>
> **Evidence**: (1) EU SaaS market size: $85B (2025), projected $142B (2028); (2) Industry CAC benchmarks for EU entry.
>
> **Confidence**: 4/5
> **Blind spots**: Limited data on vertical-specific economics; aggregated data may mask opportunity.

### User Options (1 agent failed)

```
4 of 5 perspectives complete. 1 failed: Innovator (timeout after 10:00)

Options:
(A) Retry Innovator (estimated +5 min)
(B) Proceed with 4 perspectives (reduced diversity)
(C) Abort workflow
```

User selects **(B) Proceed with 4 perspectives**. This option is available because N=4 >= 2, meeting the minimum threshold for meaningful convergence analysis.

---

## Stage 3: Convergence + Parallel Discovery (Elapsed: 21 min total, 6 min for this stage)

**State**: `[Stage 3/4 - CONVERGING]`

### Track A: Convergence Analysis (inline)

The orchestrator collects all 4 perspective outputs and performs LLM-based grouping.

**LLM Grouping Result** (grouping_method: `llm_primary`, latency: 1250ms):

```json
{
  "themes": [
    {
      "theme_id": "market_opportunity",
      "theme_description": "European B2B SaaS market presents significant growth opportunity with underserved segments",
      "insight_ids": ["optimist_insight", "analyst_insight"]
    },
    {
      "theme_id": "regulatory_cost_risk",
      "theme_description": "GDPR and regulatory compliance costs are a primary risk factor likely to be underestimated",
      "insight_ids": ["critic_insight"]
    },
    {
      "theme_id": "operational_feasibility",
      "theme_description": "Incremental phased approach reduces risk but extends timeline",
      "insight_ids": ["pragmatist_insight"]
    }
  ]
}
```

**Convergence Scoring**:

| Theme | Count | Avg Confidence | Multiplier | Research Bonus | Score |
|-------|-------|---------------|------------|----------------|-------|
| Market Opportunity | 2 | 4.0 | 1.5 | 1.2 | **7.2** |
| Regulatory Cost Risk | 1 | 5.0 | 1.0 | 1.1 | **5.5** |
| Operational Feasibility | 1 | 4.0 | 1.0 | 1.1 | **4.4** |

**Missing archetype compensation** (Innovator timed out): "Note: potentially conservative recommendations. The Innovator perspective was unavailable, which may result in underrepresentation of unconventional approaches."

**Synthesis output** (`stage-3-synthesis.md`):

**Convergent insights** (2):
1. *Market Opportunity* (score 7.2, Optimist + Analyst): EU B2B SaaS market growing at 12% CAGR with whitespace in verticals. Breakeven at 18-24 months under baseline assumptions.
2. (Below threshold -- no second convergent theme with 4 agents)

**Divergent insights** (3):
1. *Regulatory Cost Risk* (Critic, confidence 5): GDPR costs systematically underestimated; 67% exceed budget by 40-120%
2. *Operational Feasibility* (Pragmatist, confidence 4): Phased approach recommended; full pivot carries high execution risk
3. *Missing Innovator perspective*: Noted as gap; unconventional strategies not explored

### Track B: Workflow Discovery (parallel)

Discovery ran in parallel with convergence. Completed in 1.2 seconds.

**Discovered workflows** (filtered -- `brainstorming-pm` and `perspective-swarm` excluded):

```yaml
workflows:
  - skill: lit-pm
    categories: [research, literature]
    description: "Comprehensive literature review (4-24 hours)"
    scope: global
  - skill: programming-pm
    categories: [implementation]
    description: "Software implementation (2-8 hours)"
    scope: global
```

---

## Stage 4: User Review + Handoff (Elapsed: 21+ min, user-controlled)

**State**: `[Stage 4/4 - OUTPUT]`

Relevance scoring applied against synthesis content:

| Workflow | Score | Rationale |
|----------|-------|-----------|
| lit-pm | 7 | Category match (research+strategic: +3), high uncertainty (+2), low convergence (+2) |
| programming-pm | 2 | Implementation keyword partial match (+2) |

**Presented to user**:

```
[Stage 4/4 - OUTPUT] Synthesis Complete

## Executive Summary

European market expansion presents a significant growth opportunity
(12% CAGR, $85B -> $142B by 2028) with underserved vertical segments.
However, GDPR compliance costs are systematically underestimated --
67% of comparable US SaaS companies exceeded budget by 40-120%.
Unit economics suggest 18-24 month breakeven under baseline assumptions,
extending to 30+ months if compliance costs overrun by >60%.

Note: Innovator perspective unavailable (timeout). Unconventional
approaches (partnerships, platform models) were not fully explored.

**Convergent**: Market opportunity validated by Optimist and Analyst
**Key Risk**: Regulatory cost underestimation (Critic, confidence 5/5)
**Feasibility**: Phased approach recommended over full pivot (Pragmatist)

Options:
(A) Accept - Workflow complete
(B) Refine - Provide feedback, return to Stage 3

Continue with another workflow:
(C) [research] lit-pm - Comprehensive literature review (4-24 hours)
    Relevance: HIGH - 3 key uncertainties identified (GDPR costs,
    competitor response, currency impact)
(D) [implementation] programming-pm - Software implementation (2-8 hours)
    Relevance: LOW

Select an option:
```

### Timing Summary

| Stage | Duration | Notes |
|-------|----------|-------|
| Stage 1: Framing | 3 min | Inline, includes user checkpoint |
| Stage 2: Perspectives | 12 min | 4/5 agents completed; 1 timeout |
| Stage 3: Convergence | 6 min | Track A: 5.8 min; Track B: 1.2 sec (parallel) |
| Stage 4: User Review | user-controlled | Synthesis presented with options |
| **Total (to Stage 4)** | **21 min** | Within 45-min safety ceiling |
