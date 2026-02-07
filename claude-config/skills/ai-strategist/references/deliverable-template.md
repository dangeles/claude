# Deliverable Template

This reference defines the structure of the final deliverable produced by ai-strategist. The editor in Phase 6 uses this template to produce the polished output.

## Template Structure

```markdown
# AI Tool Landscape Assessment: {Scope Title}

**Date**: {YYYY-MM-DD}
**Scan Type**: {Quarterly / Deep Dive / Event-Triggered}
**Gaps Assessed**: {List of active gaps with weights}
**Tools Evaluated**: {Total count}
**Scoring Methodology**: Weighted composite (Integration {X}% / Gap Coverage {Y}% / Cost {Z}%)

---

## Executive Summary

{2-3 paragraphs summarizing:
- Key findings (what did we learn about the current AI tool landscape?)
- Top recommendations (which 3-5 tools should be prioritized?)
- Critical insights (what surprised us? what's changed since last scan?)
- Action required (what should happen next?)}

---

## Tool Landscape Overview

### MCP Servers & Claude Code Integrations
{Summary of category findings, notable tools, maturity trends}

### AI Frameworks & Agentic Workflows
{Summary of category findings, notable frameworks, pattern trends}

### Scientific & Computational Biology Tools
{Summary of category findings, notable tools, domain-specific insights}

### Community Patterns & Emerging Trends
{Summary of community signals, emerging patterns, practitioner recommendations}

---

## Gap Coverage Analysis

### Gap 1: {Gap Name} (Weight: {X}%)
**Best tool for this gap**: {Per-gap champion, score}
**Coverage summary**: {How well is this gap addressed by the evaluated tools?}
**Top tools**: {2-3 tools that score highest on this gap}
**What changes if we do nothing**: {Consequence of not addressing this gap}

### Gap 2: {Gap Name} (Weight: {X}%)
{Same structure}

{Repeat for all active gaps}

---

## Scored Tool Matrix

| Rank | Tool | Category | Integration (40%) | Gap Coverage (35%) | Cost (25%) | Composite | Radar Ring |
|------|------|----------|-------------------|-------------------|------------|-----------|------------|
| 1 | {Tool} | {Cat} | {0.00} | {0.00} | {0.00} | {0.000} | {Ring} |
| 2 | {Tool} | {Cat} | {0.00} | {0.00} | {0.00} | {0.000} | {Ring} |
{... all scored tools, ranked by composite}

### Per-Gap Champions

| Gap | Champion Tool | Gap Score | Composite Rank |
|-----|--------------|-----------|----------------|
| {Gap 1} | {Tool} | {0.00} | #{N} |
{... one row per gap}

### Sensitivity Analysis Summary

**Weight-Robust Tools** (top 5 in all 4 weight configurations):
- {Tool 1}: Stable across all weightings
- {Tool 2}: Stable across all weightings

**Weight-Sensitive Tools** (rank changes 3+ positions):
- {Tool X}: Ranks #{N} with base weights, #{M} with {config} weights
  Implication: {What this means for the recommendation}

**Top Recommendation Stability**: {Stable / Changed under {config}}

---

## Technology Radar (Optional)

### Adopt (Composite > 0.75)
{Tools ready for immediate production use}

### Trial (Composite 0.50-0.75)
{Tools ready for pilot testing this quarter}

### Assess (Composite 0.25-0.50)
{Tools worth monitoring and exploring}

### Hold (Composite < 0.25)
{Tools to avoid for now}

---

## Integration Roadmap

### Quick Wins (This Week)
| Tool | Gap(s) Addressed | Integration Approach | Effort | Risk |
|------|-----------------|---------------------|--------|------|
| {Tool} | {Gaps} | {How to integrate} | {Low/Med/High} | {Low/Med/High} |

### Short-Term (This Month)
{Same table format}

### Medium-Term (This Quarter)
{Same table format}

### Strategic (Next Quarter)
{Same table format}

---

## Risk Assessment

{From adversarial review -- challenges, concerns, and mitigations}

### Integration Risks
{Specific risks around tool integration feasibility}

### Cost Risks
{Hidden costs, vendor lock-in, sustainability concerns}

### Coverage Risks
{Gaps that remain poorly addressed, missing tool categories}

### Methodology Limitations
{Known limitations of the scoring approach, data gaps}

---

## Methodology

**Scoring Framework**: Weighted composite with 3 dimensions
- Integration Feasibility ({X}%): {Brief description}
- Workflow Gap Coverage ({Y}%): {Brief description}
- Cost/Sustainability ({Z}%): {Brief description}

**Sensitivity Analysis**: Tested 4 weight configurations: Base, Integration-heavy, Gap-focused, Cost-conscious

**Research Agents**: {N} agents covering {categories}

**Tools Evaluated**: {Total} tools, {N} passed to assessment (capped at 30)

**Limitations**: {Data freshness, search coverage gaps, known biases}

---

## Appendix: Tool Details

### {Tool Name}
- **URL**: {url}
- **Category**: {category}
- **License**: {license}
- **Integration Score**: {0.00} -- {justification}
- **Gap Coverage**: {per-gap breakdown}
- **Cost Score**: {0.00} -- {justification}
- **Composite**: {0.000}
- **Key Strengths**: {bullet list}
- **Key Risks**: {bullet list}
- **Next Step**: {specific actionable recommendation}

{Repeat for top 10-15 tools}
```

## Template Usage Notes

1. The editor should fill all `{placeholder}` fields with actual data from Phase 2-5 outputs
2. Sections marked "Optional" (Technology Radar) should be included only if the orchestrator enabled them
3. The executive summary should be written last, after all sections are populated
4. Tool details in the appendix are ordered by composite score (highest first)
5. The risk assessment section incorporates findings from the adversarial review (Phase 5)
6. If the workflow completed in degraded mode, note this in the Methodology section
