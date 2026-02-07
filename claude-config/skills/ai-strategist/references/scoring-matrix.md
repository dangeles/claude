# Weighted Scoring Matrix

This reference defines the scoring framework used by the strategist in Phase 3 to evaluate AI tools against workflow gaps.

## Scoring Dimensions

### Integration Feasibility (40%)

Measures how readily a tool integrates with existing infrastructure (Claude Code, MCP architecture, Slack, Notion, LinkedIn, X).

| Score | Level | Criteria |
|---|---|---|
| 0.0 | No integration path | No API, no MCP support, no documented integration |
| 0.25 | Theoretical | Integration path exists in theory but no documentation or community examples |
| 0.5 | Documented but untested | Integration documented by vendor or community, but no production usage reports |
| 0.75 | Community-proven | Active GitHub repo with >100 stars, documented issues resolved, >1 year maintenance, production usage reports |
| 1.0 | Official support | First-party MCP server, Anthropic-endorsed, or native Claude Code integration |

Note: Each level implicitly captures maturity -- a tool cannot score 0.75+ without demonstrated maturity and maintenance track record.

### Workflow Gap Coverage (35%)

Measures how well a tool addresses identified workflow gaps. Scored per gap, then aggregated.

| Score | Level | Criteria |
|---|---|---|
| 0.0 | Does not address | No relevance to this gap |
| 0.25 | Tangentially related | Touches the problem space but does not directly solve it |
| 0.5 | Partially addresses | Covers some use cases within the gap |
| 0.75 | Substantially addresses | Covers most use cases, proven in similar contexts |
| 1.0 | Fully solves | Comprehensive coverage, proven in comparable workflows |

**Gap Coverage Aggregation**: Weighted average across all active gaps. Gap weights are set during Phase 1 scope refinement. Default: equal weights across all active gaps.

### Cost/Sustainability (25%)

Measures financial accessibility and long-term viability for early-stage startups.

| Score | Level | Criteria |
|---|---|---|
| 0.0 | Enterprise-only | >$500/month, enterprise sales process |
| 0.25 | Expensive but available | $100-500/month, self-serve |
| 0.5 | Affordable | $20-100/month or limited free tier |
| 0.75 | Generous free tier | Free tier sufficient for startup use, or <$20/month |
| 1.0 | Open-source | Active community maintenance, no vendor dependency |

## Composite Score Formula

```
composite = 0.40 * integration + 0.35 * gap_coverage_weighted_avg + 0.25 * cost
```

Where `gap_coverage_weighted_avg` = sum(gap_weight_i * gap_score_i) / sum(gap_weight_i) for all active gaps.

## Sensitivity Analysis Protocol

After producing the base scored matrix (40/35/25 weights), recompute with three alternate configurations:

| Configuration | Integration | Gap Coverage | Cost | Bias |
|---|---|---|---|---|
| Base | 40% | 35% | 25% | Balanced |
| Integration-heavy | 50% | 25% | 25% | Favors easy adoption |
| Gap-focused | 30% | 45% | 25% | Favors problem-solving |
| Cost-conscious | 30% | 35% | 35% | Favors budget constraints |

**Analysis Outputs**:

1. **Weight-robust tools**: Tools that rank in the top 5 across all 4 weight configurations. These are safe recommendations regardless of prioritization preference.

2. **Weight-sensitive tools**: Tools whose rank changes by 3 or more positions across configurations. Flag these in the roadmap with explanation. Example: "Tool X ranks #2 with base weights but drops to #6 under gap-focused weighting."

3. **Top recommendation stability**: If the #1 recommended tool changes across configurations, escalate to user for weight validation before finalizing roadmap.

## Per-Gap Champion Analysis

For each active workflow gap, identify the top-scoring tool on that gap's dimension regardless of composite rank.

Present alongside composite rankings to prevent masking per-gap excellence. A tool that scores 1.0 on a critical gap but low on cost might be the right choice for that specific gap even if its composite rank is low.

## Differentiation Check

If the range of composite scores (max - min) is less than 0.15, flag low differentiation:
- Recommend revisiting rubric specificity (are the criteria granular enough?)
- Consider narrowing the evaluation focus to distinguish tools more clearly
- Document the low differentiation in the deliverable methodology section

## Technology Radar Ring Mapping (Optional)

Map composite scores to Technology Radar rings for stakeholder communication:

| Composite Score | Ring | Recommendation |
|---|---|---|
| > 0.75 | Adopt | Proven, safe default for immediate use |
| 0.50 - 0.75 | Trial | Ready for pilot testing, allocate time this quarter |
| 0.25 - 0.50 | Assess | Worth exploring, monitor development |
| < 0.25 | Hold | Avoid for new projects, insufficient maturity or fit |

## Example Scored Tool

For illustration, a hypothetical evaluation of Composio MCP:

```
Tool: Composio MCP
Category: MCP Servers > Workflow Automation

Integration Feasibility: 0.75
  - Native MCP server, active GitHub repo (500+ stars), 1+ year maintained
  - Community-proven integration with Claude Code

Gap Coverage (per gap, equal weights):
  - Integration Gap: 1.0 (directly connects tools)
  - Content Pipeline Gap: 0.5 (partial automation)
  - Rhythm Gap: 0.25 (tangential)
  - Knowledge Capture Gap: 0.25 (tangential)
  - Documentation/IP Trail Gap: 0.5 (partial via logging)
  - Owned Channel Gap: 0.0 (does not address)
  Gap Coverage Weighted Avg: 0.42

Cost/Sustainability: 0.75
  - Generous free tier, open-source core

Composite: 0.40 * 0.75 + 0.35 * 0.42 + 0.25 * 0.75 = 0.300 + 0.147 + 0.188 = 0.634

Sensitivity: Robust (top 5 in all 4 weight configurations)
Radar Ring: Trial
```

## Gap Parameterization

The 6 POV-expansion gaps serve as the default evaluation framework, but the scoring matrix is parameterized:

- Phase 1 scope refinement can add, remove, or re-weight gaps
- The composite formula works with any number of gaps (1 to N)
- Gap weights default to equal if not explicitly set
- New gaps added during scope refinement receive the same rubric structure (0.0-1.0 scale)

This ensures the skill remains useful even as workflow gaps evolve over time.
