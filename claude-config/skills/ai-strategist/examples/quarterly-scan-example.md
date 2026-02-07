# Example: Quarterly AI Landscape Scan

This walkthrough demonstrates a complete ai-strategist quarterly scan from invocation to final deliverable.

## User Invocation

```
ai-strategist "Q1 2026 AI landscape scan for workflow optimization"
```

Detected mode: **Quarterly scan** (broad scope, 4 agents, 15+ tools target)

---

## Phase 0: Archival Guidelines Review

```
[Phase 0/6 - Archival Guidelines Review] Starting session setup
```

1. Created session directory: `/tmp/ai-strategist-session-20260115-100000-12345/`
2. Created subdirectories: handoffs/, research/, assessment/, roadmap/, review/, final/
3. Checked for `.archive-metadata.yaml`: Not found
4. Using workflow defaults (enforcement_mode: advisory, guidelines_source: defaults)
5. Initialized workflow-state.yaml

**QG0 PASSED**: Session directory created, archival summary written, state file initialized.

```
[Phase 0/6 - Archival Guidelines Review] COMPLETE
```

---

## Phase 1: Scope Refinement

```
[Phase 1/6 - Scope Refinement] Dispatching requirements-analyst
```

Dispatched requirements-analyst via Task tool with user prompt context.

**Requirements-analyst output** (abbreviated):

> Confirmed scope for Q1 2026 quarterly scan:
>
> Active gaps (all 6, equal weight):
> 1. Integration Gap (weight: 0.167) -- tools don't talk to each other
> 2. Content Pipeline Gap (weight: 0.167) -- every post is separate creative act
> 3. Rhythm Gap (weight: 0.167) -- no weekly cadence
> 4. Knowledge Capture Gap (weight: 0.167) -- insights decay rapidly
> 5. Documentation/IP Trail Gap (weight: 0.167) -- decisions disappear
> 6. Owned Channel Gap (weight: 0.167) -- presence depends on uncontrolled platforms
>
> Tool categories to scan: All 4 (MCP servers, AI frameworks, scientific tools, community patterns)
> Known tools to include: Composio MCP, LangGraph, n8n
> Depth/breadth: Broad (quarterly mode)

User approved scope.

**QG1 PASSED**: Specific scope with measurable criteria, user approved.

```
[Phase 1/6 - Scope Refinement] COMPLETE
Summary: All 6 gaps active with equal weight, scanning all 4 categories
Next: Phase 2 - Parallel Research
Estimated time: 60-90 minutes
```

---

## Phase 2: Parallel AI Landscape Research

```
[Phase 2/6 - Parallel Research] Dispatching 5 agents in parallel
```

Dispatched 4 researcher agents + 1 brainstorming-pm agent via Task tool (parallel).

**Status board during execution**:

| Agent | Category | Status | Tools Found | Duration |
|---|---|---|---|---|
| Agent 1 | MCP Servers | Complete | 9 | 22m |
| Agent 2 | AI Frameworks | Complete | 8 | 25m |
| Agent 3 | Scientific Tools | Complete | 5 | 18m |
| Agent 4 | Community Patterns | Complete | 6 | 20m |
| Agent 5 | Creative Ideas | Complete | 3 | 38m |

**Agent 1 results** (abbreviated -- 9 tools):
- GitHub MCP Server (Adopt, Integration: 1.0)
- Notion MCP Server (Adopt, Integration: 0.75)
- Slack MCP Server (Trial, Integration: 0.75)
- Composio MCP (Trial, Integration: 0.75)
- n8n MCP Bridge (Trial, Integration: 0.5)
- Playwright MCP (Adopt, Integration: 1.0)
- Context7 (Assess, Integration: 0.5)
- PostgreSQL MCP (Adopt, Integration: 0.75)
- Linear MCP (Trial, Integration: 0.5)

**Agent 2 results** (abbreviated -- 8 frameworks):
- LangGraph (Trial, Claude-compatible)
- CrewAI (Trial, MCP via extension)
- AutoGen/Microsoft Agent Framework (Trial, no native MCP)
- DSPy (Assess, programming framework)
- Semantic Kernel (Assess, Microsoft ecosystem)
- Agno (Assess, lightweight agents)
- OpenAI Swarm (Assess, experimental)
- Claude Code Subagents (Adopt, native)

**Agent 3 results** (abbreviated -- 5 tools):
- BioAgents (Assess, AI-assisted biological reasoning)
- AlphaFold MCP Bridge (Assess, theoretical)
- Papers With Code API (Trial, literature+code discovery)
- Benchling API (Trial, lab automation)
- BioPython AI Extensions (Assess, community-driven)

**Agent 4 results** (abbreviated -- 6 patterns):
- Notion-as-knowledge-base pattern (community-proven)
- Slack-to-documentation pipeline (emerging)
- Weekly AI digest automation (emerging)
- LinkedIn content repurposing workflow (community-proven)
- GitHub Actions for content scheduling (proven)
- RSS-to-Slack monitoring pattern (proven)

**Convergence analysis** (orchestrator-owned):
- Tools found by 2+ agents: Composio MCP (Agents 1, 4), LangGraph (Agents 2, 4), n8n (Agents 1, 4), Notion MCP (Agents 1, 4), Slack MCP (Agents 1, 4)
- High-signal tools: 5
- Total unique tools/patterns: 28
- Passed to Phase 3: 28 (under 30 cap)

**QG2 PASSED**: 28 tools across all 4 categories (threshold: 15).

```
[Phase 2/6 - Parallel Research] COMPLETE
Summary: 28 tools/patterns discovered across all categories, 5 high-signal convergence tools
Key findings:
- MCP ecosystem is maturing rapidly (9 servers found, 3 at Adopt level)
- Composio and n8n appear across multiple agent results (high signal)
- Scientific AI tools have lower MCP integration maturity
Next: Phase 3 - Strategic Assessment
Estimated time: 45-60 minutes
```

---

## Phase 3: Strategic Gap Assessment

```
[Phase 3/6 - Strategic Assessment] Dispatching strategist
```

Dispatched strategist via Task tool with:
- Phase 2 research handoff (28 tools)
- Active gaps from Phase 1 (6 gaps, equal weight)
- Scoring matrix reference (references/scoring-matrix.md)
- Explicit instruction: "Evaluate against workflow gaps, NOT bioreactor project goals"

**Scored tool matrix** (top 10 of 28):

| Rank | Tool | Integration (40%) | Gap Cov (35%) | Cost (25%) | Composite |
|---|---|---|---|---|---|
| 1 | Composio MCP | 0.75 | 0.54 | 0.75 | 0.677 |
| 2 | n8n MCP Bridge | 0.50 | 0.58 | 1.0 | 0.653 |
| 3 | Slack MCP Server | 0.75 | 0.46 | 1.0 | 0.711 |
| 4 | Notion MCP Server | 0.75 | 0.50 | 0.75 | 0.663 |
| 5 | GitHub MCP Server | 1.0 | 0.29 | 1.0 | 0.652 |
| 6 | LangGraph | 0.50 | 0.42 | 1.0 | 0.597 |
| 7 | Notion-as-KB pattern | 0.50 | 0.54 | 1.0 | 0.639 |
| 8 | LinkedIn repurpose | 0.25 | 0.50 | 1.0 | 0.525 |
| 9 | CrewAI | 0.50 | 0.38 | 0.75 | 0.520 |
| 10 | Playwright MCP | 1.0 | 0.17 | 1.0 | 0.660 |

**Per-gap champions**:

| Gap | Champion | Gap Score |
|---|---|---|
| Integration Gap | Composio MCP | 1.0 |
| Content Pipeline Gap | n8n MCP Bridge | 0.75 |
| Rhythm Gap | n8n MCP Bridge | 0.75 |
| Knowledge Capture Gap | Notion MCP + KB pattern | 0.75 |
| Documentation/IP Trail | Notion MCP Server | 0.75 |
| Owned Channel Gap | LinkedIn repurpose workflow | 0.75 |

**Sensitivity analysis summary**:
- Weight-robust (top 5 in all configs): Composio MCP, Slack MCP, Notion MCP
- Weight-sensitive: GitHub MCP (drops from #5 to #9 under gap-focused weights -- high integration but low gap coverage)
- Top recommendation stable: Composio MCP remains #1 or #2 in all configurations

**QG3 PASSED**: All 28 tools scored, each gap has champion above 0.5, sensitivity analysis complete.

```
[Phase 3/6 - Strategic Assessment] COMPLETE
Summary: 28 tools scored, Composio MCP leads as weight-robust #1 recommendation
Key findings:
- MCP servers dominate top rankings (integration advantage)
- n8n emerges as strongest gap-coverage tool
- Scientific tools score lower on integration (expected)
Next: Phase 4 - Roadmap Synthesis
Estimated time: 30-45 minutes
```

---

## Phase 4: Integration Roadmap Synthesis (Orchestrator-Owned)

```
[Phase 4/6 - Roadmap Synthesis] Synthesizing integration roadmap
```

**Integration Roadmap**:

### Quick Wins (This Week)
| Tool | Gap(s) | Approach | Effort | Risk |
|---|---|---|---|---|
| Slack MCP Server | Integration, Rhythm | Install from MCP marketplace | Low | Low |
| GitHub MCP Server | Documentation/IP Trail | Install from MCP marketplace | Low | Low |

### Short-Term (This Month)
| Tool | Gap(s) | Approach | Effort | Risk |
|---|---|---|---|---|
| Notion MCP Server | Knowledge Capture, Docs | Install + configure workspace mapping | Medium | Low |
| Composio MCP | Integration (all) | Install + configure service connections | Medium | Medium |
| Notion-as-KB pattern | Knowledge Capture | Restructure Notion workspace for AI access | Medium | Low |

### Medium-Term (This Quarter)
| Tool | Gap(s) | Approach | Effort | Risk |
|---|---|---|---|---|
| n8n MCP Bridge | Content Pipeline, Rhythm | Deploy n8n instance, build automation workflows | High | Medium |
| LinkedIn repurpose | Content Pipeline, Owned Channel | Build content repurposing pipeline with n8n | High | Medium |

### Strategic (Next Quarter)
| Tool | Gap(s) | Approach | Effort | Risk |
|---|---|---|---|---|
| LangGraph | Integration (advanced) | Evaluate for complex multi-agent orchestration | High | High |

**"What changes if we do nothing" per gap**:
- Integration Gap: Manual copy-paste between tools continues, estimated 5-10 hours/week wasted
- Content Pipeline Gap: Each post remains a separate creative act, output volume stays low
- Rhythm Gap: No consistent weekly cadence, audience engagement remains sporadic
- Knowledge Capture Gap: Insights from conversations, experiments, and reading continue to decay
- Documentation/IP Trail Gap: Decisions made in Slack disappear within weeks
- Owned Channel Gap: Platform dependency risk increases as audience grows on uncontrolled platforms

**Technology Radar assignments**:
- Adopt: Slack MCP, GitHub MCP, Playwright MCP
- Trial: Composio MCP, Notion MCP, n8n MCP Bridge, Notion-as-KB pattern
- Assess: LangGraph, CrewAI, LinkedIn repurpose workflow
- Hold: Scientific AI tools (low MCP integration maturity for now)

**QG4 PASSED**: Items in 4 time horizons, each has actionable next step.

```
[Phase 4/6 - Roadmap Synthesis] COMPLETE
Summary: 8 roadmap items across 4 time horizons, 2 quick wins available this week
Next: Phase 5 - Adversarial Review
Estimated time: 20-30 minutes
```

---

## Phase 5: Adversarial Review

```
[Phase 5/6 - Adversarial Review] Dispatching devils-advocate
```

**Devils-advocate challenges** (abbreviated):

1. **Composio MCP integration risk**: "Composio scores 0.75 on integration but has only been community-maintained for 14 months. What happens if the maintainer abandons the project?"
   - **Mitigation**: Open-source (MIT license), forkable. Short-term risk is low. Add to quarterly monitoring list.

2. **n8n self-hosting burden**: "n8n MCP Bridge requires self-hosted n8n instance. This adds infrastructure overhead for a 1-person startup."
   - **Mitigation**: n8n Cloud offers hosted option ($20/month). Adjust cost score if self-hosting is not viable. Flagged in roadmap.

3. **Scientific tools Hold recommendation**: "Putting all scientific AI tools on Hold may be premature. BioAgents shows promise for the bioinformatics workflow."
   - **Mitigation**: Move BioAgents to Assess. Include in next quarterly scan for re-evaluation.

**QG5 PASSED**: All 3 critical challenges addressed with documented mitigations.

```
[Phase 5/6 - Adversarial Review] COMPLETE
Summary: 3 critical challenges raised, all addressed with mitigations
Next: Phase 6 - Editorial Polish
Estimated time: 20-30 minutes
```

---

## Phase 6: Editorial Polish

```
[Phase 6/6 - Editorial Polish] Dispatching editor
```

Editor produced final deliverable at `final/ai-tool-landscape-assessment.md` using the deliverable template.

**QG6 PASSED**: Consistent voice, no substantive errors, executive summary present, no placeholder text.

```
[Phase 6/6 - Editorial Polish] COMPLETE
Summary: Final deliverable produced with executive summary, scored matrix, and integration roadmap
```

---

## Final Deliverable Summary

The completed assessment identified 28 AI tools and patterns across 4 categories, scored them against 6 workflow gaps, and produced an integration roadmap with 2 quick wins available this week.

**Top 3 recommendations**:
1. **Slack MCP Server** (Quick Win) -- Install immediately for Integration and Rhythm gaps
2. **Composio MCP** (Short-term) -- Most robust integration tool across all weight configurations
3. **n8n MCP Bridge** (Medium-term) -- Strongest gap coverage tool for Content Pipeline and Rhythm

**Total workflow duration**: 3 hours 45 minutes

**Session directory**: `/tmp/ai-strategist-session-20260115-100000-12345/`

---

## Key Takeaways for Users

1. **Quick wins exist**: MCP marketplace servers can be installed in minutes and immediately address Integration and Documentation gaps
2. **Composio is the integration hub**: Highest signal tool across all research agents, weight-robust across sensitivity analysis
3. **Scientific AI tools need more time**: Low MCP integration maturity today, but worth monitoring quarterly
4. **n8n is the automation backbone**: Best gap coverage but requires infrastructure investment
5. **Sensitivity analysis matters**: GitHub MCP ranks high on integration but low on gap coverage -- the base composite masks this without sensitivity analysis
