# Phase 2 Agent Prompt Templates

This reference provides domain-specific Task tool prompts for each Phase 2 research agent. These prompts override the researcher skill's default methodology (which is optimized for scientific literature) to ensure appropriate search strategies for AI tool landscape research.

## Why Custom Prompts Are Required

The researcher skill defaults to PubMed-first methodology, which is appropriate for scientific literature reviews but not for AI tool discovery. Agents 1, 2, and 4 are explicitly directed to use WebSearch and WebFetch as primary search tools. Agent 3 retains PubMed access because scientific AI tools are often published in academic venues.

## Agent 1: MCP Servers and Claude Code Integrations

```
You are researching MCP servers and Claude Code integrations for an AI tool landscape assessment.

IMPORTANT: Use WebSearch and WebFetch for this research. Do NOT use PubMed, bioRxiv, or scientific literature databases. Those are not relevant for this domain.

Search targets:
- GitHub: Search for "mcp-server", "awesome-mcp-servers", Claude Code marketplace
- Product Hunt: Search for MCP-related products launched in the past 6 months
- Developer blogs: Search for "MCP server" reviews and comparisons
- Official docs: Anthropic MCP documentation, Claude Code integrations
- npm/PyPI: Search for MCP server packages

For each tool found, document:
- Name, URL, category (see tool-taxonomy.md categories)
- Integration method (native MCP, plugin, API bridge, SDK)
- Maturity indicators (GitHub stars, last commit date, open issues count, documentation quality)
- Relevance to workflow gaps: {gaps_from_scope}
- Cost model (open-source, free tier, paid tiers, enterprise-only)
- Notable limitations or risks

Minimum output: 5 tools with complete profiles.
Target: 8-12 tools.

Output format: One section per tool with all fields filled. Mark unknown fields as "Unknown - requires further investigation" rather than omitting them.
```

## Agent 2: AI Frameworks and Agentic Workflows

```
You are researching AI frameworks and agentic workflow patterns for an AI tool landscape assessment.

IMPORTANT: Use WebSearch and WebFetch for this research. Do NOT default to PubMed-first methodology. AI frameworks are documented in GitHub repos, developer blogs, and framework documentation sites.

Search targets:
- GitHub trending: Agent frameworks, multi-agent orchestration repos
- Framework comparisons: Search for "LangChain vs CrewAI vs AutoGPT 2025-2026"
- Anthropic docs: Claude Code subagent architecture, multi-agent patterns
- Developer blogs: Agentic workflow pattern discussions on dev.to, Medium, personal blogs
- Conference talks: Search for recent AI agent framework presentations

Frameworks to evaluate (starting list, expand through discovery):
LangGraph, CrewAI, AutoGen/Microsoft Agent Framework, DSPy, Semantic Kernel, OpenAI Swarm, Agno

For each framework, document:
- Name, URL, license, primary language
- Category (see tool-taxonomy.md)
- Orchestration pattern (coordinator-worker, hierarchical, sequential, graph-based, etc.)
- MCP support (native, via extension, planned, none)
- Claude Code compatibility (tested, documented, theoretical, none)
- Maturity indicators (GitHub stars, contributors, release frequency, last commit)
- Relevance to workflow gaps: {gaps_from_scope}
- Cost model (open-source, commercial, hybrid)
- Key differentiator (what makes this framework unique)

Minimum output: 5 frameworks with complete profiles.
Target: 7-10 frameworks.

Output format: One section per framework with all fields filled. Mark unknown fields as "Unknown - requires further investigation" rather than omitting them.
```

## Agent 3: AI-Powered Scientific/Computational Biology Tools

```
You are researching AI-powered scientific and computational biology tools for an AI tool landscape assessment.

For this agent, you MAY use PubMed and bioRxiv IN ADDITION to web search, as scientific tools are often published in academic venues first. Use both search strategies.

Search targets:
- PubMed/bioRxiv: "AI tools computational biology 2025-2026", "bioinformatics AI tools", "machine learning protein engineering"
- GitHub: "bioinformatics AI", "computational biology tools", "protein structure prediction", "genomics AI"
- Nature/Science reviews: AI in biosciences review articles
- Web search: "AI tools for startups computational biology", "best bioinformatics AI tools 2026"
- Specialized databases: Papers With Code (biology section), BioTools registry

For each tool, document:
- Name, URL, category (see tool-taxonomy.md)
- Scientific domain (protein structure, genomics, lab automation, literature extraction, drug discovery, etc.)
- Integration potential with Claude Code / MCP architecture
- Maturity indicators (GitHub stars, citations, user community, documentation quality)
- Relevance to workflow gaps: {gaps_from_scope}
- Cost model (open-source, academic license, commercial)
- Data requirements (what inputs does it need, what formats)

Minimum output: 3 tools with complete profiles.
Target: 5-8 tools.

Output format: One section per tool with all fields filled. Mark unknown fields as "Unknown - requires further investigation" rather than omitting them.
```

## Agent 4: Community Patterns and Emerging Trends

```
You are researching community patterns and emerging AI workflow trends for an AI tool landscape assessment.

IMPORTANT: Use WebSearch and WebFetch. Do NOT use PubMed or scientific databases. Community patterns are documented in forums, blogs, and social media.

Search targets:
- Hacker News: AI workflow discussions, Claude Code threads, MCP discussions
- dev.to / Medium: Agentic workflow pattern articles, AI tool reviews
- GitHub: Trending repos in AI/agents category, awesome lists
- Newsletter archives: Ben's Bites, The Batch, AI Engineering Weekly
- Twitter/X: AI tool recommendations from practitioners (search recent threads)
- Reddit: r/MachineLearning, r/artificial, r/ClaudeAI

Focus on:
- Workflow patterns shared by practitioners (not just products)
- Integration approaches that connect communication tools (Slack, Discord) with knowledge tools (Notion, Confluence) and social platforms (LinkedIn, X)
- Emerging patterns for knowledge capture and documentation automation
- Startup-specific AI adoption strategies and case studies
- Low-code/no-code AI integration patterns

For each pattern or tool discovered, document:
- Name or description, source URL
- Category (see tool-taxonomy.md)
- How it connects to workflow gaps: {gaps_from_scope}
- Community traction indicators (stars, shares, discussion volume, upvotes)
- Maturity assessment (Adopt/Trial/Assess/Hold)
- Key insight (why this pattern matters for startup workflows)

Minimum output: 3 tools or patterns with complete profiles.
Target: 5-8 tools or patterns.

Output format: One section per tool/pattern with all fields filled. Mark unknown fields as "Unknown - requires further investigation" rather than omitting them.
```

## Prompt Template Variables

Each prompt template contains a `{gaps_from_scope}` placeholder that the orchestrator replaces with the actual gap list from Phase 1 scope refinement. Example replacement:

```
Relevance to workflow gaps:
1. Integration Gap (weight: 0.25) - tools don't talk to each other
2. Content Pipeline Gap (weight: 0.20) - every post is separate creative act
3. Knowledge Capture Gap (weight: 0.20) - insights decay rapidly
4. Documentation/IP Trail Gap (weight: 0.15) - decisions disappear
5. Rhythm Gap (weight: 0.10) - no weekly cadence
6. Owned Channel Gap (weight: 0.10) - presence depends on uncontrolled platforms
```

## Output Quality Expectations

Each agent output is validated before proceeding to Phase 3:
- Minimum tool count met (see per-agent minimums above)
- All required fields present for each tool (no omitted fields, "Unknown" is acceptable)
- URLs provided are plausible (validation of URL accessibility is not required)
- No duplicate entries within a single agent's output
