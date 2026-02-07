# AI Tool Taxonomy

This reference provides a two-dimensional classification scheme for organizing AI tools evaluated by ai-strategist.

## Classification Axes

### Axis 1: Functional Category

#### MCP Servers

Tools that extend Claude Code capabilities through the Model Context Protocol.

| Subcategory | Description | Example Tools |
|---|---|---|
| Development/VC | Source control, CI/CD, code review | GitHub MCP, GitLab MCP |
| Project Management | Task tracking, sprint management | Notion MCP, Linear MCP |
| Communication | Messaging, notifications | Slack MCP, Discord MCP |
| Web Automation | Browser control, scraping | Playwright MCP, Puppeteer MCP |
| Database | Data access, queries | PostgreSQL MCP, SQLite MCP |
| API Integration | Multi-service connectors | Composio MCP, Zapier MCP |
| Workflow Automation | Pipeline orchestration | n8n MCP, Make MCP |
| Context/Memory | Knowledge persistence | Context7, Mem0 |
| Documentation | Doc generation, management | Notion MCP, Confluence MCP |

#### AI Frameworks

Orchestration and agent frameworks for building multi-agent systems.

| Subcategory | Description | Example Frameworks |
|---|---|---|
| Orchestration Frameworks | Multi-agent coordination | LangGraph, CrewAI, AutoGen |
| Programming Frameworks | LLM application building | DSPy, Semantic Kernel |
| Integration Platforms | Connect LLMs to services | Agno, OpenAI Swarm |

#### Scientific/Domain Tools

AI tools specific to scientific and computational biology workflows.

| Subcategory | Description | Examples |
|---|---|---|
| Bioinformatics AI | Computational biology analysis | BioAgents, BioPython AI extensions |
| Protein Structure | Structure prediction and analysis | AlphaFold tools, ESMFold |
| Genomics | Genome analysis and annotation | Genomics AI pipelines |
| Lab Automation | Experimental workflow automation | Opentrons AI, Benchling integrations |

#### Workflow Patterns

Architectural patterns for connecting tools and agents.

| Subcategory | Description | Examples |
|---|---|---|
| Sequential Pipeline | Linear stage-by-stage processing | lit-pm pattern, programming-pm pattern |
| Parallel Fan-Out | Multiple agents working simultaneously | ai-strategist Phase 2 pattern |
| Coordinator-Worker | Hub delegates to specialist workers | Claude Code subagent architecture |
| Hierarchical Teams | Nested orchestrators with team leads | Multi-level orchestration |

### Axis 2: Maturity Level

Based on Technology Radar ring methodology:

| Level | Criteria | Indicators |
|---|---|---|
| Adopt | Production-grade, official support, active community | >1000 GitHub stars, <7 days since last commit, official documentation, production case studies |
| Trial | Mature, community-proven, ready for pilot | >100 GitHub stars, <30 days since last commit, documented integration paths |
| Assess | Emerging, documented, worth exploring | >10 GitHub stars, <90 days since last commit, basic documentation |
| Hold | Experimental, limited community | <10 stars or >90 days inactive, minimal documentation |

## Classification Criteria for New Tools

When a researcher agent discovers a new tool, classify it using:

1. **Functional category**: Based on primary use case (what problem does it solve?)
2. **Subcategory**: Based on specific domain (which aspect of that problem?)
3. **Maturity level**: Based on GitHub activity, documentation quality, community size, time since first stable release

If a tool spans multiple categories (e.g., Composio serves both API Integration and Workflow Automation), list it under its primary category with a cross-reference note.

## Gap-Category Mapping

This table maps which tool categories tend to address which workflow gaps. Use this to guide researcher agent focus areas during Phase 2.

| Workflow Gap | Primary Categories | Secondary Categories |
|---|---|---|
| Integration Gap | MCP Servers (API Integration, Workflow Automation) | AI Frameworks (Integration Platforms) |
| Content Pipeline Gap | MCP Servers (Communication, Documentation) | Workflow Patterns (Sequential Pipeline) |
| Rhythm Gap | MCP Servers (Workflow Automation, Project Management) | Workflow Patterns (Coordinator-Worker) |
| Knowledge Capture Gap | MCP Servers (Context/Memory, Documentation) | AI Frameworks (Orchestration) |
| Documentation/IP Trail Gap | MCP Servers (Documentation, Database) | Workflow Patterns (Sequential Pipeline) |
| Owned Channel Gap | MCP Servers (Web Automation, Communication) | AI Frameworks (Integration Platforms) |

## Seed Tools per Category

These serve as starting points for researcher agents, not exhaustive lists:

**MCP Servers**: GitHub, Notion, Slack, Composio, n8n, PostgreSQL, Playwright, Context7, Linear, Discord

**AI Frameworks**: LangGraph, CrewAI, AutoGen/Microsoft Agent Framework, DSPy, Semantic Kernel, OpenAI Swarm, Agno

**Scientific Tools**: BioAgents, AlphaFold tools, genomics AI pipelines, Opentrons, Benchling

**Workflow Patterns**: Claude Code subagent architecture, n8n automation, Zapier MCP, Make integrations

Researcher agents should use these as starting points and expand through their own discovery process.
