# Handoff Schema Definitions

This reference defines the YAML handoff schemas used for inter-phase communication within ai-strategist and for inter-skill handoff from upstream workflows.

## Inter-Phase Handoff Schema (v1.0)

Each phase transition produces a handoff YAML file in the `handoffs/` directory.

### Common Fields

All inter-phase handoffs share these fields:

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: {N}
source_phase_name: "{phase_name}"
target_phase: {N+1}
target_phase_name: "{phase_name}"
created_at: "{ISO-8601 timestamp}"
status: "complete"
```

### Phase 0 -> Phase 1 Handoff

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: 0
source_phase_name: "Archival Guidelines Review"
target_phase: 1
target_phase_name: "Scope Refinement"
created_at: "{timestamp}"
status: "complete"
session:
  session_path: "/tmp/ai-strategist-session-{id}/"
  archival_mode: "advisory"  # advisory | enforced
  guidelines_source: "defaults"  # .archive-metadata.yaml | defaults
invocation:
  mode: "quarterly"  # quarterly | deep_dive | event_triggered
  user_prompt: "{original user prompt}"
  handoff_payload: "{path or null}"
```

### Phase 1 -> Phase 2 Handoff

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: 1
source_phase_name: "Scope Refinement"
target_phase: 2
target_phase_name: "Parallel Research"
created_at: "{timestamp}"
status: "complete"
scope:
  mode: "quarterly"
  active_gaps:
    - name: "Integration Gap"
      weight: 0.25
      description: "Tools don't talk to each other"
    - name: "Content Pipeline Gap"
      weight: 0.20
      description: "Every post is separate creative act"
    # ... additional gaps
  tool_categories:
    - "mcp_servers"
    - "ai_frameworks"
    - "scientific_tools"
    - "community_patterns"
  known_tools:
    - "Composio MCP"
    - "LangGraph"
  depth_breadth: "broad"  # broad | focused
  user_approved: true
```

### Phase 2 -> Phase 3 Handoff

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: 2
source_phase_name: "Parallel Research"
target_phase: 3
target_phase_name: "Strategic Assessment"
created_at: "{timestamp}"
status: "complete"
research:
  agents_completed: 4
  agents_total: 4
  degraded_mode: false
  total_tools_found: 22
  tools_passed_to_assessment: 22  # capped at 30
  convergence:
    tools_found_by_multiple_agents: 5
    high_signal_tools:
      - "Composio MCP"
      - "LangGraph"
  agent_outputs:
    - agent: "Agent 1 - MCP Servers"
      output_file: "research/agent-1-mcp-servers.md"
      tools_found: 8
    - agent: "Agent 2 - AI Frameworks"
      output_file: "research/agent-2-ai-frameworks.md"
      tools_found: 7
    - agent: "Agent 3 - Scientific Tools"
      output_file: "research/agent-3-scientific-tools.md"
      tools_found: 4
    - agent: "Agent 4 - Community Patterns"
      output_file: "research/agent-4-community-patterns.md"
      tools_found: 3
```

### Phase 3 -> Phase 4 Handoff

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: 3
source_phase_name: "Strategic Assessment"
target_phase: 4
target_phase_name: "Roadmap Synthesis"
created_at: "{timestamp}"
status: "complete"
assessment:
  tools_scored: 22
  weight_config: "40/35/25"
  top_tools:
    - name: "Composio MCP"
      composite: 0.634
      radar_ring: "Trial"
    - name: "LangGraph"
      composite: 0.612
      radar_ring: "Trial"
  sensitivity:
    weight_robust_tools: ["Composio MCP", "LangGraph"]
    weight_sensitive_tools: ["Tool X"]
    top_recommendation_stable: true
  per_gap_champions:
    - gap: "Integration Gap"
      champion: "Composio MCP"
      score: 1.0
  output_files:
    scored_matrix: "assessment/scored-tool-matrix.md"
    gap_analysis: "assessment/gap-coverage-analysis.md"
    sensitivity: "assessment/sensitivity-analysis.md"
```

### Phase 4 -> Phase 5 Handoff

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: 4
source_phase_name: "Roadmap Synthesis"
target_phase: 5
target_phase_name: "Adversarial Review"
created_at: "{timestamp}"
status: "complete"
roadmap:
  time_horizons_populated: 3  # out of 4
  total_items: 8
  quick_wins: 2
  short_term: 3
  medium_term: 2
  strategic: 1
  output_file: "roadmap/integration-roadmap.md"
```

### Phase 5 -> Phase 6 Handoff

```yaml
schema_version: "1.0"
workflow_id: "{session-id}"
source_phase: 5
source_phase_name: "Adversarial Review"
target_phase: 6
target_phase_name: "Editorial Polish"
created_at: "{timestamp}"
status: "complete"
review:
  critical_challenges: 2
  challenges_addressed: 2
  challenges_unresolved: 0
  methodology_challenged: false
  output_file: "review/adversarial-review.md"
```

## Inter-Skill Handoff: pov-expansion -> ai-strategist

ai-strategist accepts handoffs from pov-expansion using protocol version 2.0.

### Expected Inbound Schema

```yaml
handoff:
  protocol_version: "2.0"
  source_skill: "pov-expansion"
  target_skill: "ai-strategist"
  created_at: "{timestamp}"
context:
  original_prompt: "{user's original request}"
  problem_type: "strategy"
  gap_analysis:
    gaps:
      - name: "Integration Gap"
        severity: "high"
        description: "Tools don't talk to each other"
      - name: "Content Pipeline Gap"
        severity: "medium"
        description: "Every post is separate creative act"
insights:
  workflow_gaps:
    - "Integration between Slack, Notion, LinkedIn, X is manual"
    - "No automated content pipeline"
research_seeds:
  known_tools:
    - "Composio MCP"
    - "n8n"
```

### Ingestion Protocol

1. Validate handoff file exists and is non-empty
2. Attempt YAML parse
3. Check for `context.gap_analysis` -- if present, extract gaps as Phase 1 pre-populated scope
4. Check for `insights.workflow_gaps` -- if present, use as supplementary context
5. Check for `research_seeds.known_tools` -- if present, include in Phase 2 agent prompts
6. On any parse failure: Fall back to interactive scope refinement in Phase 1, passing the raw file content as background context for the requirements-analyst

### Validation Rules

| Field | Required | Fallback on Missing |
|---|---|---|
| handoff.protocol_version | Yes | Reject if missing or unsupported version |
| context.original_prompt | Yes | Use ai-strategist invocation prompt instead |
| context.gap_analysis | No | Interactive scope refinement in Phase 1 |
| insights.workflow_gaps | No | Proceed without supplementary context |
| research_seeds.known_tools | No | Researcher agents discover tools independently |
