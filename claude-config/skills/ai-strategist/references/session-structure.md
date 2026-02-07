# Session Directory Structure

This reference defines the full directory tree created and maintained by ai-strategist during a workflow run.

## Directory Tree

```
/tmp/ai-strategist-session-{YYYYMMDD-HHMMSS}-{PID}/
|-- archival-guidelines-summary.md        # Phase 0: Archival guidelines or defaults
|-- workflow-state.yaml                    # Persistent workflow state for resume
|-- workflow-state.yaml.bak               # Single backup (created on each update)
|-- handoffs/
|   |-- phase0-session-handoff.yaml       # Phase 0 -> Phase 1
|   |-- phase1-scope-handoff.yaml         # Phase 1 -> Phase 2
|   |-- phase2-research-handoff.yaml      # Phase 2 -> Phase 3
|   |-- phase3-assessment-handoff.yaml    # Phase 3 -> Phase 4
|   |-- phase4-roadmap-handoff.yaml       # Phase 4 -> Phase 5
|   |-- phase5-review-handoff.yaml        # Phase 5 -> Phase 6
|-- research/                              # Phase 2 parallel agent outputs
|   |-- agent-1-mcp-servers.md
|   |-- agent-2-ai-frameworks.md
|   |-- agent-3-scientific-tools.md
|   |-- agent-4-community-patterns.md     # Optional (if Agent 4 dispatched)
|   |-- agent-5-brainstorming.md          # Optional (if brainstorming-pm available)
|   |-- convergence-analysis.md           # Orchestrator-produced cross-agent analysis
|-- assessment/                            # Phase 3 strategist outputs
|   |-- scored-tool-matrix.md
|   |-- gap-coverage-analysis.md
|   |-- sensitivity-analysis.md
|-- roadmap/                               # Phase 4 orchestrator synthesis
|   |-- integration-roadmap.md
|-- review/                                # Phase 5 adversarial review
|   |-- adversarial-review.md
|-- final/                                 # Phase 6 polished deliverable
|   |-- ai-tool-landscape-assessment.md   # Final deliverable
|   |-- partial-completion-summary.md     # Only if workflow timed out or was interrupted
```

## Session ID Pattern

Session directories follow the pattern:
```
/tmp/ai-strategist-session-{YYYYMMDD-HHMMSS}-{PID}/
```

Example: `/tmp/ai-strategist-session-20260115-100000-12345/`

The session ID for resume purposes is the timestamp-PID portion: `20260115-100000-12345`

## Cleanup Policy

- **On successful completion**: Session directory is retained for reference. User may delete manually.
- **On failure or abort**: Session directory is preserved for resume capability.
- **Stale sessions** (>72 hours since last update): Resume triggers a staleness warning.

## File Ownership

| File | Created By | Phase |
|---|---|---|
| archival-guidelines-summary.md | Orchestrator | 0 |
| workflow-state.yaml | Orchestrator | 0 (created), all phases (updated) |
| handoffs/*.yaml | Orchestrator | Each phase transition |
| research/agent-*.md | Researcher agents | 2 |
| research/convergence-analysis.md | Orchestrator | 2 (post-research) |
| assessment/*.md | Strategist | 3 |
| roadmap/integration-roadmap.md | Orchestrator | 4 |
| review/adversarial-review.md | Devils-advocate | 5 |
| final/ai-tool-landscape-assessment.md | Editor | 6 |
