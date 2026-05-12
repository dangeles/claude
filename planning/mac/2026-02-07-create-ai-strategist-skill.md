# Create ai-strategist Orchestrator Skill

## Date
2026-02-07

## Objective
Create a Tier 1 orchestrator skill that coordinates parallel research agents to scan the AI tool landscape, evaluate tools against workflow gaps using a weighted scoring framework with sensitivity analysis, and produce a prioritized integration roadmap.

## Changes

### Files Created (10 files)

1. `claude-config/skills/ai-strategist/SKILL.md` -- Primary orchestrator skill (405 lines)
   - YAML frontmatter with handoff schema (protocol v2.0)
   - Delegation mandate, anti-resistance table, tool selection
   - 3 invocation modes: quarterly, deep-dive, event-triggered
   - Pre-flight validation for 5 required + 1 optional skill
   - 7-phase pipeline (Phase 0-6): archival, scope, parallel research, assessment, roadmap, review, editorial
   - Consolidated quality gate table with quality floor and override protocol
   - RACI matrix for 7 activities across 7 roles
   - Saga-style error handling with circuit breaker and graceful cancellation
   - Inter-phase status reporting with Phase 2 live status board

2. `claude-config/skills/ai-strategist/references/scoring-matrix.md` -- Weighted scoring framework
   - 3 dimensions: integration feasibility (40%), gap coverage (35%), cost (25%)
   - Detailed 5-level rubrics for each dimension
   - Composite formula and sensitivity analysis protocol (4 weight configurations)
   - Per-gap champion analysis and differentiation check
   - Technology Radar ring mapping (optional)

3. `claude-config/skills/ai-strategist/references/tool-taxonomy.md` -- Two-dimensional AI tool classification
   - Functional categories: MCP Servers, AI Frameworks, Scientific Tools, Workflow Patterns
   - Maturity levels: Adopt, Trial, Assess, Hold
   - Gap-category mapping and seed tools per category

4. `claude-config/skills/ai-strategist/references/agent-prompts.md` -- Domain-specific Phase 2 prompts
   - Agent 1: MCP Servers (WebSearch-first, no PubMed)
   - Agent 2: AI Frameworks (WebSearch-first, no PubMed)
   - Agent 3: Scientific Tools (PubMed + WebSearch)
   - Agent 4: Community Patterns (WebSearch-first, no PubMed)

5. `claude-config/skills/ai-strategist/references/error-handling.md` -- Error handling protocols
   - Saga-style compensation table
   - Circuit breaker pattern
   - Graceful cancellation protocol
   - Atomic state writes
   - State recovery protocol (--resume)
   - Global timeout handling

6. `claude-config/skills/ai-strategist/references/session-structure.md` -- Session directory tree

7. `claude-config/skills/ai-strategist/references/timeout-config.md` -- Per-phase and mode-specific timeouts

8. `claude-config/skills/ai-strategist/references/handoff-schema.md` -- Inter-phase and inter-skill handoff schemas

9. `claude-config/skills/ai-strategist/references/deliverable-template.md` -- Final deliverable structure template

10. `claude-config/skills/ai-strategist/examples/quarterly-scan-example.md` -- Complete walkthrough of a quarterly scan

## Testing

- YAML frontmatter validates (name, description, handoff schema)
- SKILL.md: 405 lines, 2871 words (target: 400-500 lines)
- All 8 reference files present
- Example file present
- No hardcoded user paths found
- Zero occurrences of MUST/NEVER/ALWAYS (threshold: <10)
- Description starts with "Use when..." (no CSO anti-pattern)
- Dry-run sync shows correct 10 files
- Actual sync successful for all 10 ai-strategist files
- Smoke test: 5 skills load correctly (no regressions)

## Outcome
Success. All files created, validated, synced, and committed.

## Analysis Reports Used
- best-practices-review.md: APPROVED with no critical blocking issues
- external-research.md: Sensitivity analysis as #1 recommendation
- edge-cases.md: 22 edge cases identified, 3 critical, 6 important
- knowledge-engineering.md: Structural completeness analysis with targeted gap-filling
