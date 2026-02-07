# Add explicit library-pm Task tool dispatch templates to archive-workflow

**Date**: 2026-02-07
**Machine**: mac
**Status**: Complete

## Objective

Apply the anti-resistance pattern to archive-workflow, the 4th orchestrator in the ecosystem. This adds explicit Task tool dispatch templates for the library-pm orchestrator role, matching the pattern established in scientific-analysis-architect (27d90c0), programming-pm (3778b4a), and lit-pm (fb53b0b).

## Changes Planned

- [x] Add Delegation Mandate section with library-pm role boundaries (after Architecture Diagram, before Pre-flight Checks)
- [x] Add anti-resistance rationalization table with archive-workflow-specific examples
- [x] Add Tool Selection decision tree (Task tool as default)
- [x] Add State Anchoring protocol ([Wave N/4] format)
- [x] Replace brief single-agent dispatch (9 lines) with full 4-wave templates (5 agents)
- [x] Include archival_context: "skip" in all dispatch prompts (circular dependency prevention)
- [x] Update code-project-organization.md example with dispatch traces
- [x] Update research-project-organization.md example with dispatch traces
- [x] Update mixed-project-organization.md example with dispatch traces

## Expected Outcome

- archive-workflow SKILL.md contains Delegation Mandate within first 100 lines
- All 5 dispatch templates include archival_context: "skip"
- Wave 2 template explicitly shows simultaneous launch
- SKILL.md line count under 650
- No aggressive language excess (< 15 MUST/NEVER/ALWAYS)
- All 3 examples show [Wave N/4] state anchoring at wave boundaries
- Existing content (Rollback Procedure, Conflict Resolution, Quality Gates, Timeout Config) unchanged
- archival-compliance-check.md reference document unchanged

## Actual Outcome

All success criteria met:
- SKILL.md: 614 lines (under 650 limit)
- Delegation Mandate at line 79
- Tool Selection at line 107
- State Anchoring at line 117
- Anti-resistance table at line 97
- Agent Dispatch Templates at line 550
- 8 occurrences of archival_context: "skip" (5 dispatch + 1 original + 1 header + 1 note)
- 7 "via Task tool" references
- 6 MUST/NEVER/ALWAYS occurrences (well under 15)
- All existing sections intact (verified by line numbers)
- archival-compliance-check.md unchanged (verified with git diff)
- No regressions in other skills (22+ skills reference archive-workflow, all via archival-compliance-check.md)

## Files Modified

- claude-config/skills/archive-workflow/SKILL.md
  - Added Delegation Mandate, Tool Selection, State Anchoring sections (after line 77)
  - Replaced Agent Dispatch section with expanded 4-wave templates
- claude-config/skills/archive-workflow/examples/code-project-organization.md
  - Added dispatch traces with [Wave N/4] state anchoring at each wave
- claude-config/skills/archive-workflow/examples/research-project-organization.md
  - Added dispatch traces with [Wave N/4] state anchoring at each wave
- claude-config/skills/archive-workflow/examples/mixed-project-organization.md
  - Added dispatch traces with [Wave N/4] state anchoring at each wave

## Related

- Planning entry: 2026-02-06-apply-anti-resistance-pattern-to-orchestrator-skills.md
- Commits: 27d90c0, 3778b4a, fb53b0b (prior anti-resistance applications)
- Next Steps section from 2026-02-06 entry: "Consider applying pattern to remaining 2 orchestrators"

## Follow-Up Items

- Apply anti-resistance pattern to remaining orchestrator (technical-pm or research-pipeline)
- Consider moving .archive-metadata.yaml schema (77 lines) to a reference file if SKILL.md grows further
- Consolidated quality gate override protocol (from knowledge-engineering analysis) -- lower priority since archive-workflow already has QG1-QG5 table
