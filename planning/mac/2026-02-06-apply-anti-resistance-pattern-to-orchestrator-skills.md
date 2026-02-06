# Apply anti-resistance pattern to orchestrator skills

**Date**: 2026-02-06
**Machine**: mac
**Status**: Complete

## Objective

Apply the anti-resistance pattern (Delegation Mandate, Self-Work Prohibition, State Anchoring, Tool Selection Decision Tree) to three highest-risk orchestrator skills -- scientific-analysis-architect, programming-pm, and lit-pm -- to eliminate delegation resistance (Active Resistance) and workflow state loss (Passive Forgetting).

Based on root cause analysis from session-20260205-232045-39202 (orchestrator-resistance-analysis.md and anti-resistance-pattern.md).

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add Delegation Mandate, Tool Selection, State Anchoring to scientific-analysis-architect
- [x] Add phase transition reminders and state file read instructions to scientific-analysis-architect
- [x] Replace Python asyncio pseudocode in phase-workflows.md with plain-language Task tool instructions
- [x] Add Delegation Mandate, Tool Selection, State Anchoring, anti-rationalization table to programming-pm
- [x] Remove self-work fallback clauses in programming-pm
- [x] Fix Tools section Task tool contradiction in programming-pm
- [x] Replace bash Wave 1/2/3 pseudocode with plain-language (preserving wave structure) in programming-pm
- [x] Add phase transition reminders and state file read instructions to programming-pm
- [x] Add Stage-adapted Delegation Mandate, Tool Selection, State Anchoring to lit-pm
- [x] Replace YAML parallel execution configs with plain-language Task tool instructions in lit-pm
- [x] Modify Owner labels that describe PM doing specialist work in lit-pm
- [x] Add non-linear stage transition reminders (6a/6b/6c, conditional 7.5) to lit-pm
- [x] Add state file read instructions to lit-pm

## Expected Outcome

- All 3 orchestrator skills have explicit delegation mandates within first 40 lines
- All 3 skills have Tool Selection decision trees (Task tool as default)
- All 3 skills have State Anchoring protocol with position prefixes
- Zero pseudocode delegation examples remain
- Zero self-work fallback clauses in programming-pm
- Wave 1/2/3 parallel structure preserved in plain language
- lit-pm uses Stage terminology consistently

## Actual Outcome

All success criteria verified programmatically:
- SC-1 through SC-13: All PASS
- QC-1 through QC-5: All PASS
- All 3 skills load from ~/.claude/ without errors
- No regressions in other skills (skill-editor, completion-verifier, technical-pm)

Line counts:
- scientific-analysis-architect: 463 lines (was 402, under 480 target)
- programming-pm: 1544 lines (was 1514, net +30)
- lit-pm: 614 lines (was 522, over 500 recommendation but acceptable for 9-stage pipeline)

## Assessment

**Result**: Success

**Improvements**:
- All 3 skills now have explicit delegation mandates within first 40 lines
- State anchoring protocol prevents Passive Forgetting
- Tool Selection decision tree eliminates Skill/Task confusion
- Phase/Stage transition reminders at every boundary
- Zero pseudocode (Python asyncio, bash scripts, YAML configs) for tool invocations
- Moderate imperative language (bold, not ALL-CAPS) per Anthropic Claude 4 guidance

**Issues**:
- lit-pm at 614 lines exceeds 500-line recommendation; future session should refactor detailed content into reference files
- programming-pm at 1544 lines is a pre-existing issue not addressed by this session

**Lessons Learned**:
- Test-first approach (scientific-analysis-architect first) was valuable for validating the pattern before batch application
- Per-skill customization (Phase vs Stage, state file format, specialist roster) is essential
- Orchestrator-owned task carve-out prevents over-prohibition of legitimate PM tasks

## Related Commits

- (pending): feat(scientific-analysis-architect): Apply anti-resistance pattern
- (pending): feat(programming-pm): Apply anti-resistance pattern
- (pending): feat(lit-pm): Apply anti-resistance pattern (Stage-adapted)

## Next Steps

- Monitor orchestrator behavior in real workflows to verify effectiveness
- Consider applying pattern to remaining 2 orchestrators (technical-pm, research-pipeline)
- Refactor lit-pm to move detailed content to reference files (reduce from 614 to <500 lines)
- Refactor programming-pm to move detailed content to reference files (reduce from 1544)
