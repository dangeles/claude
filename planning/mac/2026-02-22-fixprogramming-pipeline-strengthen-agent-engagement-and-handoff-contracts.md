# fix(programming-pipeline): strengthen agent engagement and handoff contracts

**Date**: 2026-02-22
**Machine**: mac
**Status**: Complete

## Objective

Fix four problems in the programming pipeline: (1) copilot not always engaged during Phase 5, (2) junior-developer rarely used due to narrow keyword heuristic, (3) systems-architect almost never invoked via proper Task tool delegation, (4) poor handoffs across agents due to missing dispatch templates and schemas.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Change 0: Add copilot to programming-pm provides_to frontmatter
- [x] Change 1: Add Dispatch Templates section with 6 templates and pre-conditions
- [x] Change 2: Strengthen Phase 3 systems-architect invocation with self-check and producer verification
- [x] Change 3: Broaden junior-developer keyword list (fix regex dot bug) and add Step 1b evaluation
- [x] Change 4: Make copilot review mandatory in Phase 5 with dispatch template and emergency override
- [x] Change 5: Add 4 missing handoff schemas to handoff-schema.md v1.3
- [x] Change 6: Add mandatory delegation evaluation to senior-developer with escalation path
- [x] Change 7: Add handoff frontmatter and pipeline integration section to copilot
- [x] Change 8: Add handoff frontmatter to systems-architect

## Expected Outcome

- Systems-architect reliably engaged via Task tool for every project's architecture phase
- Junior-developer considered for delegation on suitable tasks with broader criteria
- Copilot provides mandatory adversarial review before code merges
- Structured dispatch templates eliminate context loss between agents

## Actual Outcome

All changes implemented and synced successfully. YAML frontmatter validates for all 4 modified skill files. 53 skills pass smoke test with no regressions. Handoff-schema.md bumped to v1.3 with 4 new schemas.

## Assessment

**Result**: Success

**Improvements**:
- Dispatch Templates section provides structured context for all specialist invocations
- Phase 3 now has self-check, Task tool mandate, and producer field verification
- Phase 4 has mandatory Step 1b junior-developer evaluation with JSON state tracking
- Phase 5 copilot is mandatory with path verification, verdict validation, and emergency override
- Senior-developer has delegation evaluation step with escalation after 3 failed cycles
- All specialist handoffs now have defined schemas (copilot_dispatch, copilot_review_return, specialist_dispatch, specialist_return)

**Issues**:
- Pre-existing plugins/installed_plugins.json conflict required --yes flag during sync
- New handoff schemas are NOT yet validated by validate-handoff.py (documented as known gap)

**Lessons Learned**:
- Regex dot bug in keyword matching was caught by edge case analysis before implementation
- Emergency override clause for copilot unavailability prevents workflow deadlock
- Two-layer defense (grep heuristic + Step 1b evaluation) is more robust than either alone

## Related Commits

- (pending): fix(programming-pipeline): strengthen agent engagement and handoff contracts

## Next Steps

- Update validate-handoff.py to support 4 new schema types (copilot_dispatch, copilot_review_return, specialist_dispatch, specialist_return)
- Consider extracting dispatch templates to separate reference file if SKILL.md grows past 2000 lines
- Monitor real-world usage to verify systems-architect and copilot are actually invoked
