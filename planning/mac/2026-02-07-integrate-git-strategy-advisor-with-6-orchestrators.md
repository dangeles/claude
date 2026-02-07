# Integrate git-strategy-advisor with 6 orchestrators

**Date**: 2026-02-07
**Machine**: mac
**Status**: Complete

## Objective

Add optional git-strategy-advisor integration sections to six orchestrator skills (programming-pm, skill-editor, technical-pm, lit-pm, scientific-analysis-architect, research-pipeline), enabling scope-adaptive git recommendations without breaking existing workflows. Also update git-strategy-advisor itself to remove its exclusion of orchestrators that have their own git logic, add a leaf-node declaration, and reframe "Ecosystem Alignment" from adversarial ("vs.") to collaborative ("with") framing.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow (skill-editor orchestration)
- [x] Remove orchestrator exclusion from git-strategy-advisor "When NOT to Use"
- [x] Add leaf-node declaration to prevent circular invocation
- [x] Reframe Ecosystem Alignment from "vs." to "with" + add maintenance note
- [x] Add integration section to programming-pm Phase 6
- [x] Add integration section to skill-editor Phase 4 (before Step 8)
- [x] Add integration section to technical-pm (standalone after "How to Invoke Agents")
- [x] Add integration section to lit-pm (after Stage 8)
- [x] Add integration section to scientific-analysis-architect (after Phase 6)
- [x] Add integration section to research-pipeline (after completion report + handoff subsection)

## Expected Outcome

All orchestrators can optionally invoke git-strategy-advisor for scope-adaptive git recommendations. Orchestrators with existing git logic (programming-pm, skill-editor) have explicit conflict resolution rules ("existing logic takes precedence unconditionally"). Orchestrators without git logic (technical-pm, lit-pm, scientific-analysis-architect, research-pipeline) gain awareness of git workflow decisions.

## Actual Outcome

Successfully implemented all changes across 7 files. 212 lines added, 4 lines removed. All 11 structural verification checks pass. YAML frontmatter validates for all 7 files. Sync to ~/.claude/ completed. No regressions detected in smoke tests.

## Assessment

**Result**: Success

**Improvements**:
- All 6 orchestrators now know about git-strategy-advisor and can invoke it optionally
- Consistent template pattern across all integration sections (header, trigger, invocation, scope, response handling, confidence handling, advisory note)
- Circular invocation risk eliminated by leaf-node declaration
- Conflict resolution protocol explicit for programming-pm and skill-editor
- De-duplication guidance for technical-pm (coordinates sub-workflows)
- Out-of-repo limitation noted for lit-pm, scientific-analysis-architect, research-pipeline

**Issues**:
- None

**Lessons Learned**:
- The archival compliance integration precedent (commit 3d0341a) provided a reliable template for cross-cutting concern integration
- Having both a consistent template and intentional variation (context-appropriate headers) balances discoverability with natural reading

## Related Commits

- 769b80a: feat: Integrate git-strategy-advisor awareness into 6 orchestrator skills

## Next Steps

- User may optionally test functional integration (invoke an orchestrator and verify git-strategy-advisor is consulted)
- Monitor for any issues with the advisory pattern during regular usage
