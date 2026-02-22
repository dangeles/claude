# Add dynamic depth profile system to lit-pm for calibrated output length

**Date**: 2026-02-22
**Machine**: mac
**Status**: Complete

## Objective

Introduce a dynamic depth calibration system into the lit-pm pipeline so output length and detail adapt to request complexity. The existing complexity tier classification (Simple/Medium/Complex/High-Stakes) now drives a depth_profile that propagates to all specialists. Succinctness (DENSE writing) is the default for Simple and Medium tiers. No hardcoded word limits -- all parameters are qualitative directives and ranges.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Change 1: adaptive-orchestration.md -- Add Depth Profile System section with schema, tier mapping, density guidance, user presentation format, mid-pipeline adaptation
- [x] Change 2: lit-pm/SKILL.md -- Add depth_profile generation in Stage 1, replace "3-5 sections" with profile reference, add to workflow_state, propagate in all specialist handoffs (Stages 4, 5, 7, 8)
- [x] Change 3: literature-researcher/SKILL.md -- Add Depth Profile Interpretation block in Mode 2, update paper count references
- [x] Change 4: lit-synthesizer/SKILL.md -- Add depth profile blocks to Mode 2 (Introduction) and Mode 3 (Final Synthesis) with introduction_scope, augmentation_budget, conclusion_scope
- [x] Change 5: quality-gates.md -- Replace fixed ranges with depth-profile-aware ranges for Stages 3, 4, 5, 6a
- [x] Change 6: stage-specifications.md -- Update section count and paper count references, add depth profile controls note
- [x] Change 7: editor/SKILL.md -- Add DENSE compression mandate to "You DO" list, Common Edits table, and Workflow
- [x] Change 8: handoff-schema.md -- Bump version to 1.1, add depth_profile as optional field

## Expected Outcome

Simple/Medium tier lit-pm runs produce shorter, more focused documents (~40-60% of current length). High-Stakes behavior unchanged. User can override depth at Stage 1 checkpoint.

## Actual Outcome

All 8 files modified, validated, synced to ~/.claude/, and tested. YAML frontmatter intact, all key phrases present, no regressions in other skills.

## Assessment

**Result**: Success

**Improvements**:
- depth_profile system enables dynamic output calibration across entire pipeline
- DENSE writing default for Simple/Medium tiers addresses verbosity complaint
- Behavioral density_guidance text prevents inconsistent interpretation
- All changes backward compatible (absent profile = High-Stakes defaults)

**Issues**:
- Pre-existing sync discrepancy ("1 new, 1 deleted" in skills/) unrelated to these changes

**Lessons Learned**:
- Extending existing classification infrastructure (complexity tiers) is cleaner than adding new detection systems
- Behavioral translations (density_guidance text) are important to prevent abstract enums from being interpreted inconsistently

## Related Commits

- [pending]: feat(lit-pm): add dynamic depth profile system for calibrated output length

## Next Steps

- Test with actual lit-pm invocation on a Simple-tier request to verify shorter output
- Monitor first few runs to calibrate depth profile ranges if needed
- Consider adding "Executive Summary" as fifth tier below Simple if user requests ultra-brief mode
