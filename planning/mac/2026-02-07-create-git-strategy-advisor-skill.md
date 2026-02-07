# Create git-strategy-advisor skill

**Date**: 2026-02-07
**Machine**: mac
**Status**: Complete

## Objective

Create a lightweight, advisory-only skill that analyzes work context (planned or completed) and recommends git workflow strategy. The skill produces four decisions -- branch strategy, branch naming, push timing, and PR creation -- as structured YAML output with confidence calibration.

This addresses the gap where git workflow decisions are either hardcoded into orchestrators (programming-pm always branches, skill-editor always direct-commits) or left entirely to the user.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Create `claude-config/skills/git-strategy-advisor/SKILL.md` (main skill, 346 lines)
- [x] Create `claude-config/skills/git-strategy-advisor/references/decision-matrix.md` (decision tables, 181 lines)
- [x] Create `claude-config/skills/git-strategy-advisor/examples/pre-work-invocation.md` (145 lines)
- [x] Create `claude-config/skills/git-strategy-advisor/examples/post-work-invocation.md` (200 lines)

## Expected Outcome

A new leaf-level skill that:
- Supports pre-work and post-work invocation modes
- Classifies work scope (trivial/minor/moderate/major) using deterministic aggregation
- Classifies work type (code/docs/config/skill/mixed) using majority rule
- Produces four strategy decisions with rationale and confidence
- Handles edge cases: non-git repos, detached HEAD, merge/rebase in progress, clean trees
- Is advisory-only (reads git state, never modifies)

## Actual Outcome

All four files created and synced to ~/.claude/skills/ successfully. The skill appears in the available skills list as `git-strategy-advisor`. All validation checks passed:

- YAML frontmatter valid (name: git-strategy-advisor, description: 191 chars)
- All 8 required sections present in SKILL.md
- Line count within target range (346 lines, target 250-350)
- All 44 skills validate after addition (no regressions)
- Source-to-target diff: identical (all 4 files)

## Assessment

**Result**: Success

**Improvements**:
- Future orchestrators can invoke git-strategy-advisor for context-sensitive git recommendations
- Standalone users can ask for git advice directly
- Decision matrix is externalized and independently maintainable

**Issues**:
- None observed during implementation

**Lessons Learned**:
- The sync-config.py push command requires interactive conflict resolution for pre-existing divergences; plan for this when syncing
- SKILL.md line count target of 250-350 was appropriate; the implementation came in at 346 lines with all edge case handling included

## Related Commits

- [pending]: feat(git-strategy-advisor): Create git strategy advisory skill

## Next Steps

- Test in real scenarios (standalone, orchestrator pre-work, orchestrator post-work)
- Consider integration with future orchestrators
- v2 features deferred: commit message generation, PR description generation, traceability metadata
