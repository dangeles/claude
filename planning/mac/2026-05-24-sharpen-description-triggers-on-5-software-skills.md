# Sharpen description triggers on 5 software skills

**Date**: 2026-05-24
**Machine**: mac
**Status**: In Progress

## Objective

Continuation of the May 12 tier-3 description-sharpening pass (`56a667e`). That pass covered 10 overlapping skills; this one addresses 5 software-side skills that were missed or have material triggering issues identified in a follow-up audit on 2026-05-24.

Two issues are real bugs (skills under-trigger in practice):

1. `systems-architect` — description over-narrows to "bioinformatics pipelines", but the body is general (`metadata.workflow: software-development`). Would not fire on non-omics architecture work.
2. `completion-verifier` — description is a noun-phrase ("Verification system that...") rather than the verb-phrase pattern used everywhere else. Doesn't catch trigger phrases like "are we done", "is this complete", "verify completion".

Three are polish (add "NOT for X (use Y)" exclusions, matching the May 12 sharpening pattern):

3. `bioinformatician` — currently over-triggers on any -omics request even when senior-developer + a library skill would suffice.
4. `junior-developer` — no exclusion against feature-dev's unscoped work or senior-developer's autonomous scope decisions.
5. `systematic-troubleshooter` — broad description; should exclude trivial syntax/build errors and notebook-specific issues (where notebook-debugger applies).

## Changes Planned

- [ ] Follow CONFIG_MANAGEMENT.md workflow
- [ ] Rewrite `systems-architect` description to drop bioinformatics-narrowing
- [ ] Rewrite `completion-verifier` description as a verb-phrase with trigger words
- [ ] Add "NOT for ..." clause to `bioinformatician`
- [ ] Add "NOT for ..." clause to `junior-developer`
- [ ] Add "NOT for ..." clause to `systematic-troubleshooter`
- [ ] Validate frontmatter for all 5 files
- [ ] `./sync-config.py push --dry-run` and `./sync-config.py push`
- [ ] Verify all 5 skills are still discoverable at `~/.claude/skills/`

## Expected Outcome

The 5 skills behave the same in their bodies; only their description fields change. Auto-triggering improves: `systems-architect` fires on non-omics architecture work, `completion-verifier` fires on "are we done" / "is this complete", and the three with new exclusion clauses route requests to the right specialist when there's ambiguity.

## Actual Outcome

[Filled in after implementation.]

## Assessment

[Filled in after implementation.]

## Related Commits

- [pending]

## Next Steps

- Fix `sync-config.py plan` hostname detection (still creates entries under `planning/192/`; manually moved each time). Tracked as a follow-up in the plotting-advisor planning entry as well.
