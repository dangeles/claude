# Sharpen description triggers on 5 software skills

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

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

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Rewrite `systems-architect` description to drop bioinformatics-narrowing
- [x] Rewrite `completion-verifier` description as a verb-phrase with trigger words
- [x] Add "NOT for ..." clause to `bioinformatician`
- [x] Add "NOT for ..." clause to `junior-developer`
- [x] Add "NOT for ..." clause to `systematic-troubleshooter`
- [x] Validate frontmatter for all 5 files
- [x] `./sync-config.py push --dry-run` and `./sync-config.py push`
- [x] Verify all 5 skills are still discoverable at `~/.claude/skills/`

## Expected Outcome

The 5 skills behave the same in their bodies; only their description fields change. Auto-triggering improves: `systems-architect` fires on non-omics architecture work, `completion-verifier` fires on "are we done" / "is this complete", and the three with new exclusion clauses route requests to the right specialist when there's ambiguity.

## Actual Outcome

All 5 description rewrites landed verbatim. YAML frontmatter validates for each file. `./sync-config.py push` succeeded; all 5 updated descriptions appear live in the available-skills index. Sync status: clean.

## Assessment

**Result**: Success

**Improvements**:
- `systems-architect` is no longer description-locked to bioinformatics; should now fire on general Python architecture work matching its body's scope.
- `completion-verifier` description is now a verb-phrase with explicit trigger words ("are we done", "is this complete", "verify completion", "ready to ship", "ready to merge", "mark done") so the skill fires when its intended use case is mentioned.
- The three "NOT for" additions follow the May 12 tier-3 pattern (`56a667e`), keeping consistent discriminator language across the skill catalog.

**Issues**:
- `sync-config.py plan` still creates entries under `planning/192/` (hostname returns the IP-derived `192.168.1.9`). Manually moved to `planning/mac/`, same workaround as the plotting-advisor entry. Worth a real fix in `sync-config.py` (use `hostname -s` with a sane fallback).

**Lessons Learned**:
- Description-only edits remain a fast, low-risk discipline — 5 skills, 6 files (including planning entry), one commit, sync is idempotent and clean.
- The May 12 sharpening pass set a strong pattern; the model now reaches for "NOT for X (use Y)" wording automatically when asked to discriminate between sibling skills.

## Related Commits

- `8d6d4bb`: chore(skills): sharpen description triggers on 5 software skills

## Next Steps

- Fix `sync-config.py plan` hostname detection (still creates entries under `planning/192/`; manually moved each time). Tracked as a follow-up in the plotting-advisor planning entry as well.
