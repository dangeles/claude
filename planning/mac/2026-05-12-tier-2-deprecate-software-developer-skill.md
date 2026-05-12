# tier-2-deprecate-software-developer-skill

**Date**: 2026-05-12
**Machine**: mac
**Status**: Success

## Objective

Deprecate `software-developer` skill as a confirmed duplicate of `senior-developer`. Closer reading shows the only unique signal in software-developer is bioinformatics framing (CPM examples, biologist-commentator integration); the actual Python practices (type hints, docstrings, pytest, error handling, CLI) are identical. For users who need bioinformatics-specific implementation, the chain `bioinformatician` → `senior-developer` → `biologist-commentator` provides this without a duplicate skill.

Tier 2 originally had B1+B2+B3 candidates. After closer reading:
- B2 (research-pipeline): kept — 5-stage quick chain is distinct from lit-pm's 9-stage adaptive review. Trigger sharpening moved to Tier 3 C1.
- B3 (code-review vs pr-review-toolkit): kept both — `code-review` is post-PR GitHub automation, `pr-review-toolkit` is pre-PR local interactive. Distinct lanes.
- B1 (software-developer): proceed.

## Changes Planned

- [ ] Tombstone `claude-config/skills/software-developer/SKILL.md` with deprecation notice + redirect to senior-developer (preserve file for git history)
- [ ] Replace `software-developer` → `senior-developer` in `integrates-with:` metadata arrays: copilot, biologist-commentator, systematic-troubleshooter, systems-architect, principal-investigator
- [ ] Replace `software-developer` → `senior-developer` in body prose: programming-pm, principal-investigator (~10 occurrences including 3 `Skill()` invocations), perspective-swarm/references/handoff-schema.md, principal-investigator/references/research-coordination-integration.md, principal-investigator/examples/usage-examples.md
- [ ] Update `claude-config/skills/CHANGELOG.md` with deprecation entry
- [ ] Check `claude-config/skills/README.md` for index references
- [ ] `./sync-config.py push --dry-run` then `push`
- [ ] Verify live `~/.claude/skills/software-developer/SKILL.md` shows tombstone

## Expected Outcome

- One fewer overlapping skill in dispatch space; clearer trigger when user requests "implement production Python".
- All orchestrators reroute to senior-developer transparently.
- software-developer/ directory preserved for git history; tombstone description prevents auto-invocation but the file is discoverable for anyone searching.

## Actual Outcome

All changes applied. Push synced 12 files; tombstone live; references updated. The deprecated skill description now reads 'Deprecated as of 2026-05-12' in the auto-loaded skill list.

## Assessment

**Result**: Success

## Related Commits

- 574736544224a7068be9031cc13c8baeb9d88713: chore(tier-2): deprecate software-developer skill (duplicate of senior-developer)

## Next Steps

Proceed to Tier 3 trigger sharpening (C1-C4).
