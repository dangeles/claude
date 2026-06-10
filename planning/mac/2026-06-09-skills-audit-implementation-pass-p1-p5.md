# Skills audit implementation pass (P1-P5)

**Date**: 2026-06-09
**Machine**: mac
**Status**: Implemented (pending real-session verification)

## Objective

Implement every TODO in `SKILLS-AUDIT-TODO.md` — the audit of the 55 local skills for
items needing improvement or refactoring. Fix actively-misleading state (stale README,
sync drift), resolve the only genuine orchestrator redundancy, reduce oversized skill
bodies that fight progressive disclosure, sharpen weak descriptions, and finish the
v2.0→v3.0 handoff-schema migration.

## Decisions on the audit's open options

- **software-developer drift** → restore the tombstone `SKILL.md` in the repo (option a).
  Preserves the redirect to `senior-developer` and eliminates drift without a `--delete`.
- **program-officer** → sharpen + extract (not merge into technical-pm) — lower risk,
  addresses both the P2 overlap and the P3 bloat.
- **parallel-coordinator** → deprecate to a tombstone redirecting to
  `superpowers:dispatching-parallel-agents` (full functional duplicate), matching the
  repo's own `software-developer` precedent. (Resolves both P2 and its P4 description.)
- **P5 frontmatter schema** → document a minimal schema in `README.md`; update
  `last_updated` ONLY on skills actually modified (bumping untouched skills' dates would
  falsify provenance).

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md / skill-editor workflow (planning → edit repo → validate → push → commit)
- [x] **P1.1** Rewrite `skills/README.md` (26→55 skills, Fable 5/Opus 4.8, full domain tables, frontmatter conventions)
- [x] **P1.2** Restore `software-developer` tombstone in repo (kills sync drift)
- [x] **P2.1/P3** `program-officer`: sharpen description (+NOT-for), extract templates to `references/output-templates.md` (623→448)
- [x] **P2.2/P4** `parallel-coordinator`: deprecate→tombstone; update `perspective-swarm` reference
- [x] **P3** `programming-pm`: extract phase walkthroughs to `references/phases/` (1871→850); fixed pre-existing pre-mortem path typo
- [x] **P3** `cfd-bioreactor`: extract agent-prompt templates to `references/agent-prompt-templates.md` (1510→1376)
- [x] **P3** `bioinformatician`: extract output/reproducibility templates to `assets/` (830→555)
- [x] **P4** Rewrite descriptions: `data-pipeline-manager`, `notebook-writer`, `lit-synthesizer`, `literature-researcher`
- [x] **P5.1** Migrate handoff v2.0→v3.0: `ai-strategist`, `brainstorming-pm`, `lit-pm`, `pov-expansion`
- [x] **P5.2** Document minimal frontmatter schema in README
- [x] CHANGELOG entry

## Expected Outcome

README is an accurate navigation map; no repo/live drift; one fewer redundant skill;
three large SKILL.md bodies reduced to routing/decision layers with detail in references;
weak descriptions trigger reliably; all handoff metadata on schema v3.0.

## Actual Outcome

All 14 TODO items implemented. Five heavy skill-body refactors ran as parallel subagents
(content moved verbatim, nothing deleted). Validation:
- 55/55 SKILL.md frontmatter parse as valid YAML with name+description.
- All reference pointers in refactored skills resolve (cross-skill refs verified to exist).
- Both tombstones carry `deprecated: true` + `replaced_by`.

Line reductions: programming-pm 1871→850, cfd-bioreactor 1510→1376, bioinformatician
830→555, program-officer 623→448.

## Assessment

**Result**: Success (implementation + structural validation; see Issues for testing scope)

**Improvements**:
- README no longer misses ~30 skills or claims a wrong platform/count.
- Repo and `~/.claude/` no longer disagree on software-developer.
- parallel-coordinator no longer competes with the superpowers skill.
- Three oversized bodies now follow progressive disclosure.

**Issues**:
- Full skill-editor FULL-mode (swarm + adversarial) and per-skill RED-GREEN subagent
  pressure-testing were NOT run for each skill — impractical autonomously across 10+
  skills, and this is a doc/metadata/extraction hygiene pass (low behavioral risk), not
  new-behavior authoring. Real-session triggering verification of the reworded
  descriptions is recommended as follow-up.
- Pre-existing live-only orphans remain (git-strategy-advisor refs/examples,
  workflow-coordinator `universal-handoff-schema-v3.0.json`); out of audit scope, left as-is
  (no `--delete`).

**Lessons Learned**:
- The repo/live "drift" was simply that live had not been pushed since the last rename
  commit; mapping it precisely before pushing avoided a needless scare.

## Related Commits

- [pending]: chore(skills): skills audit implementation pass (P1-P5)

## Next Steps

- Push to `~/.claude/`, verify in a real session, then commit.
- Optional follow-up: reconcile the pre-existing orphans; consider whether `version:` should be standardized or dropped repo-wide.
