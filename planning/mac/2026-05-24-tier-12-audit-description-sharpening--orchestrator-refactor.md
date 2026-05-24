# Tier-1+2 audit: description sharpening + orchestrator refactor

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

Follow-up audit on 2026-05-24 identified two tiers of issues across the 57-skill catalog:

- **Tier 1** — 9 skills with description-level routing problems (real triggering bugs or routing-disambiguation gaps for overlapping skills).
- **Tier 2** — orchestrator SKILL.md files with reference-quality content embedded inline that should live in `references/` for progressive disclosure.

This entry covers both tiers in a single session.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Tier 1: principal-investigator description rewrite (noun-phrase → verb-phrase)
- [x] Tier 1: researcher / synthesizer / fact-checker / consistency-auditor / devils-advocate / scientific-analysis-architect — add "NOT for X (use Y)" exclusions
- [x] Tier 1: strategist / ai-strategist / requirements-analyst three-way disambiguation
- [x] Tier 1: bioinformatician — remove dead pointers to 4 non-existent reference files; replace with working guidance pointing at plotting-advisor + library skills
- [x] Tier 2: programming-pm — extract Dispatch Templates section (230 lines, 7 templates) to `references/dispatch-templates.md`; SKILL.md slimmed to a pointer
- [x] Tier 2: technical-pm — extract Task Breakdown / Tracking Dashboard formats (65 lines) to `references/output-formats.md`
- [x] Tier 2: cfd-bioreactor sanity check — already well-factored with 9 reference files (5441 lines extracted); no action needed
- [x] Validate frontmatter for all touched skills
- [x] `./sync-config.py push --dry-run` and `./sync-config.py push`
- [x] Verify live state at `~/.claude/skills/`

## Expected Outcome

- Triggering improves on the description-rewritten skills (`principal-investigator`, `bioinformatician`).
- Routing collisions between sibling skills resolve at description level rather than relying on body content.
- `programming-pm` SKILL.md shrinks from 2105 → 1905 lines (−200, 9.5%) without losing any orchestrator capability — dispatch templates are now a single `Read` away during a dispatch event.
- `technical-pm` SKILL.md shrinks from 1076 → 1014 lines (−62, 5.8%) — work-plan and dashboard format templates extracted.
- No behavior changes; only progressive disclosure of reference material.

## Actual Outcome

All planned changes shipped exactly as designed. 12 files modified across 2 commits, plus 2 new reference files created. Sync clean post-push; live skills index reflects all updated descriptions. Both extracted reference files (`programming-pm/references/dispatch-templates.md` 234 lines, `technical-pm/references/output-formats.md` ~80 lines) are accessible via `Read` from the live `~/.claude/skills/` paths.

cfd-bioreactor was the sanity-check target; on inspection its 1509-line SKILL.md is genuinely workflow-essential (delegation + 5 phases + protocols), and it already extracts 5441 lines of reference content across 9 files. No refactor needed.

## Assessment

**Result**: Success

**Improvements**:
- 11 sharper descriptions (Tier 1) extend the May 12 tier-3 pattern (`56a667e`) and today's earlier extension (`8d6d4bb`) to the full set of overlapping research/coordination skills.
- 2 orchestrator SKILL.md files slimmed by ~260 lines combined (Tier 2); both retain full functionality via `references/` pointers.
- The bioinformatician dead-pointer fix removes 4 broken cross-references that would have misled any agent reading the References section.
- Trio disambiguation (strategist / ai-strategist / requirements-analyst) is the first time these three have been cross-routed at description level — useful pattern for future trios that may emerge.

**Issues**:
- `sync-config.py plan` still creates entries under `planning/192/` (IP-derived hostname). Manually moved to `planning/mac/` — same workaround as the previous two entries. Worth a real fix in `sync-config.py` (use `hostname -s` with a sane fallback).

**Lessons Learned**:
- Progressive-disclosure refactors of orchestrators don't need to extract everything — the high-leverage extractions are pure templates and format blocks (dispatch templates, output formats), not workflow logic. Going further would create indirection without payoff.
- Healthy orchestrators (cfd-bioreactor) already do this naturally with 9-file `references/` structures; the audit confirms the pattern works at scale.
- Many "should be a reference" sections give themselves away by having nested H2 headings inside a parent H2 (e.g., programming-pm's Dispatch Templates section had 8 nested H2s for individual templates).

## Related Commits

- `ef792a6`: chore(skills): sharpen 10 descriptions + fix bioinformatician dead refs (Tier 1)
- `ffd7c3f`: chore(skills): extract orchestrator templates to references/ (Tier 2)

## Next Steps

- Fix `sync-config.py plan` hostname detection (worked around in three planning entries now; due for a real fix).
- Consider a programming-pm "Team Composition" section cleanup — the inline section currently duplicates content from `references/team-composition.md`; the inline copy is a quick-reference cheat-sheet but could be tightened. Low priority.
- Optional follow-up audit: skills with no `last_updated` field (~30 skills) — set the field consistently for metadata hygiene. Not a behavior issue.
