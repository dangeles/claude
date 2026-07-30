# Lighten essay-pipeline: remove state machine and quantitative goals

**Date**: 2026-07-30
**Machine**: mac
**Status**: Complete

## Objective

essay-pipeline produced more bookkeeping than prose. The skill mandated it: a 270-line
session-state schema with atomic `.tmp`/`.bak` writes, 7 quality gates (G1–G5), a state
anchor prefixed to every response, a per-paragraph annotation block (coverage tags, source
manifest, voice notes, word count), and a change-impact assessment with artifact archiving
on any backtrack.

It also defined success numerically. Stage 2 required per-section word-count targets, a
negotiated total, and a proportionality check flagging sections under 150 words. Voice was
scored 1–5 with an action table per score. The completion summary reported word count,
section count, verified-claim count, and override count. Frontmatter `success_criteria`
measured process compliance ("User completes all 4 stages"), never whether the essay was
any good.

Opus 5 follows written procedure more literally than its predecessors, so scaffolding that
earlier models ad-libbed past now gets executed verbatim. The fix is deletion, not more
instructions — `~/repos/essays/CLAUDE.md` already said "value speed and simplicity" and was
losing to the concrete procedure.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Rewrite `essay-pipeline/SKILL.md` (238 → 105): no gates, no state anchors, no session
      protocol; explicit prohibition on word counts and note-heavy tracking
- [x] Delete `references/session-state-schema.md` (270 lines)
- [x] Replace four stage files with four movement files (thesis, shape, evidence, drafting),
      dropping "Circumscribed Responsibility" blocks, checklists, and all length targets
- [x] `references/fact-check-tiers.md` → `fact-checking.md`: tier taxonomy and deferred
      queue removed, verification habit kept
- [x] Rewrite the walkthrough example to demonstrate the light flow
- [x] `essay-voice-matcher`: 1–5 rubric replaced with prose judgment
- [x] `essay-fact-checker`: tier language stripped
- [x] Update all three essay agent definitions to match

## Expected Outcome

Sessions that spend effort on the essay rather than on tracking the essay. Quality judged
semantically — one clear point, sections that earn their place, prose that sounds like the
user.

## Actual Outcome

13 files changed: 265 insertions, 2,201 deletions. essay-pipeline went from ~2,080 lines to
~460. Pushed with `--delete` to remove the 6 orphaned reference files from `~/.claude`
(verified: exactly 6, nothing else). All frontmatter validates; the live skill list reloaded
with the new descriptions, confirming the sync.

## Assessment

**Result**: Success

**Improvements**:
- Deleted the mechanisms rather than instructing against them. A retained mechanism gets
  executed; competing instructions resolve toward the specific one.
- Success criteria are now semantic in all three skills and all three agents, so the
  anti-bureaucracy rule holds whether the pipeline is entered via skill or via agent.
- Kept what was good: the Socratic thesis questions, the devil's-advocate pass, the
  source-quality hierarchy, and degradation-over-failure.

**Issues**:
- `sync-config.py push` prompts interactively on conflicts; `--yes` is required for
  non-interactive runs. `--delete` is also required whenever a change removes files, or
  stale copies linger in `~/.claude` and can still be loaded.

**Lessons Learned**:
- The same over-specification likely exists in the other orchestrators written in the same
  era — `lit-pm`, `research-pipeline`, `program-officer`. Worth a sweep.
- Numeric targets in a skill are self-fulfilling: the model optimizes whatever is counted.
  Prefer prose judgment over any score, including 1–5 quality ratings.

## Related Commits

- [commit SHA]: [commit message]

## Next Steps

- Paper-reading fix: new `paper-reader` skill (comprehension over statistics extraction)
  plus a PreToolUse hook warning when full-text PDFs are read into the main context.
- Field-test the lightened pipeline on a real essay before sweeping the other orchestrators.
