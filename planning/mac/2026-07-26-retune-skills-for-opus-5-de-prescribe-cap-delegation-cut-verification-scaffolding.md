# Retune skills for Opus 5 (de-prescribe, cap delegation, cut verification scaffolding)

**Date**: 2026-07-26
**Machine**: mac
**Status**: Success

## Objective

The 55 skills in this repo were written for older Claude models. Anthropic's current
guidance states plainly that skills built for prior models are often too prescriptive and
degrade output quality on current ones, and documents a 500-line ceiling for a SKILL.md
body. An audit against the named anti-patterns found the repo out of step in five ways:

- 162 files carried `Step N` / `Phase N` / `Stage N` choreography; 19 of 55 SKILL.md files
  exceeded 500 lines, topping out at 1,376.
- 91 files used `CRITICAL` / `MUST` / `MANDATORY` as a default register rather than as a
  targeted escalation.
- 6 orchestrators carried a verbatim-duplicated "You MUST delegate all specialist work"
  mandate, against only ~10 delegation caps.
- 20 files mandated wall-clock progress cadences ("every 15 minutes", "set timer or check
  clock") that Claude cannot execute, since it has no timer and cannot wake itself.
- 3 skills instructed use of "extended thinking" with explicit token budgets, a mechanism
  removed from the API (`budget_tokens` now returns a 400).

The model target drove two of these in the opposite direction from the initial plan. The
user runs **Opus 5** and cannot use Fable 5 at all — its safety classifiers target biology
and life-sciences content, including lab methods and molecular mechanisms, which is this
user's core work. Under Opus 5 guidance delegation gets *capped* and verification
scaffolding gets *deleted*; under Fable 5 guidance both are encouraged. Confirming the
model first prevented retuning 55 skills in the wrong direction.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Correct the model pin (it named a model the user cannot use)
- [x] Delete obsolete mechanics: extended-thinking budgets, wall-clock cadences
- [x] De-prescribe open-field skills; keep narrow-bridge sequences explicit
- [x] Fix the forced-finding quota; demote aggressive language to plain prose
- [x] Add only the Opus 5 guidance the harness does not already inject
- [x] Re-validate, dry-run, push, test in a fresh session, commit

## Expected Outcome

Every SKILL.md under 500 lines, no instructions referencing removed API mechanisms, and
delegation and verification tuned to how Opus 5 actually behaves.

## Actual Outcome

77 files changed, 8 added. SKILL.md content fell from 22,685 to 16,725 lines (−26%);
files over the 500-line limit went from 19 to 0. Largest reductions: `technical-pm`
1014→337, `cfd-bioreactor` 1376→<500, `programming-pm` 850→497, `lit-synthesizer`
590→208, `systematic-troubleshooter` 519→172, `data-pipeline-manager` 694→254,
`completion-verifier` 312→78.

Work was split: cross-cutting and high-stakes edits done directly, and 22 large skills
retuned by 8 parallel subagents against a single shared written brief so the rules were
applied consistently. Every agent respected the hard constraint on frontmatter — verified
programmatically that only one `description` changed repo-wide, the one deliberately
edited to stop advertising "extended thinking".

**Applied the degrees-of-freedom principle rather than blanket de-prescription.** Anthropic's
skill-authoring guidance still endorses explicit workflows for fragile tasks, so narrow-bridge
sequences were kept verbatim: the sync/push/validate/commit sequence, handoff-schema field
names and `validate-handoff.py` invocations, LaTeX multi-pass compilation ordering, FEniCSx
solver setup, kernel-restart procedures, `yaml.safe_load` over `yaml.load`, and every
destructive-action guard (`NEVER git add .`, force-push warnings, per-file approval gates).

**Item 5 shrank on inspection.** Claude Code's own system prompt already injects most of the
recommended Opus 5 blocks — communication style, scope discipline, autonomous operation,
"report outcomes faithfully", and the check-your-last-paragraph clause. Adding them again
would have been duplication the conciseness guidance warns against. Only the two things the
harness does not provide were added, in one shared
`claude-config/skills/references/delegation-and-scope.md` linked from the orchestrators
rather than pasted into each.

### Pre-existing bugs found along the way

- **`notebook-writer/SKILL.md` was committed as lowercase `skill.md` in the very first
  commit (`5ca1b89`).** macOS sets `core.ignorecase=true`, so it resolved locally and looked
  correct, while `sync-config.py` globs `*/SKILL.md` case-sensitively — meaning on Linux,
  including CI, that skill was invisible to both the push gate and skill discovery. CI had
  been validating 54 of 55 skills. Fixed in the index; `tests/test_skill_layout.py` added
  (3 tests) and verified to catch it against `HEAD`.
- **Stale `Co-Authored-By` trailers** pinning Claude Opus 4.6 and Sonnet 4.5 in 5 commit
  templates. Unpinned to `Claude <noreply@anthropic.com>`, since a hardcoded model name in a
  template is guaranteed to go stale.
- **Stale model-selection guidance** in `skill-editor` naming Sonnet 4.5 as default and Opus
  4.6 as premium. Rewritten against the current lineup, with the Fable 5 life-sciences
  caveat recorded so it isn't rediscovered the hard way.
- **35 broken internal links** and several stale `skill-editor` references describing agents
  and a Phase 2.5 removed per `skills/CHANGELOG.md:66`. The stale references were rewritten;
  the broken links are logged as follow-up.
- `CLAUDE.md` claimed 56 skills; there are 55.

## Assessment

**Result**: Success

**Improvements**:
- Session default is a model the user can actually use for biology work.
- Every SKILL.md is within the documented limit; ~6,000 fewer lines load on skill trigger.
- No skill instructs a removed API mechanism or an unexecutable timed cadence.
- Delegation and verification now match Opus 5's actual behaviour instead of working
  against it.
- A whole class of filesystem-case bugs is now caught by CI rather than by luck.

**Issues**:
- The forced-cadence instructions were worse than dead weight. Combined with mandated
  periodic status updates and no clock, they invited exactly the fabricated-progress
  failure mode the guidance warns about.
- Several `skill-editor` references documented a workflow that no longer existed. Reference
  files drift faster than the SKILL.md that points at them, and nothing tests them.
- The push-time gate gave false confidence: it validated 54 of 55 skills for the repo's
  entire history without reporting a discrepancy.

**Lessons Learned**:
- Confirm the model before retuning prompts for it. Opus 5 and Fable 5 guidance is opposite
  on delegation and verification; guessing would have made 55 skills worse in both.
- Check what the harness already provides before adding recommended blocks. Roughly
  two-thirds of the "add this instruction" advice was already in the system prompt.
- `core.ignorecase=true` on macOS hides an entire class of bug that only appears in CI.
  Assertions about repo layout belong against `git ls-files`, not the working tree.
- Reflexively escalating to `MUST` spends a lever that only works while it is rare.

## Related Commits

- [commit SHA]: [commit message]

## Next Steps

Housekeeping backlog, in rough priority order:

1. Decide whether `completion-verifier` should be fully deprecated (repo precedent:
   `parallel-coordinator`, `software-developer`) rather than retuned. Its 374-line
   `examples/verification-report-example.md` is now orphaned either way.
2. Repair or delete the 35 broken internal links; several skills promise `references/`
   files that were never written (`copilot`, `systems-architect`,
   `biologist-commentator`, `principal-investigator`). Consider a link-check test.
3. Delete orphaned `skill-editor` artifacts from the removed Phase 2.5:
   `phase-2-5-detection.sh` (476 lines), `mode-detection.sh`, and four report templates.
4. Remove the two deprecated stub skills (`parallel-coordinator`, `software-developer`) or
   confirm they stay as tombstones.
5. Fix `bioinformatician/scripts/`, which documents `qc_pipeline.py` and
   `differential_expression_template.py` — neither exists on disk.
6. Residual wall-clock content outside the retune's scope:
   `technical-pm` frontmatter `estimated_duration`, and the Timeout Configuration section
   of `technical-pm/references/parallel-execution.md`.
7. `project-configs/` is empty while 24 repos carry `.claude/` config; decide whether
   per-project config is tracked or the directory goes.
8. Consider whether `superpowers` 6.2.0 should stay enabled. Its `using-superpowers` skill
   is the most aggressively-worded prompt in every session ("YOU DO NOT HAVE A CHOICE",
   a 12-row rationalization table) and is exactly the pattern this pass removed everywhere
   else — but it lives in `plugins/cache/` and cannot be edited from this repo.
