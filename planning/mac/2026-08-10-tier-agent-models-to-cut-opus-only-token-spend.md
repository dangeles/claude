# Tier agent models to cut Opus-only token spend

**Date**: 2026-08-10
**Machine**: mac
**Status**: Success

## Objective

Token use was far too high. Investigation of `~/.claude/stats-cache.json` (71 sessions,
through 2026-07-29) found that **99.50% of lifetime output tokens and 99.85% of cache reads
ran on Opus**; Sonnet and Haiku had handled 0.5% of output between them. The tiering
machinery existed and was entirely unused.

Root cause was twofold:

1. 16 of 27 agents hardcoded `model: opus`; the other 11 had no `model:` key at all, so they
   inherited — and since `settings.json` pins `"model": "claude-opus-5"`, inheriting also
   meant Opus. Every agent ran Opus.
2. House guidance in `skills/skill-editor/references/anthropic-guidelines-summary.md` said
   verbatim **"All agents use Opus"**, with `skill-editor-executor` as the worked example.
   That doctrine is what produced the pins, and would have re-produced them on the next
   `skill-editor` run.

A related class of staleness: several skills carried "Current | Target" model tables built on
the premise that "the Task tool does not support explicit model selection" (as of ~2026-02).
That premise is false — the Agent/Task tool takes a `model` parameter and agent frontmatter
takes `model:`. Those skills' Sonnet/Haiku intentions had therefore never taken effect.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Tier all 27 agents deliberately: 7 `opus`, 16 `sonnet`, 4 `haiku`
- [x] Rewrite the "All agents use Opus" guidance so the doctrine does not reassert itself
- [x] Correct the false "Task tool has no model parameter" premise in 6 skill files
- [x] Replace pinned version IDs (Sonnet 4.5/4.6, Haiku 4.5, Opus 4.5/4.6) with tier aliases
- [x] Investigate the context-volume lever without changing `settings.json`

## Expected Outcome

Delegated work runs at a tier matched to the task, cutting the rate on subagent spend ~5x
where it moves. Guidance no longer instructs the next editor to re-pin Opus.

## Actual Outcome

**Tiering applied and verified live** (`~/.claude/agents`: `{haiku: 4, sonnet: 16, opus: 7}`).

- `opus` (7): essay-pipeline-orchestrator, latex-writing-expert, python-developer,
  skill-editor-adversarial-reviewer, library-pm, pov-expansion-pm, pov-synthesizer
- `sonnet` (16): coherence-manager, essay-fact-checker, essay-voice-matcher,
  latex-proofreader, seo-manager, website-designer, suggestion-engine, portfolio-manager,
  skill-editor-{request-refiner,edge-case-simulator,executor},
  pov-{perspective-analyst,abstractor-classifier,transfer-evaluator},
  archive-{decision-integrator,expandability-reviewer}
- `haiku` (4): latex-content-examiner, archive-{clutter-analyst,nomenclature-enforcer,structure-organizer}

23 agent files changed (12 `opus`->lower, 11 `model:` inserted where absent). 8 skill/doc
files corrected. `validate` passed (84 files), 80 tests passed / 1 skipped, `status` reports
no divergence between repo and `~/.claude/`.

Also pulled `plugins/installed_plugins.json` (live was newer; diff was 16 lines, all
`lastUpdated` timestamps) to clear the pre-existing drift warning before pushing.

### Context-volume investigation (no settings.json change)

Measured, since the model tier is not the dominant term:

| Fixed per-session floor | ~tokens |
|---|---|
| own skill descriptions (57) | 10,336 |
| plugin agent descriptions (34) | 4,646 |
| plugin skill descriptions (31) | 3,460 |
| own agent descriptions (27) | 1,595 |
| CLAUDE.md | 836 |
| **total before any input** | **~20,900** |

Per-session averages over 71 sessions: **890 messages**, 571K output tokens, 6.7M cache
creation, **97.5M cache read**. That works out to **~110K tokens of context re-read on every
one of ~890 messages per session**. Cache read = context size x request count, and it is 171x
the output token volume.

On-invoke skill cost is modest (median SKILL.md ~3.2K tokens, worst `programming-pm` at
7.8K); the 551K tokens sitting in `references/` load on demand, not at invoke.

## Assessment

**Result**: Success (for the stated scope)

**Improvements**:
- Delegated work now runs at a deliberate tier; fan-out agents (`pov-perspective-analyst`,
  the 5 `archive-*` analysts, 4 `skill-editor-*`) were the ones multiplying Opus cost.
- The doctrine that caused the problem is fixed, so the change is durable.
- Six skills' model intentions can now actually take effect.
- Pinned version IDs replaced with aliases, removing a documented recurring staleness class.

**Issues**:
- **Tiering alone will not fix the reported problem.** Lifetime Sonnet output was 154K tokens
  against 40.3M on Opus, which means subagents were barely used — so most spend is main-loop,
  which agent frontmatter cannot touch. This change is necessary but not sufficient.
- The dominant term is session length x context size (890 messages x ~110K), not model tier.
- `effortLevel: "high"` compounds it: more thinking -> more output -> larger context -> more
  cache read on every later message. Left unchanged per user decision (investigate first).

**Lessons Learned**:
- A "Current | Target" table describing a platform limitation is a staleness trap; it silently
  becomes wrong when the platform gains the feature. Prefer stating the real mechanism.
- Pin tiers by alias, never by version ID.
- Guidance files propagate defaults. Fixing 27 agents without fixing the sentence that
  produced them would have been undone by the next `skill-editor` run.

## Related Commits

- e29c1cb: config(agents): tier agent models and correct model-dispatch guidance

## Next Steps

1. Re-measure after a few normal sessions. `stats-cache.json` was 12 days stale at the time of
   this change, so establish a fresh baseline before drawing conclusions.
2. Decide on `effortLevel: "high"` -> `"medium"` with that fresh data in hand.
3. Prune the skill library and plugin set. 57 skills cost 10.3K tokens/session in descriptions
   alone; 12 enabled plugins add 8.1K. Anything genuinely unused is a free reduction.
4. Separate follow-up (deliberately not bundled here): `skills/README.md` is stale — claims
   "54 skills (52 active + 2 deprecated)" and "Fable 5 / Opus 4.8" at v3.0.0 dated 2026-06-09,
   against 58 directories and a `claude-opus-5` pin. `CLAUDE.md` says 58. `skills/CHANGELOG.md`
   has carried an `[Unreleased]` heading since 2026-06-09 despite later work landing.
5. Left alone on purpose: `skill-editor/references/knowledge-engineering-report-template.md:443`
   (cosmetic template placeholder) and the model-capability section of
   `anthropic-guidelines-summary.md` (legitimately version-specific, incl. the Fable 5
   life-sciences classifier note).
