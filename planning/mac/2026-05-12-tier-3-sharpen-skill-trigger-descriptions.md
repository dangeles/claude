# tier-3-sharpen-skill-trigger-descriptions

**Date**: 2026-05-12
**Machine**: mac
**Status**: Success

## Objective

Tier 3 of the improvement pass. Sharpen `description:` frontmatter on 10 skills that compete on similar-sounding prompts so dispatch picks the right one without ambiguity. Each new description adds explicit "NOT for X (use Y)" exclusion clauses pointing at sibling skills.

Skills targeted (audit findings from Tier 2 + carry-over from cancelled Tier 2 B2/B3):
- **C1 Orchestrators**: technical-pm, programming-pm, lit-pm, research-pipeline
- **C2 Quant specialists**: statistician, mathematician, calculator
- **C3 Review lanes**: copilot, senior-developer
- **C4 Edge-case scope**: edge-case-analyst

## Changes Planned

- [x] technical-pm: clarified as default fallback; excludes domain PMs
- [x] programming-pm: added Python + bioinformatics specialist chain; excludes ad-hoc fixes and non-Python work
- [x] lit-pm: emphasised 4-24h adaptive 9-stage; excludes quick chains and single lookups
- [x] research-pipeline: emphasised 2-8h fixed 5-stage; excludes adaptive reviews and custom sequences
- [x] statistician: stochastic/probabilistic; excludes deterministic algorithms and physics
- [x] mathematician: deterministic/numerical; excludes stochastic and physics
- [x] calculator: physics/engineering; excludes statistics, complexity, and cost
- [x] copilot: inline dev review only; excludes pre-PR, post-PR, and pipeline review
- [x] senior-developer: implementation + formal junior review in pipeline; excludes general review
- [x] edge-case-analyst: pre-implementation only; excludes debug and post-hoc review
- [x] Validate YAML frontmatter for all 10 (2 required quoting after revision due to internal `: `)
- [x] Sync push --yes
- [x] Commit + SHA backfill

## Expected Outcome

Cleaner dispatch behavior on ambiguous prompts:
- "Write a literature review on X" → lit-pm (comprehensive) or research-pipeline (quick) based on user signals
- "Estimate the diffusion limit" → calculator only
- "Run a power analysis" → statistician only
- "Review this code change" → copilot for inline, /pr-review-toolkit:review-pr for pre-PR, /code-review for post-PR, senior-developer for pipeline junior review
- "Identify failure modes for the new feature" → edge-case-analyst (not systematic-troubleshooter)

## Actual Outcome

All 10 description frontmatter rewrites applied. YAML validation caught 2 instances where prose `: ` (colon-space) inside description text needed quoting; both fixed with `"..."` wrapping. Sync push --yes overlaid 10 SKILL.md files; live skill list now reflects the new triggers.

## Assessment

**Result**: Success

**Improvements**:
- Explicit "NOT for X" clauses give Claude concrete redirect signals during dispatch — better than implicit semantic similarity heuristics
- Pairs that were previously discriminated only by body content (which doesn't auto-load) now discriminate at the description level (which always loads)
- Sets a template for future skill additions: every description should answer "and what should NOT trigger this?"

**Issues**:
- YAML parser caught `: ` errors at validation; this surfaces a pattern the frontmatter-validation push-time gate should also check (it might already; it caught nothing because we ran jq separately).

**Lessons Learned**:
- When writing descriptions in YAML, avoid `: ` inside unquoted strings or wrap in `"..."`. A pre-commit hook could enforce this.

## Related Commits

- [pending]: tier-3 sharpen skill trigger descriptions

## Next Steps

- Optional: add a frontmatter linter to catch `: ` in unquoted description values
- Continue refining triggers as misdispatches surface in real sessions
