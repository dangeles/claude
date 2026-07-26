# Change-Scale Signals

Signals for choosing how much review a skill change warrants. Used by skill-editor's
"Choosing How Much Review to Run" section. These are heuristics, not gates — when they
disagree with the obvious reading of the request, ask the user.

## Small — edit directly, then run Phase 4

- Documentation, typo, comment, or example changes: single file, ≤50 lines
- Minor bug fix: single file, ≤50 lines
- Description or trigger-phrase wording, where the change does not alter what the skill does

## Substantive — light review (spec if ambiguous, short plan, Phase 4)

- One skill, 2 files or fewer, roughly 50-150 lines
- A reworked workflow step, a new section, changed triggers
- Contains workflow keywords ("agent", "workflow", "phase", "quality gate") but the change
  itself is small

## Substantial — full panel

- New skill creation (architectural addition)
- More than 4 files affected
- More than about 300 lines changed
- Refactoring keywords ("refactor", "reorganize", "restructure", "migrate") combined with
  more than 2 files or more than 150 lines
- Changes to workflow phases, quality gates, or agent wiring
- The user asks for strategic review, architectural assessment, or adversarial review

## Borderline

3-4 files, or 150-300 lines, or refactoring keywords at moderate scope, sit between the
last two tiers. Say which way you are leaning and why, and let the user redirect. A change
that is large in line count but mechanical (a rename across files) is not the same risk as
a change that is small but rewires control flow — weight structural risk over size.

## Experimental / fast iteration

When the user explicitly asks for a quick or experimental pass ("try this", "prototype",
"quick change"), skip the analysis phases regardless of size, and say that you did so.
Phase 4 still runs — the sync sequence is what protects the live config.
