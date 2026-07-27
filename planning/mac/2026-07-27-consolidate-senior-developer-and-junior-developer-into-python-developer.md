# Consolidate senior-developer and junior-developer into python-developer

**Date**: 2026-07-27
**Machine**: mac
**Status**: Success

## Objective

An audit of the always-loaded skill-description block (17,515 chars across 55 skills, in
the system prompt of every session) found that 26% of it was "NOT for X, use Y"
disambiguation. That text is a symptom: it exists because skills overlap, and each
NOT-clause is context spent every session working around a taxonomy problem.

Counting which skills get named inside *other* skills' descriptions located the worst
offender precisely. `senior-developer` was cited in 8 other descriptions — no other skill
exceeded 5. Inspecting the cluster showed why:

- `senior-developer` (473 lines) and `junior-developer` (478 lines) were structurally the
  same skill twice, with near-identical headings (Overview, When to Use, Responsibilities,
  Tools, Archival Compliance, Input/Output Format, Pre-Flight Architecture Context,
  Workflow, Progress Reporting, Example).
- The only real distinction was **task scope and autonomy, not capability**: senior took
  complex multi-module components and made component-level design decisions, junior took
  single-responsibility tasks with explicit acceptance criteria and made none.
- Senior's responsibility #5 was "Reviews junior-developer code (max 3 revision cycles)".
  That review hierarchy is precisely the over-verification and over-delegation pattern
  Anthropic's Opus 5 guidance says to remove, and which commit 73b9866 removed everywhere
  else in the repo. The split's main purpose was already obsolete.

Merging them is therefore not just a token optimization — it removes a structure that now
works against the model.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Author a merged `python-developer` from the union of both skills
- [x] Replace the two-role split with one explicit scope contract
- [x] Delete the review hierarchy, revision-cycle limits, and delegation-evaluation step
- [x] Sweep every cross-reference, collapsing the relationship rather than renaming it
- [x] Merge the two agent wrappers; drop the `Task` tool
- [x] Delete the superseded skills and agents from repo and live system
- [x] Verify, push with `--delete`, commit

## Expected Outcome

One developer skill instead of two, no dangling references, and fewer always-loaded
description bytes.

## Actual Outcome

54 skills (from 55) and 27 agents (from 28). The merged `skills/python-developer/SKILL.md`
is 349 lines, replacing 951 across the two originals — a 602-line reduction. 29 files were
touched and roughly 250 mentions rewritten (176 `senior-developer`, 74 `junior-developer`).

Work was split by kind: the merge design and the frontmatter decisions were made directly;
the 951-line union and the 250-mention sweep went to two parallel subagents against written
specs, since both are mechanical but wide.

**The sweep collapsed relationships rather than renaming them**, which was the point. A
naive rename would have produced "python-developer reviews python-developer's code".
Removed outright: the `### junior-developer` optional-specialist section in
`team-composition.md`, the ">3 decomposable tasks → include junior-developer" decision
branch, Step 1b (junior-developer evaluation, ~45 lines including its JSON and `jq`
snippet) from `phase-4-implementation.md`, the junior dispatch template and its
"max 3 revision cycles" review process, the Revision Cycle Protocol in
`code-review-checklist.md`, and the optional `delegation:` block from `handoff-schema.md`.
Two 6-column RACI tables collapsed to 5. Where a dispatch needs to convey task size, it now
carries an explicit `Scope:` field instead of choosing between two roles. Code review moved
to `copilot` plus the objective gates, which is where it belonged.

**A field-name conflict would have silently broken dispatch.** The two skills used
different names for the same concepts. Checking which ones `programming-pm` actually parses
established `code_handoff` and `self_review_checklist` as canonical;
`junior_deliverable`, `self_check`, `files`, `questions`, and `coverage` were junior-only
and are now documented aliases. `junior_task` (referenced by `dispatch-templates.md`)
became `scoped_task` with the old name kept as an accepted alias.

## Assessment

**Result**: Success

**Improvements**:
- One developer skill instead of two near-identical ones; 602 fewer lines.
- The junior→senior review hierarchy is gone from the repo entirely, consistent with the
  Opus 5 retune rather than fighting it.
- Six descriptions that pointed at the deleted skills now point at `python-developer`, and
  `copilot` lost a disambiguation clause describing a hierarchy that no longer exists.
- `senior-developer` is no longer a confusion magnet, because it no longer exists.

**Issues**:
- The instruction to leave frontmatter descriptions alone was wrong for this change. Six
  descriptions named the skills being deleted, so honoring it literally would have left
  dangling pointers. The constraint should have been "don't change `name`, and don't
  gratuitously reword `description`" rather than a blanket freeze.
- Deleting skills needs `push --delete`, which is the one sync operation the repo's own
  anti-pattern list singles out as requiring a dry-run first. The dry-run confirmed exactly
  four files, which is the only reason it was safe to run.

**Lessons Learned**:
- Counting inbound citations in descriptions is a cheap, objective way to find taxonomy
  problems. The skill cited most often by others is the one whose boundary is least clear.
- Overlapping skills cost twice: once in duplicated body content, and again in the
  always-loaded disambiguation needed to tell them apart.
- When merging skills that exchange structured handoffs, resolve field names by checking
  what the *consumer* parses, not by preferring one source. Nothing would have errored —
  it would just have silently stopped matching.

## Related Commits

- [commit SHA]: [commit message]

## Next Steps

Remaining consolidation candidates, from the same citation analysis (counts are inbound
citations in other skills' descriptions):

1. Review/critique cluster — `copilot` (5), `completion-verifier`, `devils-advocate`,
   `edge-case-analyst`. These differ mainly by *when* in the lifecycle they run, which may
   not warrant four skills.
2. Fact-verification cluster — `fact-checker` (5), `essay-fact-checker`,
   `consistency-auditor`.
3. Research cluster — `synthesizer` (5), `lit-pm` (5), `researcher` (4),
   `literature-researcher`, `lit-synthesizer`, `research-pipeline`, `program-officer`.
   Seven skills is where disambiguation concentrates most heavily, but it is also the most
   entangled, so it should go last.

Unrelated follow-ups still open: retire the `software-developer` and `parallel-coordinator`
tombstones (their descriptions cost bytes every session), delete the stray
`~/.claude/skills/patent-review-workspace/` eval workspace, and the items carried over
from the 2026-07-26 entry (35 broken internal links, orphaned Phase 2.5 artifacts,
`bioinformatician/scripts/` documenting two scripts that do not exist).

**Not yet validated in a real session.** This is the third consecutive config change
pushed live without a fresh-session test.
