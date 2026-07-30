# Add paper-reader skill and PDF context guard hook

**Date**: 2026-07-30
**Machine**: mac
**Status**: Complete

## Objective

Asking Claude to read a paper produced three bad behaviors: the full text went into the main
context window, the notes were a table of extracted statistics, and the biology went
unmentioned.

None of that was Claude improvising. There was no paper-reading skill, so "read this paper"
fell through to `researcher` — a skill built for a different job. Its stated purpose is
sourcing "quantitative parameters from primary literature," its paper-notes template puts a
**Quantitative Values** table at the center with biology demoted to a one-line "Relevance to
Project," and its self-check is literally "scan your document for numbers, percentages,
rates." The model optimized exactly what the skill measured.

The context problem had the same root: `researcher` says acquire PDFs and extract tables,
with no instruction to do it out-of-context, and no delegation guidance at all. A full text
is 8–15k tokens; three of them cost more than the conversation that follows.

`biologist-commentator` is the closest thing to a "what does this mean" skill, but its
description scopes it to evaluating bioinformatics *approaches*, so it never triggers on
"summarize this paper."

## Changes Planned

- [x] New skill `paper-reader`: comprehension over extraction — what question, what they
      actually did, what they found vs. what they claim, does the evidence carry the claim,
      what it means mechanistically, what would change your mind
- [x] Explicit anti-goal in the skill: do not default to extracting statistics; numbers
      appear only when they carry the argument, always with measurement context
- [x] Context rule: one subagent per paper by default, notes back to the main thread; read
      inline only for short pieces or a known passage
- [x] New hook `paper-context-guard.py` (PreToolUse on `Read|WebFetch`): advisory warning
      when a local PDF or a publisher/preprint host URL is pulled into the main context
- [x] Hook registered in `settings.json` under a `Read|WebFetch` matcher
- [x] `tests/test_paper_context_guard.py` — 10 cases
- [x] `researcher` description gains a "NOT for" pointing at paper-reader
- [x] `skills/README.md` row; `CLAUDE.md` skill count 54 → 55

## Expected Outcome

"Read this paper" routes to a skill built for understanding rather than extraction, and full
texts stop landing in the main context.

## Actual Outcome

5 files synced, 39 tests pass (10 new). Live hook verified executable and firing correctly
after push. Skill and hook both work; the live skill list reloaded showing paper-reader and
the updated researcher description.

Hook design notes worth keeping:
- **Advisory, never blocking** (exit 0 always). Reading a paper inline is legitimately right
  sometimes; a blocking hook would be wrong and would train the user to bypass it.
- **Host matching is suffix-anchored.** `biorxiv.org.attacker.com` must not inherit
  biorxiv.org's match. Pinned by test.
- **Subagents are exempt.** A subagent reading the paper is the recommended pattern, so
  warning there would be backwards. Detection is best-effort across payload fields and
  `CLAUDE_AGENT_TYPE`, defaulting to warn — a spurious warning in a subagent is cheap, a
  missed one in the main thread is the case that matters.

## Assessment

**Result**: Success

**Improvements**:
- The routing fix matters as much as the new skill. `researcher` was capturing a request it
  was never designed for, and no amount of instruction inside it would have helped.
- Encoding the context rule in *both* the skill (as policy) and the hook (as a reminder at
  the moment of action) covers the case where the skill wasn't invoked at all.

**Issues**:
- Subagent detection relies on undocumented payload fields. If Claude Code changes them, the
  hook degrades to warning inside subagents — noisy but harmless, and the `noqa` escape
  hatch still works.
- The publisher host list is manually curated and will miss journals. Local PDFs are caught
  regardless of source, which covers most real cases.

**Lessons Learned**:
- When a skill behaves badly, check first whether the *right* skill exists. Two of the three
  complaints here were a routing problem wearing a behavior problem's clothes.
- A skill's stated success criteria are what the model optimizes. "Scan your document for
  numbers" reliably produces documents full of numbers.

## Related Commits

- [commit SHA]: [commit message]

## Next Steps

- Field-test both changes: write one real essay, read one real paper.
- Sweep the other same-era orchestrators (`lit-pm`, `research-pipeline`, `program-officer`)
  for the same over-specification found in essay-pipeline.
- Consider whether `researcher`'s own paper-notes template should drop the Quantitative
  Values table as its centerpiece, now that paper-reader owns comprehension.
