# Delegation and scope

Shared conventions for orchestrator skills. Tuned for Claude Opus 5, which delegates
to subagents more readily than earlier models and writes longer output by default.

## Contents

- Delegating to subagents
- Verification
- Length of written deliverables

## Delegating to subagents

Subagents multiply cost and time: each one re-establishes context, re-explores, and
reports back, and you then re-read its report. Delegate when the payoff clearly exceeds
that overhead.

Do delegate:

- Sizeable, genuinely independent tracks — unrelated modules, a wide multi-file
  investigation, several items that can be worked in parallel.

Do not delegate:

- Work you could finish yourself in a handful of tool calls: a few file reads, a
  handful of edits, a simple search, a quick check.
- Review or verification of your own work. That belongs in your main loop.

On parallel use:

- If one subagent can do the job, use one. Keep spawn counts low.
- Don't split one modest job across several subagents. Parallel fan-out is for
  independent, sizeable tracks.
- Brief each subagent precisely the first time rather than launching, waiting, and
  re-briefing.
- Once you delegate, commit to it. Don't redo a subagent's work or re-derive its
  findings after it reports.
- When launching several agents for independent work, send them in a single message so
  they run concurrently.

## Verification

Opus 5 verifies its own work without being told to. Instructions to double-check,
re-verify, or add a final verification pass cause over-verification without improving
correctness, so this file deliberately contains no verification checklist.

Verify where there is an objective check to run — a test suite, a validator, a
compiler, a schema. Skip the self-review prose pass.

## Length of written deliverables

Match the length of files you write — reports, Markdown documents, summaries — to what
the task needs. Cover the substance without padding documents with filler sections,
redundant summaries, or boilerplate. This applies to artifacts on disk, not just to
chat responses.
