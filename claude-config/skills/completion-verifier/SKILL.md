---
name: completion-verifier
version: 1.0
last_updated: 2026-07-26
description: "Use BEFORE marking any task complete, before user handoffs, or at quality checkpoints — verifies requirements satisfied, edge cases tested, tests pass, no regressions, deliverables ready. Triggers on 'are we done', 'is this complete', 'verify completion', 'ready to ship', 'ready to merge', 'is this ready', 'mark done'. NOT for general code review (use copilot) or pre-PR comprehensive review (use /pr-review-toolkit:review-pr)."
success_criteria:
  - All stated requirements verified as satisfied
  - Edge cases systematically identified and tested
  - Quality standards met for domain (tests pass, docs complete, etc.)
  - No regressions introduced by the work
  - Deliverables ready for handoff or production use
  - Verification checklist completed with evidence
---

# Completion Verifier

Run the objective checks that exist for this work before calling it done.

## Scope

This skill covers checks that produce evidence: a command that passes or fails, a file
that exists or doesn't, output that matches a spec or doesn't. It deliberately does not
prescribe a self-review reading pass — re-reading your own work as prose is not a check,
and current models already do it unprompted.

If there is nothing objective to run, there is nothing for this skill to do. Say what
was done and hand off.

## When to use

Before marking a task complete or handing back to the user on work that has runnable
checks — a test suite, a validator, a linter, a build, a schema, a citation that resolves
to a real source.

Skip it for trivial edits, progress updates, and work whose only verification would be
re-reading what you just wrote.

## Checks

Run what applies and report the actual result, including failures.

**Tests.** Run the suite, not a subset you expect to pass. Report the command and its
output. A test you did not run is not a test that passed.

**Regressions.** Confirm previously-passing checks still pass. For code, this usually
means the full suite rather than only tests touching changed files.

**Requirements.** Compare deliverables against what was explicitly asked for. Note
anything requested that is missing, and anything delivered that was not requested.

**Edge cases.** Where the domain has boundary conditions with a runnable check — empty
input, malformed input, missing file, zero rows — exercise them. Where they cannot be
run, note them as untested rather than assumed handled.

**Deliverables.** Confirm the artifacts the task promised exist at the paths given, and
are non-empty.

**Documentation.** Where the change alters an interface others use, confirm the
corresponding docs reflect it.

Domain-specific criteria for code, research, analysis, and documentation work:
see `references/completion-criteria-by-domain.md`.

## Reporting

State what you ran, what passed, what failed, and what you could not check. Failures and
gaps go in the report with the same prominence as passes — a verification report that
only lists passes is not a report.

If checks fail, fix them or say plainly what is broken and why. Do not describe work as
complete when a check is failing.

## Integration

Incoming: any skill or agent finishing a unit of work that has runnable checks.

Outgoing: the caller that requested verification, or the user. Report results; the
decision to ship is theirs.
