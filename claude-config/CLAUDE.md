# Working policy

Loaded in every session, in every repo. Kept short on purpose: it is re-read on
every request, so it holds policy only. Project-level `CLAUDE.md` files and direct
instructions override anything here.

## Read targeted, not whole

Reading a file puts it in context permanently, and context is re-read on every
later request in the session. A 2000-line file read once to check one function is
paid for on every turn that follows.

- Know the region you need? Read only that region — `Read` takes `offset`/`limit`.
- Looking for something? `Grep` with `output_mode: "content"` and `-A`/`-B` beats
  reading candidate files whole. `Glob` to locate, then a narrow read.
- Read the function, not the module. Read the changed hunk, not the file.
- Never re-read a file already in context, and never re-read a file you just
  edited to confirm the edit — `Edit`/`Write` error if they fail.
- Whole-file reads are right when you will genuinely use most of it: a file you
  are about to rewrite, a short config, a file under ~200 lines.

## Delegate to protect context

The same logic scales up. A subagent that reads twenty files and returns a summary
pays for those files once; doing it in the main loop pays for them on every
subsequent turn. Cache reads dominate cost.

Dispatch a subagent for:

- Answering one question that needs more than ~5 files
- Any repo-wide sweep, audit, or broad search
- Reproducing or bisecting a failure
- Reading a full-text paper or other large document

Keep in the main loop: the decision, the edit, the conversation with the user.
Do not delegate work that is a couple of tool calls, and do not delegate the
judgment itself.

## Pass `model` on every dispatch

`general-purpose` has no model of its own and inherits the session's, which means
Opus unless told otherwise. Set it explicitly:

- `sonnet` — search, review, analysis, structured generation
- `haiku` — mechanical scans and audits against explicit rules
- omit only when the subagent genuinely needs session-level capability

An agent dispatched N-at-a-time multiplies its tier N times. Weigh that first.

## Plan before multi-file work

For anything touching 3+ files, write the plan first, then dispatch the
independent pieces. Explore-then-edit incrementally in one long main-loop session
is the most expensive way to do it and the hardest to review.

## Budget

Past roughly 150K of context, finish the current step and suggest `/clear` rather
than starting something new in the same session.
