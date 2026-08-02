# genesis-design skill, token guard hook, genesis-check command

**Date**: 2026-08-02
**Machine**: mac
**Status**: Complete

## Objective

`~/repos/aesthetic_and_web_guidelines` is the source of design truth for every Genesis
surface. It ships guidance and modifies nothing, and its distribution model is — awkwardly —
copying, which is the same mechanism it exists to argue against. An agent asked to build a
Genesis surface either does not know the design system exists, or reads nineteen thousand
words to answer one question.

Close that gap with a routing skill, a nudge hook that makes the skill fire deterministically
on surface CSS, and a slash command wrapping the canonical checker.

Built from the spec at `~/repos/aesthetic_and_web_guidelines/research/skill-brief-2026-08-02.md`,
via the `skill-editor` workflow (session `session-20260802-164130-32494`).

## Changes Planned

- [x] `claude-config/skills/genesis-design/SKILL.md` — routing skill
- [x] `claude-config/hooks/genesis-token-guard.py` — non-blocking PreToolUse nudge
- [x] `claude-config/commands/genesis-check.md` — first slash command in this repo
- [x] `sync.config.yaml` — `commands/` added to `sync_rules.always`
- [x] `claude-config/hooks/claude-config-planning-required.py` — gate now covers `commands/`
- [x] `CLAUDE.md` — `commands/` and `hooks/` rows; skill count 55 → 58
- [x] `claude-config/settings.json` — hook registered on `Write` and `Edit`

## Expected Outcome

An agent editing `~/repos/papers/web/static/library.css` learns the design system governs it,
without being told, and reaches the rules by reading them rather than recalling them.

## Actual Outcome

Shipped as planned, with **one material deviation from the brief**.

Brief §4 specified seven content-inspection branches in the hook, each citing a rule ID
(COLOR-03, COLOR-02, TYPE-01, A11Y-03, A11Y-04, MOTION-01). Adversarial review measured them
against the real surface files and the design repo's own `check_surface.py`, and they did not
survive:

- The weight check produced **six false positives against six true findings** in
  `library.css` — the very file the brief names as "the single most likely place for this
  whole effort to be wasted." Cause: SN Leif genuinely is OS/2 class 600, and
  `check_surface.py:278` exempts `400` and `600` by name. A hook citing TYPE-02 on
  `font: 600 var(--t-body)/1.4 var(--serif)` would report a violation the canonical tool
  denies.
- `COLOR-02` (`--gx-` outside a token block) matched **zero** occurrences across all three
  surface repos — dead code.
- `COLOR-03` would have fired on the 38–48 hex literals in each surface `tokens.css`, all of
  which are correct there.
- `TYPE-01` matched `font-size: var(--t-*)` 44× in white_papers and 46× in kol — ~97% false.

The deeper problem was not the regexes. It was that seven branches reimplement rules owned by
another repository, in a third repository, with no shared fixtures pinning the two together —
which is verbatim the failure the design system exists to prevent, and which brief §7's own
first bullet forbids ("must not restate the guidelines"). Two independent attempts at the
weight regex produced two different wrong answers; that was the signal.

**Resolution (user-approved): cut to two branches.** `tokens.css` → ADOPT-01, which is a
filename test and cannot drift; and the once-per-session-per-repo pointer, which routes to
`check_surface.py`. The checker is canonical and runs in ~125 ms on a whole repo, so there is
no performance argument for pre-judging. The hook lost two thirds of its size, all of its
measured false positives, and none of its purpose — the deterministic skill trigger is intact,
since every message ends with "Invoke the `genesis-design` skill before continuing."

Three smaller review findings, all fixed:

- The path gate matched `genesis` as a **substring**, which captures
  `~/repos/organogenesis_notion/` (185 files). Now matched as a path component, with a
  stylesheet extension required on every clause.
- The skill was to pin the design repo's HEAD SHA. That repo has **one commit** and 32 dirty
  entries; `docs/rules.md`, `docs/INDEX.md`, `check_surface.py` and `genesis.rules.json` are
  all **untracked**, so rule count at HEAD is 0, not 42. The pin would have cried wolf on
  every invocation while asserting freshness about a tree containing none of the documents.
  Replaced with a dated rule-count sentence phrased as an instruction, plus a preflight
  warning that those files are untracked.
- `check_surface.py` handed a single **file** globs zero stylesheets and reports
  `{"ok": true, "error": 0}` — indistinguishable from clean. `/genesis-check` now resolves to
  a repository root via `git rev-parse --show-toplevel`.

## Assessment

**Result**: Success

**Improvements**:
- The design system is now reachable from a Claude session without knowing it exists.
- `commands/` is a synced surface for the first time, and the planning gate covers it — the
  gap was created and closed in the same change.
- `session-state/` exclusion (committed separately, `70be339`) fixes a real pre-existing
  sync bug found while deciding where hook state could live.

**Issues**:
- **Nothing gates surface CSS today.** None of `papers`, `white_papers` or `kol` has a
  pre-commit hook, and the guidelines repo's own hook runs only the inward checks
  (`check_tokens.py`, `palette_audit.py --selftest`) — never `check_surface.py`. The brief
  assumed otherwise. All 37 current findings are in surfaces.
- `noqa-genesis` reads `tool_input.description`, which `Write`/`Edit` may not carry. It is
  the house pattern so it stays, but the escape hatch may be unreachable in practice. Low
  cost now that the hook only ever emits a routing pointer.
- `font-variation-settings: "wght" 500` evades both this hook and `check_surface.py`. Worth
  reporting upstream.

**Lessons Learned**:
- A detailed, confident brief is still a hypothesis. Every quantitative claim in this one was
  checkable in about a minute, and the two that mattered most — which rules actually fire,
  and what the canonical checker exempts — were the ones that changed the design.
- When a hook would duplicate logic owned elsewhere, ask what pins the copies together. If
  the answer is "nothing," the hook should route instead of judge. `skill-frontmatter-validator.py`
  duplicates a rule deliberately and `tests/test_frontmatter_contract.py` pins it; no such
  test was possible across repository boundaries here.
- `push --dry-run` validates skill/agent frontmatter and JSON syntax only. It does not check
  the exec bit, that a hook is valid Python, or that its path resolves. Those need explicit
  pre-push commands.

## Related Commits

- `70be339`: fix(sync): exclude session-state/ from project config sync
- `bb4582d`: feat(config): add genesis-design skill, token guard hook, and /genesis-check
- `92d883e`: feat(config): register genesis-token-guard on Write and Edit

## Next Steps

- **Correct long-term intervention, out of scope here**: generalise
  `~/repos/aesthetic_and_web_guidelines/tools/install_hooks.py` (it hardcodes `ROOT` to its own
  repo) to install a `check_surface.py` pre-commit hook into the three surface repos. Right
  layer, canonical checker, zero machine-wide cost. That is a change to the guidelines repo and
  the surface owner's decision.
- Revisit the `Stop` hook idea only if the nudge turns out to be ignored.
- Re-run brief §8 acceptance tests 1, 2 and 5 in a genuinely fresh session — triggering
  behaviour cannot be self-assessed from the session that authored the skill.
