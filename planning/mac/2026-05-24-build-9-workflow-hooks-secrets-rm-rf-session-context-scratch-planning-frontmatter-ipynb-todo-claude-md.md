# Build 9 workflow hooks: secrets/rm-rf/session-context/scratch/planning/frontmatter/ipynb/todo/claude-md

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

After landing the `session-scope-guard` hook, the user asked for additional workflow hooks tailored to their actual repo inventory (28 repos spanning research / writing / business / finance / personal). I inspected `~/repos/` to ground recommendations, then built three tiers of nine hooks total.

## Inventory signals that shaped the design

- 28 active repos; 44 notebooks; 63 data files >1MB versioned with `!path` gitignore negations (sophisticated data-versioning workflow — don't fight it)
- No `.env`-style files at top level but ~10 source files contain secret-shaped patterns (some prose, some real risk surface)
- A dedicated `~/repos/scratch/` repo exists as the deliberate dumping ground
- `claude-config/` has a planning-entry workflow that's currently advisory (CLAUDE.md prose, not enforced)
- `sync-config.py push` already gates SKILL.md frontmatter at push time — worth catching upstream too

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Tier 1 — `secret-leak-guard.py` (PreToolUse Write/Edit/NotebookEdit), `rm-rf-sanity-guard.py` (PreToolUse Bash), `session-start-context.py` (UserPromptSubmit, once per session)
- [x] Tier 2 — `scratch-redirect-suggest.py` (PostToolUse Write), `claude-config-planning-required.py` (Stop, scoped to claude-config repo), `skill-frontmatter-validator.py` (PreToolUse Write/Edit)
- [x] Tier 3 — `ipynb-output-warn.py` (PreToolUse Bash), `stale-todo-tracker.py` (Stop), `claude-md-improver-suggest.py` (Stop)
- [x] Wire all 9 hooks in `claude-config/settings.json` under appropriate matchers
- [x] Smoke-test each hook with synthetic stdin payloads
- [x] `./sync-config.py push` and verify all 9 hooks live under `~/.claude/hooks/`

## Implementation

### Tier 1 — high-value safety

**`secret-leak-guard.py`** (PreToolUse on Write / Edit / NotebookEdit, ~100 lines)
- 12 patterns: Anthropic / OpenAI / GitHub PAT+fine-grained+OAuth / AWS access+secret / Slack / Google / Stripe / JWT / PEM private key / hard-coded password assignments
- The password-assignment heuristic excludes obvious placeholders (`<...>`, `password`, `changeme`, `TODO`, `REPLACE`, etc.)
- Exit 2 on match with stderr listing each finding (excerpt truncated for long matches); exit 0 otherwise
- Bypass: `SKIP_SECRET_LEAK_GUARD` in description OR env var

**`rm-rf-sanity-guard.py`** (PreToolUse on Bash, ~110 lines)
- Matches `rm -r*` / `rm -*r*` / `find ... -delete`
- Refuses paths that are: empty, bare `*` / `.` / `./`, `$HOME` exactly, `~/repos` exactly, a top-level repo directory (e.g. `~/repos/claude`), filesystem root, or outside the current git repo's root
- Splits multi-statement commands on `\n` / `&&` / `||` / `;` and inspects each independently
- Bypass: `SKIP_RM_GUARD`

**`session-start-context.py`** (UserPromptSubmit, ~90 lines)
- One-shot per session via a marker file in `.claude/session-state/<session_id>-start.marker`
- On the first prompt, prints to stdout (which Claude Code prepends to the prompt): cwd, branch, upstream + ahead/behind counts, dirty-files count, last 2 commits, time-since-last-commit
- Exit 0 always — informational only

### Tier 2 — workflow-specific

**`scratch-redirect-suggest.py`** (PostToolUse on Write, ~60 lines)
- Non-blocking. When a file with a scratch-y name (`temp_*`, `scratch_*`, `wip_*`, `debug_*`, `untitled_*`, `throwaway_*`, etc.) is written to a repo's root, emit a stderr suggestion pointing at `~/repos/scratch/`
- Only fires for top-level writes (not buried in subdirs), and skips when the repo IS `~/repos/scratch`

**`claude-config-planning-required.py`** (Stop hook, ~110 lines)
- Scoped to the `claude` repo (the claude-config repository) via `git rev-parse --show-toplevel | basename`
- Walks the session transcript for any Write/Edit/NotebookEdit on paths matching `claude-config/(skills|agents|hooks|plugins)/` or `claude-config/settings.json`
- If touched but no `planning/<host>/*.md` file was Written this session, exit 2 (Stop hook semantics: refuse to stop, force continuation) with a reminder to run `./sync-config.py plan`
- Bypass: `SKIP_PLANNING_REQUIRED` env var or substring in any user message

**`skill-frontmatter-validator.py`** (PreToolUse on Write/Edit, ~110 lines)
- Triggers only for paths matching `claude-config/skills/*/SKILL.md`
- For Write: validates the content's frontmatter. For Edit: reads the on-disk file, applies the substitution, validates the result.
- Requires `name:` and `description:` fields. Uses PyYAML if available for full validation; falls back to a regex scan for key names if not.
- Exit 2 if invalid, with the specific reason
- Bypass: `SKIP_FRONTMATTER_VALIDATION`

### Tier 3 — observational

**`ipynb-output-warn.py`** (PreToolUse on Bash, ~35 lines)
- Non-blocking. Detects `git add ... *.ipynb` patterns and emits a stderr nudge toward the jupytext markdown workflow
- Suppress: `noqa-ipynb-warn` in description

**`stale-todo-tracker.py`** (Stop hook, ~95 lines)
- Walks the transcript for every Write/Edit/NotebookEdit and counts `TODO|FIXME|XXX|HACK` markers added (in `content` / `new_string` / `new_source`) and removed (in Edit `old_string`)
- Subtracts overlapping (file, line) tuples — Edits that move a marker count as no net change
- Emits a session summary to stderr with top 10 additions and top 5 removals
- Exit 0 always

**`claude-md-improver-suggest.py`** (Stop hook, ~90 lines)
- Scans user-role messages in the transcript for durable-preference phrases: `from now on`, `next time`, `remember to`, `always (use|prefer|do)`, `never (use|do|commit|push)`, `don't (use|do|forget)`, `please (always|never|don't)`, `going forward`, `make (a|this a) (rule|habit|convention)`
- If any hits, suggest `claude-md-management:revise-claude-md` to capture the convention
- Dedupes by excerpt; shows up to 3 examples
- Exit 0 always

### Wiring (`settings.json`)

Hooks added under four hook types:

| Type | Hooks (in order) |
|---|---|
| PreToolUse Bash | session-scope-guard, rm-rf-sanity-guard, ipynb-output-warn |
| PreToolUse Write | secret-leak-guard, skill-frontmatter-validator |
| PreToolUse Edit | secret-leak-guard, skill-frontmatter-validator |
| PreToolUse NotebookEdit | secret-leak-guard |
| PostToolUse Write | scratch-redirect-suggest |
| UserPromptSubmit | session-start-context |
| Stop | claude-config-planning-required, stale-todo-tracker, claude-md-improver-suggest |

Total: 10 hook instances across 7 unique scripts plus the existing session-scope-guard (carried over from `dee3a86` / `740204a`).

## Expected Outcome

- Three blocking hooks (secret-leak / rm-rf-sanity / skill-frontmatter / session-scope-guard) intercept common destructive or risky operations.
- Three informational hooks (scratch-redirect / ipynb-output-warn / scratch suggestions) nudge without enforcing.
- Three end-of-session hooks (claude-config-planning-required / stale-todo-tracker / claude-md-improver-suggest) catch process gaps that are hard to notice mid-flight.
- `session-start-context` gives the model situational awareness from the first prompt without explicit explanation.

## Actual Outcome

All 9 hooks tested cleanly with synthetic stdin payloads. Live state in `~/.claude/hooks/`:

```
$ ls ~/.claude/hooks/
claude-config-planning-required.py    secret-leak-guard.py
claude-md-improver-suggest.py         session-scope-guard.py
ipynb-output-warn.py                  session-start-context.py
rm-rf-sanity-guard.py                 skill-frontmatter-validator.py
scratch-redirect-suggest.py           stale-todo-tracker.py
```

Smoke-test results (verified independently from `head`-truncated pipelines, which mis-reported exit codes):

- `secret-leak-guard`: blocks fake Anthropic key (exit 2), allows benign content (exit 0)
- `rm-rf-sanity-guard`: blocks `rm -rf ~` (exit 2), allows `rm -rf scratch/build` inside a repo (exit 0)
- `session-start-context`: emits compact `<session-context>` block with cwd / branch / dirty / last commit (exit 0)
- `scratch-redirect-suggest`: emits suggestion when `temp_test.py` is written to repo root (exit 0)
- `claude-config-planning-required`: blocks Stop (exit 2) when claude-config touched but no planning entry created
- `skill-frontmatter-validator`: blocks missing `description:` (exit 2), allows valid frontmatter (exit 0), no-op on non-SKILL paths (exit 0)
- `ipynb-output-warn`: emits warning on `git add foo.ipynb` (exit 0)
- `stale-todo-tracker`: reports new TODO markers from synthetic transcript
- `claude-md-improver-suggest`: matches `remember to` and suggests claude-md update

## Assessment

**Result**: Success

**Improvements**:
- The catalog has three distinct hook flavors now: structural blocks (session-scope, secrets, rm-rf), upstream validators (frontmatter), and observational nudges (session-context, scratch, ipynb, todo, claude-md, planning). The pattern is generalizable.
- Each hook is self-contained — no shared module — which makes them independently testable and reviewable. Tradeoff is some duplication of patterns like "walk transcript for tool_use events"; if a third use case appears, a `_hook_lib.py` extraction makes sense.
- Sync infrastructure already in place from the previous hooks work (`hooks/` in `sync_rules.always`), so adding the 9 was a no-op for `sync-config.py`.

**Issues**:
- `session-start-context.py` writes a marker file under `$CLAUDE_PROJECT_DIR/.claude/session-state/`. If `CLAUDE_PROJECT_DIR` isn't set, it falls back to `cwd`. The marker filename uses `session_id` so concurrent sessions don't collide. Worth keeping an eye on.
- The `secret-leak-guard.py` `password` heuristic regex can hit overlapping matches with the API-key patterns (a fake Anthropic key matches both "Anthropic API key" and "OpenAI API key (sk-)" and "hard-coded password"). Cosmetic; the block decision is correct.
- `claude-md-improver-suggest.py` triggers on phrases that might appear conversationally without being a durable preference. The dedupe + 3-example cap keeps it from being noisy, but expect occasional false positives.

**Lessons Learned**:
- macOS `echo` and `zsh` interpret `\n` inside single-quoted JSON test strings, which breaks shell-piped tests. Use `python3 -c "import json; print(json.dumps(...))"` to generate test JSON instead.
- `head` in a pipeline captures the exit code of `head`, not the upstream command. For exit-code-sensitive smoke tests, drop the `head` or use `${PIPESTATUS[0]}`.
- Hooks that emit stderr-only output with exit 0 are a clean way to do informational notices that don't block but still surface in the model's view.

## Related Commits

- `<pending>`: feat(hooks): build 9 workflow hooks (secret-leak, rm-rf-sanity, session-start-context, scratch-redirect, planning-required, skill-frontmatter, ipynb-output-warn, stale-todo, claude-md-improver-suggest)

## Next Steps

- If two more hooks adopt the transcript-walking pattern, extract a shared `_hook_lib.py` helper.
- Consider a `--diagnose` flag on each blocking hook that prints what it would do given an explicit `--command` arg, useful for debugging false positives without crafting JSON.
- The `claude-md-improver-suggest` patterns could be tightened or expanded based on real-world false-positive rate; revisit after a week of use.
- If sessions across multiple projects pile up `session-state/*.marker` files, add a periodic cleanup (>30 days old).
