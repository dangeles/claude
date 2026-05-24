# Slim git-strategy-advisor + add session-scope-guard hook

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

User pushback on the prior change to `git-strategy-advisor`: the 405-line skill (4 phases, decision matrices, confidence calibration, YAML schemas, error-handling protocol) was overengineered for solo work where complex git decisions are rare. They wanted a lightweight tool that quickly decides what git actions to take.

Two parts:

1. **Slim the skill aggressively** (~405 → ~80 lines): drop the phase machinery, decision matrices, output schemas, error handling, and integration notes. Keep the sharpened description (it earns its weight in the always-loaded index) and replace the body with a heuristic table + a terse output shape + a worked example.
2. **Add a structural pre-commit guard hook**: a Python script wired as a PreToolUse hook on the Bash tool. Refuses bulk-stage patterns (`-A` / `.` / `-u`) outright. Refuses `git add <path>` or `git commit` when the staged path wasn't modified during the current session (derived from the conversation transcript). Replaces the advisory session-scope discipline that lived in the skill body with a structural enforcement layer.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Rewrite `claude-config/skills/git-strategy-advisor/SKILL.md` (405 → 99 lines): drop 3-phase workflow, decision matrices, confidence calibration, YAML output schema, error-handling protocol, integration notes; keep heuristic table + worked example + terse output shape
- [x] Delete `claude-config/skills/git-strategy-advisor/references/` and `examples/` (supported the deleted machinery)
- [x] Write `claude-config/hooks/session-scope-guard.py` (executable Python, PreToolUse hook)
- [x] Wire the hook in `claude-config/settings.json` under `hooks.PreToolUse` matching `Bash`
- [x] Add `hooks/` to `sync.config.yaml`'s `always` sync rules (without it the hook file doesn't reach `~/.claude/`)
- [x] Smoke-test the hook locally with synthetic stdin payloads (7 cases)
- [x] `./sync-config.py push --yes` and verify live state at `~/.claude/hooks/session-scope-guard.py`

## Implementation

### Slim skill (99 lines, was 405)

Kept the sharpened description (783 chars — same triggering surface, same four NOT-for clauses routing to `/commit-commands:commit`, `/commit-commands:commit-push-pr`, `superpowers:finishing-a-development-branch`, `superpowers:using-git-worktrees`). The description still mentions session-scope but now defers structural enforcement to the hook.

Body now consists of:

- One-paragraph framing of the skill as "lightweight nudge for non-routine git decisions"
- Brief When-to-Use (3 bullets) and When-NOT-to-Use (5 bullets) lists
- One-line session-scope rule pointing at the hook for enforcement
- A 6-row heuristic table covering the realistic decision space (small commit on main, single-topic on main, multi-topic, WIP, feature-ready, out-of-session dirty tree)
- Branch-naming rule (when branching is the answer)
- Terse output template (no YAML schema)
- One worked example
- Three closing notes

Dropped entirely: `Phase 1-3` workflow, decision-matrix references, confidence-calibration table, `OUT_OF_SESSION_CHANGES` / `STAGED_OUT_OF_SESSION` warning schema (now lives in the hook's stderr output), `Integration Notes` mentioning six orchestrators, error-handling protocol.

### Hook (`session-scope-guard.py`)

Python 3 script invoked as a PreToolUse hook on the `Bash` matcher. Reads the standard PreToolUse JSON payload from stdin (`tool_name`, `tool_input.command`, `tool_input.description`, `transcript_path`).

Logic:

1. **Bypass checks**: `SESSION_SCOPE_BYPASS=1` env var or `SKIP_SESSION_SCOPE_GUARD` in description → exit 0.
2. **Non-Bash or non-git Bash** → exit 0 (no opinion).
3. **Bulk-stage patterns** (`-A` / `--all` / `.` / `-u` / `--update` on `git add`) → exit 2 (block) with stderr message, regardless of session-set status. These are unsafe by construction.
4. **Session-set derivation**: walk the transcript JSONL, collect paths from `Write` / `Edit` / `NotebookEdit` tool inputs and best-effort regex parsing of Bash mutations (`mv`/`cp`/`rm`/`tee`/`sed -i`/`>`/`>>`).
5. **`git add <paths>`**: check each path against the session set. Block (exit 2) if any are out-of-session. Render the offending paths and the session set in the stderr message.
6. **`git commit`**: query `git diff --cached --name-only` for currently-staged paths, check each. Block with a `git reset HEAD <paths>` remediation if any are out-of-session.
7. **Uncertain session set** (empty / no transcript): fail OPEN for specific-path operations. Bulk-stage refusal still applies.

### Sync rule + settings

- `sync.config.yaml`: added `hooks/` under `sync_rules.always` so `claude-config/hooks/` syncs to `~/.claude/hooks/`.
- `settings.json`: added `hooks.PreToolUse` with a single entry — `matcher: Bash`, `command: $HOME/.claude/hooks/session-scope-guard.py`. Uses `$HOME` so the path is stable across the dev's filesystem.

## Expected Outcome

- Triggering: the skill remains discoverable (description unchanged in triggering surface), but invocation now lands on a brief, scannable body instead of 400 lines of machinery.
- Safety: bulk-stage and out-of-session commits are structurally rejected by the hook regardless of whether the skill is invoked — the hook fires on every Bash tool call, the skill only fires when the model deems it useful.
- Reversibility: bypass via `SKIP_SESSION_SCOPE_GUARD` in the Bash description for intentional cases, or `SESSION_SCOPE_BYPASS=1` env var for full session bypass.

## Actual Outcome

All seven smoke tests pass:

| # | Scenario | Expected | Actual |
|---|---|---|---|
| 1 | Non-Bash tool | exit 0 | exit 0 |
| 2 | Non-git Bash | exit 0 | exit 0 |
| 3 | `git add -A` | exit 2 (block) | exit 2, clear stderr |
| 4 | `git add -A` + `SKIP_SESSION_SCOPE_GUARD` description | exit 0 | exit 0 |
| 5 | `git add . && git commit` | exit 2 (block on `.`) | exit 2 |
| 6 | `git commit` with no transcript (uncertain) | exit 0 (fail open) | exit 0 |
| 7 | `git add <session-path>` vs `git add <out-of-session>` | allow / block respectively | matches |

Live state verified: `~/.claude/hooks/session-scope-guard.py` exists, is executable (chmod +x), and runs cleanly when invoked with a synthetic stdin payload. Live `~/.claude/settings.json` includes the `PreToolUse` hook config under the `Bash` matcher. Live `~/.claude/skills/git-strategy-advisor/SKILL.md` is 99 lines.

## Assessment

**Result**: Success

**Improvements**:
- Skill body is 4× smaller (405 → 99 lines). Reading the skill no longer floods the model with orchestrator-style scaffolding it doesn't need.
- Session-scope discipline is now enforced by the hook regardless of whether the skill fires. The advisory layer ("here's a recommendation") and the structural layer ("this commit is rejected") are properly separated.
- The hook is reusable infrastructure — any future skill that mutates git can rely on the same guard without re-implementing it.
- Bulk-stage patterns are now structurally unreachable for Claude unless the user explicitly bypasses.

**Issues**:
- The session-set derivation depends on transcript parsing, which is best-effort. If Claude Code's transcript format changes, the regex/JSON walk may need updates. Tested against the current `message.content[].tool_use` shape (both with and without the outer `message` wrapper).
- The hook fails OPEN when uncertain (no transcript, malformed JSON, no session set). Trade-off is "occasionally allow a bad commit" rather than "block legitimate work" — appropriate for advisory infrastructure.
- Bypass mechanisms (`SESSION_SCOPE_BYPASS=1` and `SKIP_SESSION_SCOPE_GUARD`) put the user in control but could be misused. Documented in the SKILL.md "Notes" section: do not bypass silently.

**Lessons Learned**:
- The right separation: lightweight advisory in the skill, structural enforcement in a hook. Over-fitting a skill body with rules the model is supposed to "follow" gets you neither — the skill bloats and the rules are still discretionary.
- For session-scope, the transcript is the authoritative source of "what did Claude touch this session" — beats both `git status` (too coarse, includes pre-session dirt) and explicit caller declaration (requires a contract every caller honors).
- Smoke-testing a hook with synthetic stdin payloads is fast and catches the obvious failures before any real workflow hits it.

## Related Commits

- `dee3a86`: chore(git-strategy-advisor): slim to advisory layer + add session-scope-guard hook

## Next Steps

- If `weekly-digest`, `meeting-notes-to-action`, or other future skills mutate git state, no additional work is needed — the hook covers them automatically.
- Consider a small `--diagnose` mode on the hook that prints what it would block (without actually blocking) given a stdin payload, useful for debugging when the user is surprised by a refusal.
- If transcript format changes (Claude Code version bump), update the `_collect_from_event` walker in the hook.
