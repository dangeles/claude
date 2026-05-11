# Add frontmatter validation to sync-config push gate

**Date**: 2026-05-11
**Machine**: mac
**Status**: Success

## Objective

Earlier today, 6 archive-workflow agents shipped to `~/.claude/agents/` without a `description:` frontmatter field and were silently rejected by the Claude Code agent loader. The breakage was only discovered at agent-load time — by which point the bad state was already live.

This entry closes that loop at the repo boundary: add a frontmatter linter to `sync-config.py` that runs automatically before every `push` (and is also exposed as a standalone `validate` subcommand), so files missing required loader keys can never reach `~/.claude/` via the canonical sync path.

Follow-up to the entry at `2026-05-11-add-missing-description-frontmatter-to-archive-workflow-agents.md`.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add `_parse_frontmatter()` helper to `ConfigSync` (extracts and parses the YAML block; returns `(dict, error_message)` tuple)
- [x] Add `validate_frontmatter()` method that scans `claude-config/agents/*.md` and `claude-config/skills/*/SKILL.md`, requires `name` + `description` keys, prints per-file errors, returns bool
- [x] Wire `validate_frontmatter()` as the first step of `push_config()` — aborts with `sys.exit(1)` on failure, before any user prompt or file modification
- [x] Add `validate` to the `argparse` `command` choices and main() dispatch
- [x] Update `--help` epilog with the new subcommand
- [x] Module-level constant `REQUIRED_FRONTMATTER_FIELDS = ('name', 'description')` so the contract is one-line to amend if Claude Code adds required keys later

## Expected Outcome

- `./sync-config.py validate` runs in <1s on the current tree, exits 0, prints a success line
- `./sync-config.py push` (and `push --dry-run`) abort with exit 1 if any agent or SKILL.md has missing/empty `name` or `description`, or unparseable YAML frontmatter
- The push prompt is never shown when validation fails — the user doesn't get to confirm a doomed operation

## Actual Outcome

Matches expectations. Verified manually:

| Scenario | Expected | Observed |
|---|---|---|
| `validate` on current (just-fixed) tree | exit 0, success line | exit 0, success line |
| `validate` after deleting `description:` from `library-pm.md` | exit 1, per-file error | exit 1, error names the file and the missing field |
| `push --dry-run` with same broken frontmatter | abort before any sync, exit 1 | abort before any sync, exit 1 |
| `validate` after restoring `library-pm.md` | exit 0 | exit 0 |

The linter scanned 82 files (agents + SKILL.md) — fast enough that running it on every push is essentially free.

## Assessment

**Result**: Success

**Improvements**:
- The agent-loader rejection class of error is now caught at push time, not at load time. The repo→live propagation path is gated.
- Adding `validate` as a first-class subcommand also makes the check usable from a git pre-commit hook later, without further code changes.
- Module-level `REQUIRED_FRONTMATTER_FIELDS` keeps the contract narrow and easy to extend if the loader adds requirements (e.g., `tools`, `model`).

**Issues**:
- Pyright initially flagged a type-narrowing issue on the `(data, err)` tuple return from `_parse_frontmatter` because Pyright can't infer "if err is None then data is not None" from the convention. Added an explicit `data is None` guard to narrow.
- Validation runs *before* `--dry-run` short-circuits the actual copy, which is intentional but worth flagging: dry-run will now fail loudly if frontmatter is broken, where previously it would have appeared to "succeed" with garbage about to be synced.

**Lessons Learned**:
- Keep the validator's scope tight to what the loader actually checks. Validating other plausible-but-unrequired fields (e.g., model names) creates false positives the first time the loader's contract drifts.
- Wiring validators to abort *before* the confirmation prompt is the right ordering. Aborting after a "Continue? [y/N]" is user-hostile.

## Related Commits

- 7fe3357: feat(sync-config): add frontmatter validation as push-time gate

## Next Steps

1. Reconverge the repo with live state (`installed_plugins.json` + skill subdir drift). The linter is in place to validate the result.
2. Consider a git pre-commit hook that runs `./sync-config.py validate` — catches the error one step earlier, before commit. Trade-off is one more moving part to maintain across machines.
3. If the loader adds new required fields in the future, amend `REQUIRED_FRONTMATTER_FIELDS` in `sync-config.py` (one-line change).
