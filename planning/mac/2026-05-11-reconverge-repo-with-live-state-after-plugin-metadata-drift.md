# Reconverge repo with live state after plugin metadata drift

**Date**: 2026-05-11
**Machine**: mac
**Status**: Success

## Objective

Earlier today, `./sync-config.py push --dry-run` (while fixing the archive-workflow agent descriptions) surfaced unrelated drift between `~/.claude/` and `claude-config/`:

1. `plugins/installed_plugins.json` — live had newer plugin versions (e.g., `superpowers/5.1.0` vs repo's `4.3.1`), some plugins shifted to `version: "unknown"`. Last-updated timestamps were newer in live.
2. A case-conflict on `skills/notebook-writer/SKILL.md` (repo, uppercase) vs `skill.md` (live, lowercase). Same content (identical SHA1), just naming-mismatched. APFS is case-insensitive on macOS, so the loader works fine, but `diff` and any future sync sees them as different files.
3. Twelve empty placeholder subdirs (`references/`, `assets/`, `scripts/`) present only in live; three empty placeholder subdirs present only in repo. None of them git-tracked (empty dirs).

Push was blocked because `--yes` (source/repo wins) would have clobbered live `installPath` values and broken plugin resolution. Workaround at the time: `cp` the 6 fixed agent files to live directly. This entry closes the loop by reconverging properly.

Follow-up to `2026-05-11-add-frontmatter-validation-to-sync-config-push-gate.md` — the new linter is now in place to validate the reconverged state.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Resolve case-conflict first: rename live `skill.md` → `SKILL.md` (via two-step `mv` through a `.tmp` intermediate because APFS is case-insensitive). Convention across every other skill is uppercase.
- [x] Confirm `./sync-config.py status` no longer reports the skills diff after rename
- [x] Run `./sync-config.py pull --yes` (no `--delete`) to bring live's plugin metadata into the repo. `--yes` means source/live wins, which is what we want for this one-way reconverge; omitting `--delete` preserves the 3 repo-only placeholder subdirs in case they were authored content
- [x] Run `./sync-config.py validate` against the pulled state to confirm the new linter passes (would catch any pulled file that lacks `name`/`description` frontmatter)
- [x] Review `git diff` on `installed_plugins.json` to ensure the changes are bounded plugin metadata, not anything more invasive
- [x] Commit

## Expected Outcome

- `./sync-config.py status` reports "No changes detected"
- `./sync-config.py validate` passes (all 82 agent + SKILL.md files have required frontmatter)
- Repo plugin metadata matches live exactly
- `diff -rq` still shows some empty-subdir asymmetry, but this does not affect any sync-tracked state — empty dirs aren't git-tracked and aren't compared by sync-config's file-hash logic

## Actual Outcome

Matches expectations:

| Check | Expected | Observed |
|---|---|---|
| `./sync-config.py status` after pull | "No changes detected" | "No changes detected" |
| `./sync-config.py validate` after pull | success, exit 0 | success, exit 0 |
| Files modified in git | only `installed_plugins.json` (54 lines, 27/27 insertions/deletions) | confirmed via `git diff --stat` |
| `installed_plugins.json` diff content | bounded plugin metadata (version strings, installPath, lastUpdated) | confirmed; no surprises |

Key data points from the diff:

- `superpowers`: 4.3.1 → 5.1.0 — a real upgrade
- `frontend-design`, `commit-commands`, `pr-review-toolkit`, `code-review`, `feature-dev`: `version` and `installPath` shifted from a specific commit-sha string (`55b58ec6e564`) to `"unknown"`. All five share the same `lastUpdated` timestamp of `2026-04-17T22:57:27`, suggesting a single sync/refresh event reset them.

The `"unknown"` version strings are an information regression (specific sha → opaque label), but they reflect the truth of what Claude Code is actually using on disk. The repo should match reality, not retain stale-but-prettier metadata.

## Assessment

**Result**: Success

**Improvements**:
- `./sync-config.py push` is now unblocked. Future fixes that touch the global config can flow through the canonical sync workflow without the `cp` workaround.
- The frontmatter linter added in the previous entry was exercised by this reconverge and confirmed working — it validated the pulled state cleanly.
- `notebook-writer/SKILL.md` casing is now uniform with every other skill, removing a future foot-gun on case-sensitive filesystems (Linux, CI).

**Issues**:
- The five plugins with `"version": "unknown"` is a known curiosity — would be worth a separate investigation if it matters (does the marketplace re-resolve them? does it block updates?), but is not blocking for sync purposes.
- Empty placeholder subdirs (12 live-only, 3 repo-only) remain asymmetric. Sync ignores them, git ignores them, but `diff -rq` shows them. Could be unified later by adding `.gitkeep` markers in canonical locations, or by deleting all empty placeholders. Not urgent.

**Lessons Learned**:
- When `push` is blocked by unrelated drift, the right move is often to `pull` the drift first to reconverge, *then* `push` whatever original fix triggered the discovery. Trying to surgical-merge in both directions simultaneously (the `cp` workaround used earlier today) is correct as an emergency move but accumulates entropy.
- Case-conflicts on case-insensitive filesystems need a two-step `mv` via a temp name (`mv skill.md SKILL.md.tmp && mv SKILL.md.tmp SKILL.md`) — a direct `mv skill.md SKILL.md` is a no-op on APFS.
- `diff -rq` and `sync-config.py status` disagree about what counts as drift. The latter is the authoritative signal for sync purposes (it compares file hashes, skips empty dirs); `diff -rq` is more useful for hunting filesystem-level inconsistencies.

## Related Commits

- 0020ad4: chore(sync): reconverge repo with live state; record validate SHA

## Next Steps

1. Optional cleanup: unify the empty placeholder subdirs (either add `.gitkeep` to all expected ones, or `rmdir` the inconsistent ones). Low priority.
2. Investigate why five plugins show `version: "unknown"` — is that a marketplace registration issue or a Claude Code state quirk? Not blocking.
3. Push remaining commits to origin.
