# Pull superpowers 6.2.0 plugin update

**Date**: 2026-07-26
**Machine**: mac
**Status**: Success

## Objective

The live system had drifted ahead of the repo: Claude Code's plugin manager
auto-updated `superpowers` from 5.1.0 to 6.2.0 on 2026-07-24, and refreshed the
`lastUpdated` timestamps on eight other official plugins. This is the inverted
case for this repo — `~/.claude/` held the newer truth, so the version pin in
`claude-config/plugins/installed_plugins.json` needed to catch up rather than be
pushed back down.

Pushing instead of pulling would have tried to force the live install back to
superpowers 5.1.0 and fought the plugin manager.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Confirm the divergence is live-side only, not repo drift
- [x] `pull --dry-run` to confirm blast radius before applying
- [x] `pull` to capture the new version pin
- [x] Re-run `validate` and the test suite

## Expected Outcome

Repo pin matches the live install, `sync-config.py status` reports no
divergence, and no other tracked config is disturbed.

## Actual Outcome

One file changed: `claude-config/plugins/installed_plugins.json`, 11 insertions
and 11 deletions — the superpowers `version` and `installPath` bump to 6.2.0,
plus nine refreshed `lastUpdated` timestamps.

Notably, `skills/`, `agents/`, and `hooks/` were untouched. The superpowers 6.x
skill content lives under `plugins/cache/`, which `sync.config.yaml` excludes, so
a major-version plugin bump does not sweep skill bodies into the repo. Worth
knowing: the tracked pin records *which* version is installed, not its contents.

`validate` passed (frontmatter + JSON gates), the 67-test suite still passes, and
`status` no longer reports divergence.

## Assessment

**Result**: Success

**Improvements**:
- Repo and live system are back in agreement; the session-start divergence
  warning is cleared.
- The 5.1.0 → 6.2.0 jump is now recorded in git history, so a fresh machine
  provisioned from this repo gets 6.2.0 instead of silently reinstalling 5.1.0.

**Issues**:
- `pull` cannot run non-interactively without `--yes`: the conflict resolver
  calls `input()` and dies with `EOFError` on a closed stdin, even under
  `--dry-run`. `--yes` resolves it correctly here (source wins = live wins on a
  pull), but the failure mode is a traceback rather than a clear message.
- `--dry-run` still writes a real backup file into `claude-config/.backups/`.
  Harmless (the directory is gitignored) but it is a side effect in a mode that
  claims to have none.
- Plugin `lastUpdated` timestamps churn on essentially every Claude Code launch,
  so `status` will keep reporting drift on this file for reasons unrelated to
  any real version change. Signal-to-noise on that warning is poor.

**Lessons Learned**:
- Read a divergence warning's *direction* before reacting. The reflex in this
  repo is "repo is source of truth, therefore push," and here that reflex would
  have been actively wrong.
- Timestamp-only churn in `installed_plugins.json` is not a change worth
  chasing; check whether a `version` field actually moved before doing anything.

## Related Commits

- a3903f9: chore(plugins): pull superpowers 6.2.0 version pin

## Next Steps

Consider stripping or normalizing `lastUpdated` when syncing
`installed_plugins.json`, so that `status` flags only genuine version changes.
That would make the session-start divergence warning trustworthy instead of
routine background noise.
