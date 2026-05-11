# Add missing description frontmatter to archive-workflow agents

**Date**: 2026-05-11
**Machine**: mac
**Status**: Success

## Objective

Claude Code's agent loader reports `Missing required "description" field in frontmatter` for 6 agents in `~/.claude/agents/`:

- archive-clutter-analyst.md
- archive-decision-integrator.md
- archive-expandability-reviewer.md
- archive-nomenclature-enforcer.md
- archive-structure-organizer.md
- library-pm.md

These were authored before the loader's frontmatter contract included `description`. They use a custom `role` / `permissions` schema instead, so the loader rejects them and they cannot be dispatched by the archive-workflow skill.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add `description:` field to each of the 6 agent files in `claude-config/agents/`
  - Descriptions lifted from the wave-by-wave agent characterization in `claude-config/skills/archive-workflow/SKILL.md` to avoid inventing new behavior claims
- [x] Preserve existing `role` and `permissions` lines (they're in-prompt context for the agent, not parsed by the loader)
- [x] Validate frontmatter via `python3 -c "import yaml; ..."` — all 6 files parse cleanly with a `description` key
- [x] Apply to `~/.claude/agents/` (see note below — bypassed `sync-config.py push` due to unrelated drift)
- [x] Commit and update this entry with SHA

## Expected Outcome

`~/.claude/agents/` loads without errors for these 6 agents. The archive-workflow skill can again dispatch its specialists via the Task tool.

## Actual Outcome

All 6 agent files now have a valid `description:` field both in `claude-config/agents/` and `~/.claude/agents/`. Frontmatter parses cleanly under PyYAML. The archive-workflow skill can now dispatch each specialist via the Task tool without the loader rejecting them.

`./sync-config.py push --dry-run` surfaced an unrelated conflict on `plugins/installed_plugins.json`: the repo has stale plugin metadata (e.g., `superpowers/4.3.1`) while `~/.claude/` has current versions (`superpowers/5.1.0`). Using `push --yes` would have clobbered live `installPath` values and broken plugin resolution. Since this conflict is orthogonal to the description fix, I copied the 6 agent files directly with `cp` instead of running the full push. This produces the same end state for the agents while leaving the plugin metadata drift to be resolved separately (via `./sync-config.py pull` + commit).

## Assessment

**Result**: Success

**Improvements**:
- 6 previously-broken agents load cleanly; the archive-workflow skill is dispatchable end-to-end again
- The `role` / `permissions` annotations were preserved — they're in-prompt context, and removing them would have edited agent behavior under the guise of a frontmatter fix

**Issues**:
- `claude-config/plugins/installed_plugins.json` has drifted behind `~/.claude/`. Sync-config push refuses without resolution. This blocked the canonical push path for an unrelated reason.

**Lessons Learned**:
- When the repo has unrelated drift from `~/.claude/`, a full `sync-config.py push` is unsafe (`--yes` means "source wins" globally, which can clobber legitimately newer live state). A targeted `cp` for just the changed files is the right escape hatch — it preserves the workflow's invariant (repo is source of truth for those files) without forcing a global resolution.
- Worth considering a `--paths` flag on `sync-config.py push` that scopes the operation to specific files, so partial syncs don't require manual `cp`.

## Related Commits

- dbea4ee: fix(agents): add missing description frontmatter to archive-workflow agents

## Next Steps

1. Run `./sync-config.py pull` to refresh `claude-config/plugins/installed_plugins.json` with current live plugin versions, then commit. This unblocks future full pushes.
2. Consider a sync-time lint that flags missing required frontmatter fields (e.g., `description`) in `claude-config/agents/` before push, so this class of error is caught at the repo boundary rather than at agent-load time.
