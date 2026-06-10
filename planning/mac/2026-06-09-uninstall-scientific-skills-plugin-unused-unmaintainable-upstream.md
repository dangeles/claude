# Uninstall scientific-skills plugin (unused, unmaintainable upstream)

**Date**: 2026-06-09
**Machine**: mac
**Status**: Implemented

## Objective

Remove the `scientific-skills@claude-scientific-skills` plugin (146 skills, K-Dense-AI).

## Rationale

- **Unused**: across all transcripts, exactly 2 invocations ever (`scientific-writing`,
  `research-lookup`) — both with near-equivalents in the local curated skill set.
- **Standing context cost**: 146 skill descriptions injected into every session's menu
  (~6-9k tokens) for ~zero use.
- **Unmaintainable**: upstream removed `marketplace.json` (repo commit 398fc37; was v2.34.2)
  and migrated to the `gh skill` / agentskills.io distribution model, so the plugin could no
  longer be updated via `claude plugin update` (resolves "Plugin not found"). Installed copy
  was frozen at the 2026-01-27 snapshot (sha a31cf4dd).

## Changes

- `claude plugin uninstall scientific-skills@claude-scientific-skills` (user scope).
- `claude plugin marketplace remove claude-scientific-skills` (now-orphaned marketplace).
- Captured the resulting live state into the repo (copy live→repo, NOT push — a push would
  have re-added it): `claude-config/settings.json` (dropped enabledPlugins entry) and
  `claude-config/plugins/installed_plugins.json` (13→12 plugins).

## Outcome

`sync-config.py status` → "No changes detected" (repo == live). 12 plugins remain, all
confirmed at latest version in the same pass.

**Result**: Success

## Re-adding later

The old plugin path is gone upstream. To get the current skills, use the new model:
`gh skill install K-Dense-AI/scientific-agent-skills --agent claude-code` (gh CLI ≥ 2.90.0).

## Related Commits

- [pending]: chore(plugins): uninstall scientific-skills plugin
