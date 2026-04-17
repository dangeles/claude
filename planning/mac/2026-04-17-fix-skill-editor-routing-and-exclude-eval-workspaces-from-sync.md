# fix skill-editor routing and exclude eval workspaces from sync

**Date**: 2026-04-17
**Machine**: mac
**Status**: Success

## Objective

Two related fixes:

1. **skill-editor routing**: A single line in skill-editor/SKILL.md explicitly redirects new-skill
   creation to skill-creator ("Use skill-creator for new skills from scratch without quality gates").
   This bypasses the repo protocol (CONFIG_MANAGEMENT.md, planning journal, sync-config.py) entirely.
   Fix: remove the redirect, tighten the description to be the default entry point for all repo
   skill work.

2. **Workspace exclusion**: skill-creator eval workspaces land in `~/.claude/skills/*-workspace/`
   as siblings of skill dirs. sync.config.yaml syncs the entire `skills/` path, so running
   `sync-config.py pull` would drag ephemeral eval artifacts into the repository. Fix: add
   `skills/*-workspace/` to sync exclusions.

## Changes Planned

- [ ] Edit `claude-config/skills/skill-editor/SKILL.md`:
  - Frontmatter description: rewrite to third-person, remove subjective qualifier ("require"),
    cover all repo skill work (new + modify + refactor)
  - "When NOT to Use" section: remove the skill-creator redirect for new skills from scratch
  - Add one sentence about skill-creator's eval pipeline being optionally available in Phase 4
- [ ] Edit `sync.config.yaml`: add `skills/*-workspace/` to exclusions

## Expected Outcome

- Claude routes to skill-editor for ALL skill work in this repo, including new skills from scratch
- skill-creator is understood as the eval/benchmark layer, optionally invoked from within
  skill-editor Phase 4 — not as a competing standalone workflow
- Eval workspaces never appear as divergences in `sync-config.py status`

## Actual Outcome

- skill-editor description rewritten; trigger phrases now explicit
- skill-creator redirect removed from "When NOT to Use" section
- sync.config.yaml exclusion added; sync-config.py updated with fnmatch glob
  support (previously only exact prefix matching existed)
- `sync-config.py status` now clean of workspace orphans

## Assessment

**Result**: Success

**Improvements**:
- skill-editor is now the clear entry point for all repo skill work
- skill-creator and plugin-dev:skill-development have defined supporting roles
- Eval workspaces no longer pollute sync status

**Issues**:
- None

**Lessons Learned**:
- sync.config.yaml exclusion patterns didn't support globs — needed a code
  change alongside the config change
- Exclusion paths are relative to the *rule's subdirectory*, not the repo
  root (e.g., `*-workspace` not `skills/*-workspace`)

## Related Commits

- 69db61f: fix(skill-editor): clarify routing and exclude eval workspaces from sync

## Next Steps

- Consider whether skill-editor Phase 4 should offer an optional eval step via skill-creator
  infrastructure in a future iteration
