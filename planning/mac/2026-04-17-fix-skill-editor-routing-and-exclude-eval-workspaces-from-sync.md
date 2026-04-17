# fix skill-editor routing and exclude eval workspaces from sync

**Date**: 2026-04-17
**Machine**: mac
**Status**: In Progress

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

[After implementation: What actually happened?]

## Assessment

**Result**: [Success / Partial / Failed]

**Improvements**:
- [What got better?]

**Issues**:
- [What problems emerged?]

**Lessons Learned**:
- [What would you do differently?]

## Related Commits

- [pending]

## Next Steps

- Consider whether skill-editor Phase 4 should offer an optional eval step via skill-creator
  infrastructure in a future iteration
