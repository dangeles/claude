# Project Organization Report

**Project**: `claude` (Claude Code config repo)
**Session**: archive-workflow-session-20260524-160609-57202
**Date**: 2026-05-24
**Type**: Cleanup pass (incremental — baseline from 2026-02-21)

## Summary

This pass reconciled documentation drift, normalized a handful of uppercase filenames to kebab-case, removed one deprecated skill, and refreshed `.archive-metadata.yaml` to reflect the current structure (including the `claude-config/hooks/` directory added since the February baseline). Scope was deliberately narrow — no large structural moves.

**Before**: 55 skills, `project-configs/` referenced inconsistently (sync.config.yaml said `claude-config/project-configs`, docs said `project-configs/` at repo root), 2 uppercase test filenames + 1 uppercase bibliography file + 1 version-suffixed schema filename in violation of project naming conventions, deprecated `software-developer` tombstone still present.

**After**: 54 skills, `claude-config/project-configs/` consistent everywhere, all renamed to kebab-case (with schema version moved from filename to JSON body), tombstone removed, metadata regenerated.

---

## Changes Made

### Renames (4 files via `git mv`)
| Original | New |
|---|---|
| `claude-config/skills/programming-pm/references/testing/END-TO-END-TEST-PLAN.md` | `.../end-to-end-test-plan.md` |
| `claude-config/skills/programming-pm/references/testing/TEST-RESULTS.md` | `.../test-results.md` |
| `claude-config/skills/plotting-advisor/references/SOURCES.md` | `.../sources.md` |
| `claude-config/skills/workflow-coordinator/references/universal-handoff-schema-v3.0.json` | `.../universal-handoff-schema.json` |

Internal references updated across `implementation-summary.md`, `plotting-advisor/SKILL.md`, `workflow-coordinator/SKILL.md` + 4 reference files + 1 example JSON.

### Documentation reconciliation (C1)
`project-configs/` references in 4 doc files (9 lines total) updated to `claude-config/project-configs/` to match the actual location used by `sync.config.yaml`. Empty `project-configs/` directory at repo root removed via `rmdir`.

### Deletions (C5)
- Removed `claude-config/skills/software-developer/` (1 file: deprecated SKILL.md tombstone). Pre-flight grep confirmed no live orchestrator dispatches to it; remaining references (CHANGELOG entries, planning files) are historical and intentionally preserved for archaeology.

### Schema versioning convention
Per user override, `universal-handoff-schema-v3.0.json` → `universal-handoff-schema.json`. The `"version": "3.0"` constant inside the JSON body was already present (lines 13–17), so it now serves as the sole source of truth for schema version.

### Metadata refresh
`.archive-metadata.yaml` regenerated. Notable additions:
- `claude-config/hooks/` entry (the hooks directory was added between the Feb baseline and now)
- Naming rule for `claude-config/hooks/*.py` (kebab-case entrypoints, `_hook_lib.py` snake_case shared)
- Naming rule for `*.schema.json` (version belongs in JSON body, not filename)
- `project.skill_count: 54`
- `history[]` array tracking both archive-workflow runs

---

## Decisions Made (Phase 1 conflict resolution)

| Item | Decision | Source |
|---|---|---|
| C1 (project-configs location) | Option 1: nest under claude-config/ (docs follow code) | User approval |
| C2 (SOURCES.md case) | Rename to sources.md (kebab-case wins over emphasis caps) | User override |
| C3 (schema version) | Move v3.0 from filename to JSON body field | User override |
| C4 (local-disk orphan) | Skip — out of scope | User instruction |
| C5 (software-developer tombstone) | Remove (no live references) | User override, pre-flight verified |

---

## Files Touched (full list)

**Renamed** (4):
- `claude-config/skills/programming-pm/references/testing/end-to-end-test-plan.md`
- `claude-config/skills/programming-pm/references/testing/test-results.md`
- `claude-config/skills/plotting-advisor/references/sources.md`
- `claude-config/skills/workflow-coordinator/references/universal-handoff-schema.json`

**Modified** (11):
- `CLAUDE.md`
- `README.md`
- `.archive-metadata.yaml`
- `docs/CLAUDE_CONFIG_GUIDE.md`
- `docs/SYNC_WORKFLOW.md`
- `claude-config/skills/plotting-advisor/SKILL.md`
- `claude-config/skills/programming-pm/references/implementation-summary.md`
- `claude-config/skills/workflow-coordinator/SKILL.md`
- `claude-config/skills/workflow-coordinator/examples/skill-editor-to-programming-pm.json`
- `claude-config/skills/workflow-coordinator/references/distributed-tracing.md`
- `claude-config/skills/workflow-coordinator/references/handoff-registry.md`
- `claude-config/skills/workflow-coordinator/references/handoff-validation.md`

**Deleted** (1):
- `claude-config/skills/software-developer/SKILL.md`

**Directory removed** (1):
- `project-configs/` (was empty at repo root)

---

## Before / After Counts

| Metric | Before | After |
|---|---|---|
| Skills (SKILL.md count) | 55 | 54 |
| Uppercase filename violations | 4 | 0 |
| Version-suffixed schema files | 1 | 0 |
| Doc references to `project-configs/` (orphan path) | 9 | 0 |
| Top-level dirs in `.archive-metadata.yaml` | 3 | 7 (added hooks/, skills/, agents/, project-configs/ subdirs) |

---

## Next Steps

1. Run `./sync-config.py push --dry-run` to verify the live `~/.claude/skills/` will receive the rename + deletion correctly.
2. After push, run a fresh session and sanity-check that `workflow-coordinator` still loads handoff schema via the new filename.
3. Optional follow-up: review whether the `principal-investigator/references/CHANGELOG.md` historical references to `software-developer` should be marked-up (e.g., `[deprecated 2026-05-12]` annotations) for future readers. Not done in this pass to keep the diff narrow.

## Rollback Instructions

If you need to undo this pass:
```bash
git reset --hard 14ef56e9cd14adb99bbc3cf93c9fda3ca89b3622
./sync-config.py push  # re-push the prior state to ~/.claude/
```

## Session Files

Full session logs and reports:
- `/tmp/archive-workflow-session-20260524-160609-57202/execution-log.md`
- `/tmp/archive-workflow-session-20260524-160609-57202/clutter-report.md`
- `/tmp/archive-workflow-session-20260524-160609-57202/naming-violations.md`
- `/tmp/archive-workflow-session-20260524-160609-57202/structure-proposal.md`
- `/tmp/archive-workflow-session-20260524-160609-57202/expandability-assessment.md`
