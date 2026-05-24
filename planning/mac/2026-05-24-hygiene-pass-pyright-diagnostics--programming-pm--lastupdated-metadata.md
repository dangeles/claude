# Hygiene pass: pyright diagnostics + programming-pm + last_updated metadata

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

Three follow-ups carried over from the Tier-1+2 audit and the hostname-fix work:

1. **Pyright diagnostics** — 6 pre-existing static-analysis warnings in `sync-config.py` (unused imports, unused parameters, structurally unreachable code).
2. **programming-pm Team Composition cleanup** — inline section duplicated content already in `references/team-composition.md`; the inline copy was 46 lines when it could be a 7-line cheat-sheet pointing at the reference.
3. **Metadata hygiene** — 41 of 57 skills were missing the `last_updated` field, leaving the catalog inconsistent.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Fix sync-config.py pyright diagnostics: `Dict` unused → remove from imports; `label1/label2` unused → drop the dead parameters from `show_diff`; `pattern` unused → drop the assignment; `abs_path` unused → use plain `_`; `sys.platform == 'darwin'` unreachable-branch warning → use a module-level `_IS_DARWIN: bool` constant so pyright can't narrow
- [x] Slim programming-pm Team Composition inline section (46 lines → 7 lines: heading + Default Team table + one-line pointer to the full reference)
- [x] Add `last_updated: 2026-05-24` to all 41 skills missing it
- [x] Validate frontmatter for all 55 SKILL.md files
- [x] `./sync-config.py push --yes` and verify live state

## Expected Outcome

- Pyright diagnostics down from 6 → 0 in `sync-config.py`.
- programming-pm SKILL.md drops another ~40 lines while preserving the at-a-glance Default Team table.
- All 55 skills carry consistent `last_updated` metadata, making it easy to see catalog freshness at a glance.

## Actual Outcome

All three changes shipped cleanly.

**Pyright fixes** (`sync-config.py`):
- Removed unused `Dict` from `from typing import ...`
- Removed unused `label1: str = "source", label2: str = "target"` parameters from `show_diff` (no callers pass them — verified by grep)
- Removed unused `pattern` assignment in the wildcard-project-discovery block
- Renamed unused `abs_path` loop variable to `_` (Pyright recognizes the bare underscore but not the underscore-prefix convention)
- Introduced module-level `_IS_DARWIN: bool = sys.platform == 'darwin'` and replaced the inline platform check with a reference to it. Pyright narrows `sys.platform` to a literal on the host running the type checker (treating one branch as unreachable); routing through a `bool` constant defeats the narrowing while keeping both runtime branches alive.

**programming-pm Team Composition**:
- Replaced 46-line inline section with a 7-line cheat-sheet: heading, Default Team table (5 specialists × 3 columns), and one-line pointer listing what `references/team-composition.md` covers (optional specialists with inclusion criteria, decision tree, RACI matrix, team-size guidelines, override flags).
- SKILL.md drops another ~40 lines (1905 → ~1865).

**`last_updated` metadata**:
- Python script identified 41 SKILL.md files missing the field, inserted `last_updated: 2026-05-24` immediately after the last of `name:` / `version:` (preserves the standard YAML ordering used by sibling skills like calculator and notebook-debugger).
- Final scan: 55/55 SKILL.md files have valid frontmatter; all 55 now carry `last_updated`.
- `sync-config.py push --yes` reported 83 successful copy/overwrite operations (consistent with the 41 file updates × 2 log lines each).

## Assessment

**Result**: Success

**Improvements**:
- `sync-config.py` has zero pyright diagnostics (confirmed via re-run after final edit).
- The `_IS_DARWIN` constant pattern is reusable — any future platform-specific branch can follow it to avoid the narrowing artifact.
- programming-pm SKILL.md is now closer to a pure orchestrator body (workflow logic) with the cheat-sheet only where it earns its weight inline.
- The catalog has a single observable "freshness" field per skill — useful for future audits ("which skills haven't been touched in 90 days?").

**Issues**:
- None observed.

**Lessons Learned**:
- Pyright's `_var` underscore-prefix convention isn't a recognized "unused-OK" signal; only plain `_` is. Worth knowing for future tuple unpacks.
- Pyright narrows `sys.platform` to a literal on the dev host — adding `# pyright: ignore[reportUnreachable]` or `# type: ignore` on the apparently-unreachable line didn't suppress the warning; the cleanest fix is to route through a typed variable.
- A bulk metadata pass can be safely done with a Python script that uses `yaml.safe_load` to validate every file after the edit. 41 files updated in one pass with zero validation failures.

## Related Commits

- `<pending>`: chore: pyright cleanup + programming-pm cheat-sheet + last_updated metadata pass

## Next Steps

None outstanding from the prior carry-over list. Future work surfaced by this session:

- Consider extending `_IS_DARWIN` to a small `_PLATFORM` enum if more platform branches accumulate in `sync-config.py`.
- The `last_updated` field is now a clean signal for future audits — could power a `./sync-config.py audit --stale-days 90` style report later.
