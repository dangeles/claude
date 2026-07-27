---
name: software-developer
last_updated: 2026-05-24
description: Deprecated as of 2026-05-12. Use python-developer for production-quality Python implementation; chain through bioinformatician for bioinformatics-specific analysis and biologist-commentator for biological validation.
deprecated: true
deprecated_date: 2026-05-12
replaced_by: python-developer
---

# software-developer (deprecated)

This skill was deprecated on 2026-05-12 as a duplicate of `python-developer`. Both skills provided production-quality Python implementation; the only meaningful difference was bioinformatics framing in this skill's prose examples. The actual Python practices (type hints, docstrings, pytest, error handling, CLI scaffolding) are identical.

## What to use instead

| Need | Use |
|---|---|
| Production-quality Python implementation | `python-developer` |
| Bioinformatics-specific data analysis | `bioinformatician` |
| Biological validation of bioinformatics code | `biologist-commentator` |
| Bioinformatics production code | `bioinformatician` → `python-developer` → `biologist-commentator` |

`python-developer` integrates with `mathematician`, `statistician`, `copilot`, and `programming-pm` for full pipeline coordination.

## Rationale

The closer-reading audit (2026-05-12, planning entry `2026-05-12-tier-2-deprecate-software-developer-skill.md`) found:
- Identical responsibilities (production code with tests, docs, type hints, error handling)
- Identical tooling expectations (pytest, type checking, coverage)
- `software-developer` had bioinformatics framing in examples; `python-developer` had broader integration (mathematician/statistician handoffs)
- Bioinformatics-specific work is better served by chaining `bioinformatician` (analysis) → `python-developer` (productionize) → `biologist-commentator` (validate) than by maintaining a separate "production bioinformatics" skill

## Migration

Orchestrators that previously dispatched `software-developer` (programming-pm, principal-investigator, biologist-commentator, systematic-troubleshooter, systems-architect, copilot) now reference `python-developer` directly. See the deprecation planning entry for the full file list updated in this change.
