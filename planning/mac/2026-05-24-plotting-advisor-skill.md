# Add plotting-advisor skill

**Date**: 2026-05-24
**Machine**: mac
**Status**: In Progress

## Objective

Add a new `plotting-advisor` skill (rules engine, advisor pattern) that recommends chart type, palette, axes, annotation, and accessibility choices for Python plotting. Delegates library mechanics to existing `scientific-skills:matplotlib`/`seaborn`/`plotly` skills. Includes passive lint flow (`scripts/style_lint.py`) for checking existing figures.

## Changes Planned

- [ ] Follow CONFIG_MANAGEMENT.md workflow
- [ ] Create `claude-config/skills/plotting-advisor/` skill scaffolding + SKILL.md
- [ ] Write seven reference docs (`SOURCES`, `chart-selection`, `palettes`, `anti-patterns`, `accessibility`, `interactive-adaptation`, `scientific-conventions`)
- [ ] Write three scripts (`palettes.py`, `figure_spec.py`, `style_lint.py`) with unittest coverage
- [ ] Run `./sync-config.py push --dry-run` and `./sync-config.py push`
- [ ] Verify skill appears at `~/.claude/skills/plotting-advisor/` and YAML parses

## Expected Outcome

`plotting-advisor` is invocable, fires on plotting verbs, returns structured checklist + YAML decision card, and provides lint via `python3 scripts/style_lint.py`.

## Actual Outcome

[Filled in after implementation.]

## Assessment

[Filled in after implementation.]

## Related Commits

- [pending]

## Next Steps

- Update `bioinformatician/SKILL.md` to point its missing visualization reference at this skill's `references/scientific-conventions.md` (separate PR).
