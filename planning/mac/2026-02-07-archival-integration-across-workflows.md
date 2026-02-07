# Archival Integration Across Workflows

**Date**: 2026-02-07
**Machine**: mac
**Status**: Complete

## Objective

Integrate archival guidelines compliance across all 16 file-generating workflow skills using a centralized, machine-readable .archive-metadata.yaml configuration produced by archive-workflow (single producer) and consumed by all other workflows. Eliminates the DRY violation where each skill independently parsed CLAUDE.md for guidelines.

## Changes Planned

- [x] Create centralized archival-compliance-check.md reference document
- [x] Create naming-conventions-mixed.md and structure-template-mixed.md references
- [x] Extend archive-workflow SKILL.md with Wave 4 .archive-metadata.yaml generation
- [x] Replace lit-pm Stage 0 with .archive-metadata.yaml primary source
- [x] Replace programming-pm Phase 0 with .archive-metadata.yaml primary source
- [x] Add archival check to research-pipeline Step 1
- [x] Add archival awareness check to skill-editor pre-workflow
- [x] Add archival check to scientific-analysis-architect Phase 0
- [x] Add archival compliance to 7 specialist skills
- [x] Add archival compliance to 3 developer skills
- [x] Validate all YAML frontmatter (Quality Gate 4)
- [x] Sync to ~/.claude/ (Quality Gate 5)

## Architecture

Hub-and-spoke data-mediated coordination pattern:
- Single producer: archive-workflow generates .archive-metadata.yaml
- 16 consumers: all file-generating skills read and comply
- Centralized compliance check: one reference document, all skills reference it
- Showstopper fix: all consumers use absolute paths (~/.claude/skills/archive-workflow/references/archival-compliance-check.md) with graceful degradation

## Key Decisions

1. .archive-metadata.yaml as machine-readable replacement for CLAUDE.md parsing
2. CLAUDE.md retained as deprecated fallback (removal in v2.0)
3. Advisory enforcement mode as default (never blocks workflow)
4. Orchestrators pass archival_context in handoffs; specialists check handoff first
5. skill-editor gets awareness-only check (writes to claude-config/, not project dirs)
6. Circular dependency prevention via archival_context: "skip" in handoffs

## Files Modified

### New Files (3)
- claude-config/skills/archive-workflow/references/archival-compliance-check.md
- claude-config/skills/archive-workflow/references/naming-conventions-mixed.md
- claude-config/skills/archive-workflow/references/structure-template-mixed.md

### Modified Files (16)
- claude-config/skills/archive-workflow/SKILL.md (Wave 4 extension)
- claude-config/skills/lit-pm/SKILL.md (Stage 0 replacement)
- claude-config/skills/programming-pm/SKILL.md (Phase 0 replacement)
- claude-config/skills/research-pipeline/SKILL.md (Step 1 addition)
- claude-config/skills/skill-editor/SKILL.md (pre-workflow awareness)
- claude-config/skills/scientific-analysis-architect/SKILL.md (Phase 0 addition)
- claude-config/skills/technical-pm/SKILL.md (archival compliance section)
- claude-config/skills/researcher/SKILL.md (archival compliance section)
- claude-config/skills/literature-researcher/SKILL.md (archival compliance section)
- claude-config/skills/editor/SKILL.md (archival compliance section)
- claude-config/skills/synthesizer/SKILL.md (archival compliance section)
- claude-config/skills/lit-synthesizer/SKILL.md (archival compliance section)
- claude-config/skills/fact-checker/SKILL.md (archival compliance section)
- claude-config/skills/software-developer/SKILL.md (archival compliance section)
- claude-config/skills/senior-developer/SKILL.md (archival compliance section)
- claude-config/skills/junior-developer/SKILL.md (archival compliance section)

## Testing

- All 41 skills pass YAML frontmatter validation
- Dry-run sync succeeds (only pre-existing installed_plugins.json divergence)
- Smoke test: all modified skills load correctly

## Known Limitations

- CLAUDE.md fallback not yet removed (planned for v2.0)
- No automated end-to-end test for .archive-metadata.yaml round-trip
- Schema versioning policy defined but v1.1 enhancements not yet planned
