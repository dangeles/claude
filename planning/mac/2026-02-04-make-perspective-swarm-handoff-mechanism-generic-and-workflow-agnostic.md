# Make perspective-swarm handoff mechanism generic and workflow-agnostic

**Date**: 2026-02-04
**Machine**: mac
**Status**: Complete

## Objective

Replace the hardcoded lit-pm-specific handoff mechanism in perspective-swarm with a generic, workflow-agnostic discovery and selection system. Enable perspective-swarm to hand off to any skill that declares handoff eligibility via frontmatter metadata.

## Changes Made

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Update perspective-swarm/SKILL.md - Generic handoff UX in Stage 4
- [x] Update perspective-swarm/references/handoff-schema.md - v2.0 workflow-agnostic schema
- [x] Create perspective-swarm/references/workflow-discovery.md - Discovery algorithm documentation
- [x] Create perspective-swarm/references/handoff-payload.schema.json - JSON Schema validation
- [x] Update perspective-swarm/examples/product-decision-example.md - Multi-option handoff example
- [x] Add handoff metadata to lit-pm/SKILL.md
- [x] Add handoff metadata to programming-pm/SKILL.md
- [x] Add handoff metadata to pov-expansion/SKILL.md

## Key Features Implemented

### Architecture
- **Hybrid discovery**: File-system scan + optional registry
- **Extensible category system**: Core recommendations (research, implementation, analysis, architecture, verification, creative) with free-form extensions
- **Smart relevance scoring**: Based on synthesis content, problem type, uncertainties, and convergence level
- **Project-scoped workflow support**: Git root priority ordering (project > git-root > global)

### Technical Details
- Discovery algorithm with session-scoped caching (30-minute TTL)
- JSON Schema validation for handoff payloads
- Top 3-5 multi-option handoff UX
- Health check and lifecycle hooks support (future)
- Session state preservation for rollback

### Edge Cases Handled
- EC-01: No workflows found (differentiated error messages)
- EC-02: Malformed metadata (schema validation + lenient defaults)
- EC-03: Git root detection failures (submodule handling, timeout handling)

### CRITICAL Fixes Applied
1. JSON Schema workflow_id pattern: `^swarm-session-` (matches session directory format)
2. Category enum: Removed fixed constraint, now extensible
3. Git submodule logic: Fixed inverted conditional for `--show-superproject-working-tree`
4. Handoff chain: Clarified that chain is built from session history, not discovery-time

## Validation Results

### Quality Gate 4 - Pre-Sync Validation
- YAML frontmatter: 4/4 skills PASS
- JSON Schema syntax: PASS
- sync-config.py push: PASS (with conflict resolution)

### Quality Gate 5 - Post-Execution Verification
- Original requirement met: Stage 4 now generic (no hardcoded lit-pm in handoff logic)
- All 8 files modified: PASS
- Edge cases documented: EC-01, EC-02, EC-03 all addressed
- No regressions: PASS

## Actual Outcome

Successfully implemented generic handoff mechanism. perspective-swarm can now:
1. Discover handoff-eligible skills at runtime via frontmatter scan
2. Present relevance-ranked options (top 3-5) based on synthesis content
3. Generate v2.0 handoff payloads validated by JSON Schema
4. Support any future skill that adds `handoff:` metadata

## Assessment

**Result**: Success

**Improvements**:
- perspective-swarm is no longer coupled to lit-pm
- New handoff targets can be added without modifying perspective-swarm
- Relevance scoring provides intelligent recommendations
- JSON Schema ensures payload compatibility

**Issues**:
- None observed during implementation

**Lessons Learned**:
- The adversarial review caught 4 critical issues that would have caused runtime failures
- Line number references in implementation plans drift quickly - verify before editing
- Interactive sync scripts need special handling in non-TTY environments

## Related Commits

- dc84724: feat(perspective-swarm): Make handoff mechanism generic and workflow-agnostic

## Next Steps

1. Monitor handoff usage to refine relevance scoring weights
2. Consider implementing health checks for target skills
3. Add acknowledgment phase for handoff confirmation (future v2.1)
4. Create metadata validation tool: `claude skill validate ./SKILL.md`
