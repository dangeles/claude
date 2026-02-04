# Create lit-pm orchestrator skill for literature review pipeline

**Date**: 2026-02-03
**Machine**: mac
**Status**: Complete

## Objective

Create an 8-stage literature review orchestrator skill (lit-pm) that coordinates comprehensive literature reviews with adaptive checkpoints, parallel execution, and robust operational mechanisms. The skill manages handoffs between literature-researcher, lit-synthesizer, fact-checker, and editor skills to produce decision-useful literature reviews.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Create lit-pm skill directory structure
- [x] Create SKILL.md with 8-stage pipeline specification (~400 lines)
- [x] Create references/stage-specifications.md (~200 lines)
- [x] Create references/handoff-schema.md (~100 lines)
- [x] Create references/error-handling.md (~150 lines)
- [x] Create references/adaptive-orchestration.md (~100 lines)
- [x] Create references/quality-gates.md (~100 lines)
- [x] Create examples/hepatocyte-review-example.md (~300 lines)

## Expected Outcome

A fully specified orchestrator skill that:
- Coordinates 8-stage literature review pipeline
- Implements adaptive checkpoints based on complexity and stakes
- Supports parallel execution for review discovery (Stage 2) and section writing (Stage 5)
- Provides robust error handling with Saga-style compensation
- Manages workflow state for resumption after interrupts
- Uses YAML handoff protocols between stages

## Actual Outcome

Successfully created lit-pm skill with all planned files:

**Files Created**:
- `claude-config/skills/lit-pm/SKILL.md` (361 lines)
- `claude-config/skills/lit-pm/references/stage-specifications.md` (16KB)
- `claude-config/skills/lit-pm/references/handoff-schema.md` (8.3KB)
- `claude-config/skills/lit-pm/references/error-handling.md` (9.3KB)
- `claude-config/skills/lit-pm/references/adaptive-orchestration.md` (9.6KB)
- `claude-config/skills/lit-pm/references/quality-gates.md` (9.6KB)
- `claude-config/skills/lit-pm/examples/hepatocyte-review-example.md` (21KB)

**Total**: ~74KB across 7 files

**Key Features Implemented**:
1. 8-stage pipeline (Scope, Reviews, Outline, Intro, Sections, Fact-Check, Synthesis, Polish)
2. 4-tier complexity detection (Simple, Medium, Complex, High-Stakes)
3. Adaptive checkpoint plans based on complexity
4. Parallel execution protocols for Stages 2 and 5
5. Timeout configuration per stage and per agent
6. Resource limits (max 4 concurrent agents)
7. Saga-style compensation matrix for rollback
8. Circuit breaker configuration for failure handling
9. Atomic state write protocol for crash recovery
10. Escape hatch protocol for infinite loop prevention
11. Quality floor checks that cannot be overridden

## Assessment

**Result**: Success

**Improvements**:
- New orchestrator skill enables comprehensive literature reviews
- Progressive disclosure keeps SKILL.md under 500 lines while providing full detail in references
- Robust error handling addresses all 8 critical edge cases identified in analysis phase
- Example walkthrough demonstrates complete workflow

**Issues**:
- None encountered during implementation

**Lessons Learned**:
- Progressive disclosure (SKILL.md + reference files) is effective pattern for complex skills
- Pre-flight validation catches dependency issues before workflow begins
- Explicit timeout and resource limit configuration prevents runaway execution

## Related Commits

- [pending]: feat(lit-pm): Create literature review pipeline orchestrator

## Next Steps

- Create literature-researcher skill (enhanced from researcher)
- Create lit-synthesizer skill (senior scientific author personality)
- Integration test full pipeline after dependencies exist
- Calibrate timeouts based on actual execution durations
