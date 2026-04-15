# Add patent-review skill

**Date**: 2026-04-15
**Machine**: mac
**Status**: In Progress

## Objective

Add a patent-review skill that guides scientists through structured review of lawyer-drafted patent
drafts. The skill covers technical accuracy checking, coverage gap identification, and interactive
proposal of new or revised claim language. Literature review is optional, triggered only at explicit
user request.

## Changes Planned

- [x] Design spec written and approved (3 review iterations; committed fab7177, 733adf3, 0f18d91)
- [x] Eval fixtures created: DMD gene therapy, KRAS inhibitor, SCD base editing
- [x] Benchmark run: with-skill 100% (16/16), without-skill 68% (11/16), Δ+32%
- [ ] Create `claude-config/skills/patent-review/SKILL.md` following plugin-dev:skill-development standards
- [ ] Create `claude-config/skills/patent-review/references/report-template.md`
- [ ] Copy evals fixtures into `claude-config/skills/patent-review/evals/`
- [ ] Run `sync-config.py push` to apply to `~/.claude/`
- [ ] Commit and push to git

## Expected Outcome

Skill available in all Claude Code sessions; correctly triggers on patent review requests; produces
structured interactive review with per-claim status table and concrete revised claim language.

## Actual Outcome

[After implementation: What actually happened?]

## Assessment

**Result**: [Success / Partial / Failed]

**Improvements**:
- Benchmark shows +32% pass rate vs. baseline; two discriminating patterns: per-claim status table
  and concrete revised claim text (baseline produces prose analysis only)

**Issues**:
- Initial skill creation violated CONFIG_MANAGEMENT workflow: created directly in `~/.claude/` instead
  of `claude-config/` first; used skill-creator:skill-creator instead of plugin-dev:skill-development;
  never committed to git. This run corrects those violations.

**Lessons Learned**:
- Follow CONFIG_MANAGEMENT.md before any skill creation; use plugin-dev:skill-development protocol.

## Related Commits

- [pending]

## Next Steps

- Iteration 2 evals: tighten assertions (domain-identification assertion doesn't discriminate; add
  assertion for correctness of scientific claims, not just presence)
