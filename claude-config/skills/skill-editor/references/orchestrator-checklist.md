# Orchestrator Best Practices Evaluation Checklist

**Version**: 1.0
**Usage**: fill in this checklist during Phase 2 analysis when the target skill is detected
as an orchestrator. Pattern definitions and templates live in the orchestrator
best-practices reference listed in SKILL.md.

## Target Information

- **Skill Name**: ___
- **Target Type**: [ ] New orchestrator | [ ] Existing orchestrator modification
- **Architecture Type**: [ ] Phase-based | [ ] Stage-based | [ ] Event-driven | [ ] Other: ___

---

## Core Patterns (widely adopted)

| # | Pattern | Status | Evidence |
|---|---------|--------|----------|
| 1 | Delegation Mandate | [ ] PRESENT [ ] PARTIAL [ ] ABSENT | _citation_ |
| 2 | State Anchoring | [ ] PRESENT [ ] PARTIAL [ ] ABSENT | _citation_ |
| 3 | Tool Selection Table | [ ] PRESENT [ ] PARTIAL [ ] ABSENT | _citation_ |
| 4 | Error Handling (Structured) | [ ] PRESENT [ ] PARTIAL [ ] ABSENT | _citation_ |
| 5 | Timeout Configuration | [ ] PRESENT [ ] PARTIAL [ ] ABSENT | _citation_ |
| 6 | Session Management | [ ] PRESENT [ ] PARTIAL [ ] ABSENT | _citation_ |

**Core coverage**: ___/6

---

## Common Patterns

| # | Pattern | Status | Evidence |
|---|---------|--------|----------|
| 7 | Handoff Schema | [ ] PRESENT [ ] PARTIAL [ ] ABSENT [ ] N/A | _citation_ |
| 8 | Pre-Flight Validation | [ ] PRESENT [ ] PARTIAL [ ] ABSENT [ ] N/A | _citation_ |
| 9 | Mode Selection | [ ] PRESENT [ ] PARTIAL [ ] ABSENT [ ] N/A | _citation_ |
| 10 | Handoff Frontmatter | [ ] PRESENT [ ] PARTIAL [ ] ABSENT [ ] N/A | _citation_ |

**Common coverage**: ___/4 (excluding N/A)

N/A is valid when an architectural mismatch exists (e.g., event-driven orchestrator lacks sequential phases for Mode Selection).

---

## Emerging Patterns (noted, not scored)

| # | Pattern | Status | Notes |
|---|---------|--------|-------|
| 11 | Quality Gate Override | [ ] PRESENT [ ] ABSENT | _optional notes_ |

---

## Assessment

Coverage is descriptive, not a grade. Report the two counts, then say which absent patterns
actually matter for this orchestrator's architecture and why. A pattern that does not fit
the architecture is not a gap.

**Overall Assessment**: ___

---

## Notes

- For **existing orchestrators**: PARTIAL with a working variant is acceptable. Avoid recommending that working implementations be replaced with standard templates unless the variant has known issues.
- For **new orchestrators**: note absent core patterns so the implementation plan can consider them.
- **Sanity check**: if none of the six core patterns are present, reconsider whether the target is actually an orchestrator. This may indicate a false positive in detection.
- **Evidence citations**: one sentence per pattern citing the specific section, line range, or feature that demonstrates presence/absence.
