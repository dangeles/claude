# Pre-Mortem Template

## Project Information

**Project Name**: [Project Name]
**Date**: [YYYY-MM-DD]
**Facilitator**: programming-pm
**Participants**: [List of involved skills/stakeholders]

---

## Frame

> "It is 90 days from now and the project has failed. What went wrong?"

Use prospective hindsight to identify risks before they materialize. Each participant independently brainstorms potential failure modes before group discussion.

---

## Risk Identification

### Technical Risks

| Risk ID | Description | Likelihood (1-5) | Impact (1-5) | Risk Score | Category |
|---------|-------------|------------------|--------------|------------|----------|
| T-001 | | | | | |
| T-002 | | | | | |
| T-003 | | | | | |

**Categories**: Architecture, Performance, Security, Integration, Data Quality, Dependencies

**Common Technical Risks**:
- Underestimated complexity of core algorithm
- Performance bottleneck in critical path
- Security vulnerability in input handling
- Integration failure with external API
- Data quality issues in input sources
- Dependency version conflicts

### Process Risks

| Risk ID | Description | Likelihood (1-5) | Impact (1-5) | Risk Score | Category |
|---------|-------------|------------------|--------------|------------|----------|
| P-001 | | | | | |
| P-002 | | | | | |
| P-003 | | | | | |

**Categories**: Communication, Scope, Resources, Timeline, Quality

**Common Process Risks**:
- Scope creep from unclear requirements
- Underestimated task duration
- Code review bottleneck
- Missing specialist expertise
- Inadequate test coverage

### External Risks

| Risk ID | Description | Likelihood (1-5) | Impact (1-5) | Risk Score | Category |
|---------|-------------|------------------|--------------|------------|----------|
| E-001 | | | | | |
| E-002 | | | | | |
| E-003 | | | | | |

**Categories**: Dependencies, Environment, Requirements Change, Resource Availability

**Common External Risks**:
- Third-party API changes or downtime
- Hardware/environment constraints discovered late
- Requirements change mid-implementation
- Key resource unavailability

---

## Risk Scoring Guide

### Likelihood Scale
| Score | Label | Description |
|-------|-------|-------------|
| 1 | Rare | < 10% chance |
| 2 | Unlikely | 10-25% chance |
| 3 | Possible | 25-50% chance |
| 4 | Likely | 50-75% chance |
| 5 | Almost Certain | > 75% chance |

### Impact Scale
| Score | Label | Description |
|-------|-------|-------------|
| 1 | Negligible | Minor inconvenience, easily resolved |
| 2 | Minor | Some rework required, < 1 day delay |
| 3 | Moderate | Significant rework, 1-3 day delay |
| 4 | Major | Critical path affected, > 3 day delay |
| 5 | Severe | Project failure or fundamental redesign |

### Risk Score Interpretation
- **Risk Score** = Likelihood x Impact
- **1-4**: Low risk (monitor)
- **5-9**: Medium risk (mitigation plan recommended)
- **10-14**: High risk (mitigation plan required)
- **15-25**: Critical risk (must address before proceeding)

---

## Risk Disposition

For each identified risk, select one disposition:

| Disposition | Definition | When to Use |
|-------------|------------|-------------|
| **Mitigate** | Take action to reduce likelihood or impact | Risk is controllable, cost-effective to address |
| **Accept** | Acknowledge and monitor, no action | Low risk OR high cost to mitigate |
| **Transfer** | Assign ownership outside project | Risk is outside project control |
| **Avoid** | Change scope to eliminate risk | Risk is too high, can redesign around it |

### Risk Response Table

| Risk ID | Disposition | Rationale | Owner | Mitigation/Notes | Review Trigger |
|---------|-------------|-----------|-------|------------------|----------------|
| | | | | | |
| | | | | | |
| | | | | | |

---

## Mitigation Plans

### For each risk with Mitigate disposition:

**Risk ID**: T-001
**Risk**: [Description]
**Mitigation Strategy**:
1. [Action 1]
2. [Action 2]
3. [Action 3]

**Success Criteria**: [How we know mitigation worked]
**Contingency**: [What to do if mitigation fails]
**Owner**: [Who is responsible]
**Timeline**: [When actions occur]

---

## Accepted Risks Summary

Document all accepted risks for the record:

```yaml
accepted_risks:
  - id: "RISK-XXX"
    description: ""
    severity: "High" | "Medium" | "Low"
    user_decision: "accept"
    rationale: ""
    accepted_date: ""
    review_trigger: ""  # When to re-evaluate
```

### Accepted Risk Acknowledgment

| Risk ID | Description | Severity | Accepted By | Date | Review Trigger |
|---------|-------------|----------|-------------|------|----------------|
| | | | | | |

---

## Pre-Mortem Completion Checklist

### Minimum Requirements (Quality Gate 2)

- [ ] At least 3 risks identified across all categories
- [ ] All risks have likelihood rating (1-5)
- [ ] All risks have impact rating (1-5)
- [ ] All risks have disposition (mitigate/accept/transfer/avoid)
- [ ] Critical risks (score >= 15) have contingency plans

### Recommended

- [ ] At least 1 risk from each category (Technical, Process, External)
- [ ] All mitigations have owners assigned
- [ ] All mitigations have success criteria
- [ ] Accepted risks documented with rationale
- [ ] Review triggers set for accepted risks

---

## Sign-Off

**Pre-Mortem Complete**: [ ] Yes / [ ] No (explain)

**Date**: [YYYY-MM-DD]

**Notes**:
[Any additional context or concerns]

---

## Post-Implementation Review

After project completion, revisit this pre-mortem:

| Risk ID | Occurred? | Actual Impact | Mitigation Effective? | Lessons Learned |
|---------|-----------|---------------|----------------------|-----------------|
| | | | | |

This feedback improves future pre-mortem accuracy.
