# Risk Matrix for Technical PMs

**Last Updated**: 2026-01-28
**Purpose**: Framework for identifying, assessing, and mitigating project risks

---

## Likelihood × Impact Grid

### Risk Scoring

| | **Low Impact** | **Medium Impact** | **High Impact** |
|---|---|---|---|
| **High Likelihood** | 🟨 Monitor | 🟧 Mitigate | 🟥 Act Now |
| **Medium Likelihood** | 🟩 Accept | 🟨 Monitor | 🟧 Mitigate |
| **Low Likelihood** | 🟩 Accept | 🟩 Accept | 🟨 Monitor |

**Legend**:
- 🟥 **Act Now**: Immediate mitigation required, escalate to user if needed
- 🟧 **Mitigate**: Proactive action needed, include mitigation plan in work plan
- 🟨 **Monitor**: Watch closely, prepare mitigation if likelihood or impact increases
- 🟩 **Accept**: Document risk, revisit if conditions change

### Likelihood Definitions

- **High**: >50% chance of occurring during project
- **Medium**: 20-50% chance
- **Low**: <20% chance

### Impact Definitions

- **High**: Delays critical path >2 days OR requires major scope change OR blocks deliverable
- **Medium**: Delays critical path 1-2 days OR requires rework >4 hours
- **Low**: Minor inconvenience, <1 day impact, no scope change needed

---

## Common Project Risks and Mitigations

### 1. Scope Creep
**Symptom**: "Quick review" becomes comprehensive analysis; tasks expand 50-200% beyond initial estimate
**Likelihood**: High (especially on novel research)
**Impact**: Medium-High (timeline slips, resources overcommitted)

**Mitigations**:
- ✅ Define scope boundaries explicitly: "ONLY X, exclude Y"
- ✅ Include "Out of Scope" section in task descriptions
- ✅ Set clear deliverable format (e.g., "2-page summary, max 10 citations")
- ✅ Use MVP thinking: What's minimum needed to answer question?
- ✅ Time-box exploratory work: "Spend max 2 hours on tangent, then decide"

**Escalation trigger**: If scope expands >30%, ask user if new scope is priority or revert to original

---

### 2. Agent Spinning / Timeout
**Symptom**: Agent works >30 min without progress update; stuck in loop or overly broad task
**Likelihood**: Medium (especially on ambiguous or open-ended tasks)
**Impact**: Medium (wastes time, blocks downstream work)

**Mitigations**:
- ✅ Set up progress monitoring: Check progress files every 30-60 min for tasks >1 hour
- ✅ Define concrete milestones: Not "research topic" but "find 3 papers with parameter X"
- ✅ Narrow scope when agent spins: Reduce parameters, shorten timeline, focus on subset
- ✅ Use timeout intervention protocol: Analyze → Notify → Options → Execute

**Escalation trigger**: Agent >30 min without update → use timeout intervention; if scope unclear, ask user for prioritization

---

### 3. Dependency Blocking
**Symptom**: Task B waiting for Task A, which is delayed; Task B agent idle or reassigned
**Likelihood**: Medium (in multi-agent workflows with sequential dependencies)
**Impact**: Medium (timeline slip, resource inefficiency)

**Mitigations**:
- ✅ Create parallel tracks when possible: Identify tasks that DON'T depend on each other
- ✅ Add buffer tasks: If critical path task delayed, assign agent to lower-priority task instead of idle
- ✅ Explicit handoff timing: "Task B starts AFTER you've reviewed Task A output", not just "depends on A"
- ✅ Front-load uncertainty: Do risky/uncertain tasks first to minimize late surprises

**Escalation trigger**: Critical path task delayed >1 day → notify user of timeline impact, present options (narrow scope, defer features, accept delay)

---

### 4. Quality Issues (Insufficient Evidence, Citation Errors)
**Symptom**: Fact-Checker finds 5+ errors; Devil's Advocate requires 2+ iterations; values don't match primary sources
**Likelihood**: Medium (especially on fast-turnaround work or low-quality sources)
**Impact**: Medium-High (rework required, delays publication, credibility damage)

**Mitigations**:
- ✅ Use Devil's Advocate pairing for substantive documents (mandatory workflow)
- ✅ Fact-Checker verification before handoff to integration/synthesis
- ✅ Enforce source quality hierarchy: Prefer high-impact journals, cross-reference claims
- ✅ Citation format checks: Verify DOI, quote accuracy, table/figure references

**Escalation trigger**: If >3 citation errors or major claim unsupported, pause downstream work and fix evidence base first

---

### 5. Timeline Pressure (Deadline at Risk)
**Symptom**: Critical path tasks delayed; completion date slipping toward or past deadline
**Likelihood**: Medium (especially on multi-week projects with external deadlines)
**Impact**: High (missed deadline, rushed work, quality compromise)

**Mitigations**:
- ✅ Build 20-30% buffer at project level (not communicated as deadline, but internal target)
- ✅ Identify must-have vs. nice-to-have: Ruthlessly prioritize core deliverables
- ✅ Narrow scope aggressively: Drop analysis depth, reduce parameter sweeps, defer tangential questions
- ✅ Increase parallelization: Can multiple agents work simultaneously if we split task?
- ✅ Trade breadth for depth: Better to answer one question well than three questions poorly

**Escalation trigger**: If deadline at risk despite all mitigations, notify user EARLY (1 week+ before deadline) with options: (A) Extend deadline, (B) Reduce scope, (C) Accept reduced quality

---

### 6. Resource Conflict (Agent Double-Booked)
**Symptom**: Researcher needed for Task A and Task B simultaneously; work stalls
**Likelihood**: Low (with good PM coordination)
**Impact**: Medium (delays one or both tasks)

**Mitigations**:
- ✅ Stagger tasks: Schedule Task A and B sequentially instead of parallel
- ✅ Agent substitution: Can Calculator handle literature parameter extraction? Can Synthesizer do light research?
- ✅ Priority ranking: If conflict unavoidable, which task is critical path?

**Escalation trigger**: If both tasks are critical path and staggering isn't viable, ask user which takes priority

---

### 7. Missing Information / Data Unavailable
**Symptom**: Agent can't find parameter values, papers paywalled, user input needed but unresponsive
**Likelihood**: Medium (especially on novel research or niche topics)
**Impact**: Medium-High (task blocked, potentially indefinitely)

**Mitigations**:
- ✅ Document search dead-ends: "Searched PubMed with X keywords, found nothing" → prevents repeated failed searches
- ✅ Use proxies or ranges: If exact value unknown, can we bound it? (e.g., "literature reports 5-15, we'll use 10")
- ✅ Escalate information needs early: Don't spin for 2 hours, ask user after 30 min if they have data source
- ✅ Compile paywall list proactively: Send to user in batch for access, don't block per-paper

**Escalation trigger**: If information unavailable after exhaustive search (30+ min), use AskUserQuestion with summary of what you tried and what's needed

---

### 8. Agent Output Quality Mismatch (Doesn't Meet Downstream Needs)
**Symptom**: Synthesizer receives Calculator output but it lacks context; Editor receives draft but doesn't know which sections are critical
**Likelihood**: Low-Medium (higher on first handoff in new workflow)
**Impact**: Medium (rework, communication overhead, delays)

**Mitigations**:
- ✅ Handoff ceremony: Include (1) deliverable, (2) context (why this matters), (3) what to check/focus on, (4) known gaps
- ✅ Pre-handoff review: PM quickly scans output before passing to next agent—does it answer downstream needs?
- ✅ Explicit format requirements: In task description, specify deliverable structure ("include parameter table with citations")

**Escalation trigger**: If handoff breakdown requires >2 hours rework, pause and clarify requirements before proceeding

---

## Risk Escalation Criteria

**When to escalate to user** (use AskUserQuestion):

| Situation | Escalate If | Example |
|---|---|---|
| **Timeline slip** | Critical path delayed >1 day | "Detailed model taking 2 days instead of 1, synthesis delayed to Feb 5" |
| **Scope change** | Expansion >30% effort | "Oxygen analysis now includes CO₂ removal (adds 8 hours)" |
| **Resource conflict** | Two critical-path tasks need same agent | "Researcher needed for both literature review and parameter extraction simultaneously" |
| **Quality vs. speed tradeoff** | Unclear user priority | "Comprehensive analysis takes 2 weeks, quick assessment takes 2 days—which matters more?" |
| **Missing information** | Unavailable after 30+ min search | "Can't find hepatocyte ammonia clearance rates in literature; do you have experimental data?" |
| **Risk mitigation requires user decision** | Options have different strategic implications | "Narrow oxygen model to 3 parameters (faster) or full 6-parameter sweep (comprehensive)?" |

**When NOT to escalate** (handle internally):

- Low/Low risks (accept and document)
- Routine process decisions (which agent does task, when to schedule)
- Minor timeline adjustments (<1 day on non-critical path)
- Technical implementation details (unless user expertise needed)

---

## Risk Review Cadence

### Daily (for active sprints)
- Review progress dashboard for status changes (✅ → ⚠️ → 🔴)
- Check for new blockers or timeline slips
- Update risk table if likelihood/impact changes

### Weekly (for ongoing projects)
- Comprehensive risk review: Are mitigations working?
- Velocity check: Are estimates accurate? Adjust if needed
- User communication: Dashboard update with risks/mitigations

### Monthly (for long-term initiatives)
- Post-mortem on completed work: What risks materialized? What didn't?
- Update common risk library: Add new patterns, refine mitigations
- Estimation model refinement: Adjust duration ranges based on actuals

---

## Risk Documentation Template

Use this format in work plans and dashboards:

```markdown
## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation | Owner | Status |
|------|-----------|--------|-----------|-------|--------|
| [Risk name] | [H/M/L] | [H/M/L] | [What you're doing to reduce risk] | [Who's monitoring] | [Current state] |
```

**Example**:

| Risk | Likelihood | Impact | Mitigation | Owner | Status |
|------|-----------|--------|-----------|-------|--------|
| Sensitivity analysis overruns timeline | Medium | Medium | Narrow scope to 3 critical parameters (saves 1 day with minimal information loss) | Technical PM | Monitoring - decision pending |
| Synthesis document requires multiple Devil's Advocate iterations | Low | Medium | Ensure Calculator deliverables are complete with clear conclusions before handoff | Technical PM | Mitigated - Calculator output reviewed |

---

## References

- **Risk Management Framework**: https://www.pmi.org/learning/library/risk-analysis-project-management-7070
- **RAID Log** (Risks, Assumptions, Issues, Dependencies): Common PM tool for tracking
