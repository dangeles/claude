# Estimation Frameworks for Technical PMs

**Last Updated**: 2026-01-28
**Purpose**: Quick reference for task estimation, velocity tracking, and timeline planning

---

## T-Shirt Sizing (Quick Complexity Assessment)

Use for initial scoping when precise hours aren't critical.

| Size | Complexity | Typical Duration | Examples |
|------|-----------|------------------|----------|
| **XS** | Trivial | <30 min | Update index, fix typo, single-value fact-check |
| **S** | Simple | 30 min - 2 hours | Single paper notes, back-of-envelope calc, short document edit |
| **M** | Moderate | 2-6 hours | Comprehensive review (3-5 papers), synthesis document, detailed model |
| **L** | Complex | 6-16 hours (1-2 days) | Multi-track work plan, extensive literature search (10+ papers), multi-agent coordination |
| **XL** | Very complex | 16+ hours (2+ days) | Novel research area exploration, full system design, multi-week sprint planning |

**When to use**: Initial work breakdown, prioritization discussions, rough timeline sketches.

**When NOT to use**: Final timeline commitments, resource allocation, detailed scheduling.

---

## PERT Estimation (Three-Point Estimates)

Use when uncertainty is high and you need confidence intervals.

**Formula**: Expected time = (Optimistic + 4×Most Likely + Pessimistic) / 6

**Example - Literature Review Task**:
- **Optimistic** (best case): 3 hours - papers are highly relevant, well-cited, easy to find
- **Most Likely**: 5 hours - typical search complexity, standard paper quality
- **Pessimistic** (worst case): 8 hours - need to access paywalled papers, conflicting data, extensive citation chaining

**Expected time** = (3 + 4×5 + 8) / 6 = **5.2 hours**

**When to use**: Tasks with significant uncertainty, first-time activities, external dependencies (paywall access, user input).

**How to communicate**: "Expected 5 hours, could range 3-8 hours depending on paper accessibility."

---

## Common Task Duration Ranges by Agent

### Researcher
- **Single paper notes**: 1-3 hours (depends on paper length, complexity, citation density)
- **Focused review (3-5 papers)**: 4-8 hours
- **Comprehensive review (10+ papers)**: 12-20 hours (typically 2-3 days)
- **Parameter extraction**: 1-2 hours (from existing review)
- **PDF acquisition**: 15-30 min per batch of 5-10 papers

**Estimation factors**:
- Paper length (10-page paper = 1 hr, 50-page review = 3-4 hrs)
- Citation chaining depth (1 level = +30%, 2 levels = +60%)
- Paywall issues (+50% if extensive paywall navigation needed)

### Calculator
- **Back-of-envelope calculation**: 30 min - 2 hours (depends on parameter availability)
- **Detailed model (1D)**: 3-6 hours (model setup, validation, parameter sweeps)
- **Detailed model (2D/3D)**: 8-16 hours (computational complexity, meshing, sensitivity analysis)
- **Sensitivity analysis**: +50% to detailed model time (per 3-parameter sweep)

**Estimation factors**:
- Parameter availability (known = fast, need literature search = +50%)
- Model complexity (algebraic = 1x, ODEs = 2x, PDEs = 4x)
- Validation needs (order-of-magnitude check = fast, detailed validation = +50%)

### Synthesizer
- **Synthesis document (2-3 sources)**: 2-3 hours
- **Synthesis document (4-6 sources)**: 4-6 hours
- **Comprehensive synthesis (7+ sources)**: 8-12 hours

**Estimation factors**:
- Source agreement (convergent = 1x, divergent with conflicts = 1.5x)
- Citation density requirements (+30% for high-citation contexts)

### Editor
- **Light edit (clarity, flow)**: 30 min per 1000 words
- **Medium edit (restructure sections)**: 1 hour per 1000 words
- **Heavy edit (major rewrite)**: 2 hours per 1000 words

**Estimation factors**:
- Draft quality (clean = light, rough = heavy)
- Technical density (prose-heavy = fast, equation-heavy = slow)

### Fact-Checker
- **Quick verification (5-10 citations)**: 15-30 min
- **Thorough verification (20+ citations)**: 1-2 hours
- **Deep verification (quote accuracy, table values)**: 2-4 hours

**Estimation factors**:
- Citation accessibility (DOI links = fast, need PDF hunt = +50%)
- Quantitative vs. qualitative claims (numbers = slower, need exact table/figure match)

### Devil's Advocate
- **Single-pass review**: 1-2 hours (depends on document length and claim density)
- **Iterative exchange**: +1 hour per round (typically 2 rounds max)

**Estimation factors**:
- Document quality (well-supported = fast, weak evidence = longer)
- Claim novelty (standard claims = fast, novel claims need deeper challenge)

### Technical PM
- **Work breakdown (simple, 3-5 tasks)**: 30 min - 1 hour
- **Work breakdown (complex, multi-track, 8+ tasks)**: 1-2 hours
- **Daily progress monitoring**: 15-30 min per day
- **Dashboard updates**: 30 min per week

---

## Story Points vs. Hours

### When to Use Hours
- Short-term tasks (<1 week)
- Single-agent work
- Familiar task types with historical data
- Timeline-critical work requiring precision

### When to Use Story Points
- Long-term planning (multi-week sprints)
- Multi-agent coordination (relative complexity more useful than hours)
- Novel tasks with high uncertainty
- Capacity planning across team

**Conversion (rough guideline)**:
- 1 point = 2-4 hours (small task)
- 2 points = 4-8 hours (medium task)
- 3 points = 8-16 hours (large task)
- 5 points = 16+ hours (very large, consider breaking down)

---

## Velocity Tracking

### Metrics to Track

**Completion rate**: Tasks completed on time / Total tasks assigned
- **Target**: >80% (indicates realistic estimation)
- **If <70%**: Estimates too optimistic, increase task durations by 25%
- **If >95%**: Estimates too conservative, might be over-planning

**Estimation accuracy**: Actual time / Estimated time
- **Target**: 0.8 - 1.2 (within 20% of estimate)
- **If <0.8**: Overestimating (faster than expected)
- **If >1.5**: Underestimating (tasks taking 50% longer)

**Blocker frequency**: Blocked task-days / Total task-days
- **Target**: <10% (most tasks progress smoothly)
- **If >20%**: Dependency issues, need better sequencing or parallel work

### How to Improve Estimates Over Time

1. **Track actuals**: Record actual completion time for each task type
2. **Calculate ratios**: Actual/Estimated for each agent type
3. **Adjust future estimates**: If Calculator tasks average 1.3x estimates, multiply future Calculator estimates by 1.3
4. **Review monthly**: Recalculate adjustment factors every 10-20 tasks

**Example**:
- Researcher literature reviews: 5 completed, average 1.2x estimate → multiply future reviews by 1.2
- Calculator back-of-envelope: 3 completed, average 0.9x estimate → estimates are good, no adjustment

---

## Estimation Errors to Avoid

### 1. Anchoring Bias
**Symptom**: User says "this should be quick" → you estimate 1 hour, actually takes 4 hours
**Fix**: Ignore user time expectations; estimate based on task complexity and historical data

### 2. Planning Fallacy
**Symptom**: "Best case is 2 hours, so I'll estimate 2 hours" → actually takes 5 hours
**Fix**: Use PERT estimation (factor in most likely and pessimistic scenarios)

### 3. Sunk Cost Fallacy
**Symptom**: Task is 80% done but scope expanded; reluctant to re-estimate upward
**Fix**: Re-estimate remaining work based on current state, ignore past effort

### 4. Scope Creep (Undetected)
**Symptom**: "Quick review" estimate, but agent finds interesting tangent and spends 3x time
**Fix**: Define scope boundaries explicitly ("ONLY oxygen transport, exclude clinical outcomes")

### 5. Dependency Blindness
**Symptom**: Estimate task in isolation, but actually blocked waiting for handoff
**Fix**: Add dependency wait time (if Task B needs Task A output and A takes 2 days, B can't start until Day 3)

### 6. No Buffer for Unknowns
**Symptom**: Sum individual task estimates to get project timeline, no slack for unexpected issues
**Fix**: Add 20-30% buffer at project level (not task level) for unknowns

---

## Quick Estimation Checklist

Before committing to timeline, verify:

- [ ] Task scope clearly defined (what's IN scope, what's OUT)
- [ ] Dependencies identified (what must complete first?)
- [ ] Agent availability confirmed (not double-booked)
- [ ] Historical data consulted (how long did similar tasks take?)
- [ ] Uncertainty quantified (best/likely/worst case)
- [ ] Buffer included (20-30% for unknowns at project level)
- [ ] User expectations managed (communicate range, not single number)

---

## References

- **PERT Estimation**: https://en.wikipedia.org/wiki/Program_evaluation_and_review_technique
- **Story Points**: https://www.atlassian.com/agile/project-management/estimation
- **Velocity Tracking**: https://www.scrum.org/resources/what-is-velocity
