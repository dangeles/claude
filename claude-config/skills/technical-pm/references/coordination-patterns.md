# Coordination Patterns for Technical PMs

**Last Updated**: 2026-01-28
**Purpose**: Common agent workflows, handoff best practices, and blocker resolution strategies

---

## Common Agent Workflows

### Pattern 1: Literature Pipeline (Standard Research)

**Flow**: Researcher → Devil's Advocate → Fact-Checker → Editor → archive-workflow

**Use when**: Producing literature review, analysis document, synthesis report

**Stages**:
1. **Researcher** writes draft (review, analysis, paper notes)
2. **Devil's Advocate** challenges assumptions, identifies gaps (1-2 iterations)
3. **Fact-Checker** verifies citations, checks quoted values
4. **Editor** polishes prose, enforces CLAUDE.md style
5. **archive-workflow** organizes project structure, ensures proper naming/location

**Parallel opportunities**: None (sequential by design for quality)

**Estimated duration**:
- Single paper notes: 2-4 hours total (1.5hr research + 30min fact-check + 30min edit)
- Comprehensive review: 8-12 hours total (6hr research + 2hr DA + 1hr fact-check + 1hr edit)

---

### Pattern 2: Analysis Pipeline (Quantitative Feasibility)

**Flow**: Researcher (parameters) → Calculator → Synthesizer → Devil's Advocate → Editor

**Use when**: Assessing technical feasibility with quantitative models

**Stages**:
1. **Researcher** extracts parameters from literature
2. **Calculator** runs back-of-envelope calculation (quick feasibility)
   - If feasible → proceed to detailed model
   - If infeasible by >10x → pivot strategy (don't waste time on detailed model)
3. **Calculator** builds detailed model (if back-of-envelope showed feasibility)
4. **Synthesizer** integrates literature + calculations into decision document
5. **Devil's Advocate** challenges assumptions and model validity
6. **Editor** polishes final document

**Parallel opportunities**:
- While Calculator works on detailed model, Researcher can start related literature work for next phase

**Estimated duration**:
- Parameter extraction: 1-2 hours
- Back-of-envelope: 1-2 hours
- Detailed model: 4-8 hours
- Synthesis + review + edit: 4-6 hours
- **Total: 10-18 hours** (2-3 days)

---

### Pattern 3: Multi-Source Synthesis (Integrating Diverse Evidence)

**Flow**: Researcher (source 1) || Researcher (source 2) || Researcher (source 3) → Synthesizer → Devil's Advocate → Editor

**Use when**: Need to integrate findings from multiple independent sources (e.g., clinical data + in vitro studies + animal models)

**Stages**:
1. **Researchers** (can be same agent, sequential) gather evidence from different source types
2. **Synthesizer** integrates across sources, identifies convergence/divergence
3. **Devil's Advocate** challenges integration logic, flags weak synthesis
4. **Editor** polishes narrative flow

**Parallel opportunities**:
- Multiple researchers can work on different sources simultaneously if independent
- Example: Researcher A handles clinical papers, Researcher B handles in vitro papers (parallel)

**Estimated duration**:
- 3 sources × 2hr each = 6 hours research (or 6hr sequential, 2-3hr parallel)
- Synthesis: 3-4 hours
- Review + edit: 2-3 hours
- **Total: 11-13 hours** (can reduce to 7-10 hours if parallel research)

---

### Pattern 4: Rapid Iteration (Debugging, Refinement)

**Flow**: Calculator → Fact-Checker (quick check) → Calculator (revise) → repeat until converged

**Use when**: Model or calculation has errors, needs iterative refinement

**Stages**:
1. **Calculator** produces initial model/calculation
2. **Fact-Checker** quickly verifies units, order of magnitude, logic
3. **Calculator** fixes issues
4. Repeat until Fact-Checker confirms correctness (usually 1-2 cycles)
5. Then proceed to standard quality pipeline (Devil's Advocate, Editor)

**Parallel opportunities**: None (iterative by nature)

**Estimated duration**:
- 1-2 cycles × 1-2hr each = 2-4 hours
- Then +3-4hr for full quality pipeline

---

## Parallel vs. Sequential Work: Decision Tree

```
Does Task B need Task A's output?
│
├─ YES → Sequential (Task A must complete first)
│   │
│   └─ Can Task B do preparatory work while waiting?
│       ├─ YES → Start Task B "prep phase" in parallel, pause for handoff, resume
│       └─ NO → Task B waits (assign agent to buffer task)
│
└─ NO → Parallel (Task A and B independent)
    │
    └─ Do Task A and B need the same agent?
        ├─ YES → Sequential (agent can't clone) OR find substitute agent
        └─ NO → Fully parallel
```

**Examples**:

| Task A | Task B | Decision | Reasoning |
|--------|--------|----------|-----------|
| Literature review on oxygen transport | Literature review on ammonia clearance | **Parallel** (if different agents) | Independent topics, outputs don't depend on each other |
| Parameter extraction | Back-of-envelope calculation | **Sequential** | Calculation needs parameter values |
| Back-of-envelope calculation | Detailed model | **Sequential** | Detailed model should only run if back-of-envelope shows feasibility |
| Detailed oxygen model | Detailed ammonia model | **Parallel** (if different agents or independent parameters) | Independent models (unless they share parameters) |
| Researcher writing draft | Editor polishing draft | **Sequential** | Editor needs completed draft |
| Researcher writing Section 1 | Researcher writing Section 2 | **Sequential** (same agent) | Can't parallelize single agent's work |

---

## Handoff Best Practices

### Handoff Ceremony (Minimum Viable)

Every handoff must include:

1. **Deliverable** (the file/output)
2. **Context** (why this matters, how it fits into larger project)
3. **Focus areas** (what recipient should check or emphasize)
4. **Known gaps** (what's incomplete or uncertain)

**Example handoff** (Calculator → Synthesizer):

```
DELIVERABLE:
- models/oxygen-transport/back-of-envelope-oxygen-feasibility.md
- models/oxygen-transport/spatial-oxygen-profile-model.md

CONTEXT:
We're assessing feasibility of 10-billion-cell BAL device. Oxygen delivery is the most
critical constraint (cells consume O₂ 10x faster than nutrients diffuse). These models
determine if hollow fiber architecture can support target cell mass.

FOCUS AREAS:
- Back-of-envelope conclusion: 12 m² membrane area, 48k fibers → FEASIBLE but large
- Detailed model: Sensitivity analysis shows cell density and fiber spacing are critical parameters
- Key finding: O₂ concentration drops to zero at 150 μm from fiber (limits cell packing density)

KNOWN GAPS:
- Detailed model is 1D (assumes radial symmetry); 2D model would capture end effects but
  we decided 1D sufficient for architecture decision (discussed with user, approved Jan 31)
- Flow rate optimization not included (defer to V2 if this architecture selected)

WHAT YOU SHOULD DO:
Integrate these calculations with the hollow fiber literature review (docs/literature/hollow-fiber-membranes/)
to recommend oxygen delivery architecture. Specifically:
- Cite the 12 m² / 48k fiber numbers as design requirement
- Connect 150 μm limit to literature's max cell packing distances
- Recommend fiber configuration (bundle size, cartridge design) based on these constraints
```

**Why this format works**:
- **Deliverable**: Recipient knows exactly what files to read
- **Context**: Recipient understands why this work matters (not just "here's a file")
- **Focus areas**: Recipient knows what's important (saves re-reading to figure out key points)
- **Known gaps**: Recipient doesn't waste time on out-of-scope questions or assume completeness
- **What you should do**: Recipient knows next steps (prevents "what do I do with this?" confusion)

---

### Handoff Acknowledgment Pattern

Recipient should confirm handoff received and understood:

**Synthesizer response** (to above handoff):
```
Handoff received. I understand:
- Use 12 m² / 48k fibers as design requirement
- Connect 150 μm O₂ diffusion limit to literature
- Recommend specific fiber configuration

I'll integrate with hollow fiber literature review and produce synthesis document at
docs/reports/analysis-oxygen-delivery-system-design.md. Estimated completion: 4 hours (Feb 3).

Clarifying question: Should I include cost implications (cartridge size → manufacturing complexity)
or defer to Economist agent?
```

**Why this works**:
- Confirms understanding (prevents misinterpretation)
- States deliverable location and timeline (sets expectations)
- Asks clarifying questions BEFORE starting work (prevents rework)

---

## Blocker Resolution Playbook

### Type 1: Missing Information

**Symptom**: Agent can't find parameter, data unavailable, literature search yields nothing

**Diagnosis questions**:
- Have we searched exhaustively? (PubMed, Google Scholar, citation chaining)
- Are we using correct keywords? (biological vs. clinical terminology)
- Is information actually published, or do we need experimental data?

**Resolution paths**:
1. **Broaden search**: Try related terms, older literature, adjacent fields
2. **Use proxies**: If exact value unknown, can we bound it or use related parameter?
3. **Escalate to user**: Ask if they have unpublished data or alternate source
4. **Document gap and proceed with range**: "Literature reports 5-15, we'll use 10 ± 50%"

**When to escalate**: After 30 minutes of unsuccessful search, use AskUserQuestion

---

### Type 2: Technical Complexity (Agent Stuck on Hard Problem)

**Symptom**: Calculator can't solve equation, model won't converge, parameter sweep taking forever

**Diagnosis questions**:
- Is this actually required, or can we simplify? (Do we need 2D model or is 1D sufficient?)
- Are we solving the right problem? (Maybe back-of-envelope already answers question)
- Is this a tool limitation or fundamental complexity?

**Resolution paths**:
1. **Simplify**: Reduce dimensions (3D → 2D → 1D), fewer parameters, analytical instead of numerical
2. **Break into sub-problems**: Solve pieces separately, integrate later
3. **Defer**: If not critical path, mark as "future work" and proceed with simpler approach
4. **Consult user**: If technical decision has strategic implications (accuracy vs. speed), escalate

**When to escalate**: If simplification changes conclusion or accuracy significantly, ask user if acceptable

---

### Type 3: Ambiguous Requirements

**Symptom**: Agent asks "Should I include X?", "How detailed should this be?", "Which approach?"

**Diagnosis questions**:
- Did task description specify scope and deliverable format?
- Are success criteria clear?
- Is this actually ambiguous, or is agent overthinking?

**Resolution paths**:
1. **Clarify scope**: Add "In scope" and "Out of scope" sections to task description
2. **Define deliverable format**: "2-page summary with parameter table, max 10 citations"
3. **Set decision criteria**: "If finding takes >1 hour, document and move on"
4. **Escalate if strategic**: If ambiguity reflects unclear user priorities, use AskUserQuestion

**When to escalate**: If decision affects project direction or resource allocation significantly

---

### Type 4: Time Constraint (Not Enough Time)

**Symptom**: Task estimated 4 hours, agent at 80% after 6 hours, deadline approaching

**Diagnosis questions**:
- Is task actually larger than estimated? (scope creep or bad estimate?)
- Is agent spinning or making steady progress?
- Can we narrow scope and still meet core objective?

**Resolution paths**:
1. **Narrow scope aggressively**: Cut nice-to-haves, reduce parameters, shorten document
2. **Trade breadth for depth**: Answer one question well instead of three poorly
3. **Accept "good enough"**: Does it need to be perfect, or is 80% sufficient for decision?
4. **Defer non-critical elements**: Move to "future work" section, complete core deliverable
5. **Extend deadline**: If scope reduction compromises core objective, negotiate timeline with user

**When to escalate**: If deadline at risk despite narrowing scope, notify user 1+ week early with options

---

## Agent Timeout Response Flowchart

```
Agent working >30 min without progress update
│
├─ Check progress file: Is agent making progress?
│   ├─ YES (recent updates, milestones advancing) → Continue monitoring, check again in 30 min
│   └─ NO (no updates, or repeated failed attempts) → INTERVENTION NEEDED
│       │
│       ├─ Diagnose issue:
│       │   ├─ Overly broad scope → Narrow focus
│       │   ├─ Ambiguous requirements → Clarify scope and deliverable
│       │   ├─ Missing information → Escalate to user or use proxy
│       │   └─ Technical complexity → Simplify or break into sub-tasks
│       │
│       └─ Execute intervention:
│           ├─ Update task description with narrower scope
│           ├─ Provide additional context or decision criteria
│           └─ If user input needed → AskUserQuestion with options
│
└─ Agent completes or reaches new milestone → Update work plan, check for handoffs
```

**Key principle**: Intervene early (30 min) to prevent waste, but don't micromanage steady progress.

---

## Coordination Anti-Patterns (What NOT to Do)

### ❌ "Fire and Forget" (Assign task, never check progress)
**Problem**: Agent spins for hours, blocker undetected until too late
**Fix**: Progress monitoring every 30-60 min for tasks >1 hour

### ❌ Implicit Dependencies ("Task B depends on Task A")
**Problem**: Agent doesn't know WHEN to start Task B (after A starts? After A completes? After PM reviews A?)
**Fix**: Explicit handoff timing: "Task B starts AFTER Task A deliverable reviewed and approved"

### ❌ No Handoff Context (Just pass file)
**Problem**: Recipient doesn't know why file matters, what to focus on, what's incomplete
**Fix**: Use handoff ceremony format (deliverable + context + focus + gaps)

### ❌ Over-Planning (2-hour Gantt chart for 3-hour task)
**Problem**: Planning overhead exceeds task duration; analysis paralysis
**Fix**: "Plan enough to start" principle—simple checklist for short tasks, detailed plan for multi-day work

### ❌ Scope Creep Tolerance ("Oh, they found something interesting, let them explore")
**Problem**: 1-hour task becomes 4 hours, timeline slips, resources overcommitted
**Fix**: Enforce scope boundaries ruthlessly; tangents are fine IF time-boxed

### ❌ No Decision Escalation ("I'll just decide this myself")
**Problem**: PM makes strategic call above authority level, user surprised by direction
**Fix**: Use Decision Escalation Framework—Major/Medium decisions always escalate

---

## Quick Reference: Common Handoffs

| From | To | Deliverable | Focus Areas | Estimated Duration |
|------|-----|-------------|-------------|-------------------|
| Researcher | Fact-Checker | Literature review draft | Citation accuracy, quote verification | 30-60 min |
| Researcher | Synthesizer | Multiple paper notes | Parameter convergence, contradictions | 3-5 hours synthesis |
| Calculator | Synthesizer | Model results | Key conclusions, parameter sensitivities | 3-4 hours synthesis |
| Synthesizer | Devil's Advocate | Synthesis document | Assumption validity, gap identification | 1-2 hours review |
| Devil's Advocate | Editor | Reviewed document | Prose clarity, CLAUDE.md style | 30-60 min per 1000 words |
| Editor | archive-workflow | Polished document | Project organization, naming conventions | 15-30 min |

---

## References

- **Handoff Protocols**: Based on surgical team coordination literature
- **Blocker Resolution**: Adapted from Agile "impediment removal" practices
- **Coordination Patterns**: Inspired by software engineering workflows (CI/CD pipelines, code review)
