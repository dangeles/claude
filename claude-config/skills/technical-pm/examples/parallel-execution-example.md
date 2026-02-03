# Parallel Execution Example

## Scenario

User requests a comprehensive analysis requiring both literature research and calculations on independent topics.

## User Goal

> "I need to understand hepatocyte oxygenation requirements from the literature AND calculate the oxygen delivery capacity of our bioreactor design. These will inform whether we need to redesign the oxygenation system."

## Step 1: Dependency Detection

technical-pm analyzes the goal:

**Task A**: Research hepatocyte oxygenation requirements (researcher)
**Task B**: Calculate bioreactor oxygen delivery capacity (calculator)

**Dependency analysis**:
- Explicit markers: None ("AND" connects tasks, no "based on", "using")
- Shared terms: "oxygenation", "oxygen" (2 terms, at threshold)
- Natural dependencies: researcher and calculator have no inherent dependency
- Duration estimate: researcher (45 min) + calculator (20 min) = 65 min sequential vs ~45 min parallel

**Decision**: PARALLEL (borderline - confirm with user for first-time)

## Step 2: User Confirmation

```
I plan to run these tasks in PARALLEL:

Task A: Researcher - review hepatocyte oxygenation literature
Task B: Calculator - compute bioreactor oxygen delivery capacity

These appear independent because:
- No explicit dependency markers
- 2 shared terms (at threshold)

However, parallel execution means:
- Calculator will NOT have access to literature values
- If calculation needs literature data (e.g., specific oxygen consumption rates),
  calculator will use general assumptions instead

Confirm parallel execution? (Y/N/Help me decide)
```

**User response**: "Y - the calculator can use standard textbook values for now. If we need more precision, we'll iterate."

## Step 3: Parallel Execution Setup

### Generate Batch ID
```
batch_id: batch-20260203-1430
```

### Create Output Directories
```bash
mkdir -p scratchpad/researcher/batch-20260203-1430/
mkdir -p scratchpad/calculator/batch-20260203-1430/
```

### Record in Workflow State
```yaml
parallel_batches:
  - batch_id: "batch-20260203-1430"
    workflow_id: "workflow-hepatocyte-analysis"
    tasks:
      - task_id: "batch-20260203-1430-researcher"
        skill: researcher
        status: launching
        output_dir: "scratchpad/researcher/batch-20260203-1430/"
        goal: "Review hepatocyte oxygenation requirements from literature"
      - task_id: "batch-20260203-1430-calculator"
        skill: calculator
        status: launching
        output_dir: "scratchpad/calculator/batch-20260203-1430/"
        goal: "Calculate oxygen delivery capacity of bioreactor design"
    started_at: "2026-02-03T14:30:00Z"
    status: in_progress
```

## Step 4: Launch Parallel Tasks

### Task 1: Researcher

**Template substitution**:
- {OUTPUT_DIR} -> scratchpad/researcher/batch-20260203-1430
- {TASK_ID} -> batch-20260203-1430-researcher
- {GOAL} -> "Review hepatocyte oxygenation requirements from literature"

**Invocation**:
```
Task(general-purpose, "[substituted researcher template]")
```

### Task 2: Calculator

**Template substitution**:
- {OUTPUT_DIR} -> scratchpad/calculator/batch-20260203-1430
- {TASK_ID} -> batch-20260203-1430-calculator
- {GOAL} -> "Calculate oxygen delivery capacity of bioreactor design"

**Invocation**:
```
Task(general-purpose, "[substituted calculator template]")
```

## Step 5: Monitor and Wait

**Progress updates**:
- 14:35 - Both tasks in progress
- 14:45 - Researcher: checkpoint.md updated ("Screening abstracts, 12 papers found")
- 14:55 - Calculator: checkpoint.md updated ("Completed mass transfer calculation")
- 15:05 - Calculator: COMPLETE (output.md exists, 520 chars)
- 15:15 - Researcher: checkpoint.md updated ("Reading Jiang 2024 review")
- 15:20 - Researcher: COMPLETE (output.md exists, 1450 chars)

**Elapsed time**: 50 minutes (vs estimated 65 minutes sequential)

## Step 6: Quality Validation

### Researcher Output Validation

**File**: scratchpad/researcher/batch-20260203-1430/output.md

**Quality checklist**:
- [x] Output exists at expected location
- [x] Output >= 500 characters (1450 chars)
- [x] Contains citations (5 found: [1], [2], [3], [4], [5])
- [x] Has "Executive Summary" section
- [x] No placeholder text
- [x] Template sentinel echoed ("TEMPLATE_RECEIVED: researcher_template_v1_complete")

**Result**: PASS

### Calculator Output Validation

**File**: scratchpad/calculator/batch-20260203-1430/output.md

**Quality checklist**:
- [x] Output exists at expected location
- [x] Output >= 100 characters (520 chars)
- [x] Contains numeric result with units ("Oxygen delivery capacity = 2.4 mL O2/min")
- [x] Shows work ("Step 1: Calculate membrane surface area...")
- [x] Result is plausible (2.4 mL/min is reasonable for small bioreactor)
- [x] Template sentinel echoed ("TEMPLATE_RECEIVED: calculator_template_v1_complete")

**Result**: PASS

## Step 7: Cross-Check for Contradictions

**Analysis**:
- Researcher output mentions: "Hepatocyte oxygen consumption: 0.4-0.8 nmol/s per 10^6 cells (Jiang 2024)"
- Calculator output uses: "Assuming oxygen consumption rate of 0.5 nmol/s per 10^6 cells (standard textbook value)"

**Assessment**: Values align within expected range. No contradiction detected.

## Step 8: Output Aggregation

```markdown
## Aggregated Parallel Outputs

**Batch ID**: batch-20260203-1430
**Tasks Completed**: 2/2
**Duration**: 50 minutes (parallel)

### Output 1: Researcher
- File: scratchpad/researcher/batch-20260203-1430/output.md
- Summary: Comprehensive review of hepatocyte oxygenation found oxygen consumption
  rates of 0.4-0.8 nmol/s per 10^6 cells, with critical thresholds below 0.3 nmol/s
  causing viability issues. Key factors: perfusion rate, membrane permeability.
- Key findings:
  - Oxygen consumption varies by hepatocyte source (primary > immortalized)
  - Hypoxia triggers metabolic shift at <50% saturation
  - Review of 5 key papers on bioreactor oxygenation
- Known gaps: Limited data on long-term (>7 day) cultures

### Output 2: Calculator
- File: scratchpad/calculator/batch-20260203-1430/output.md
- Summary: Bioreactor oxygen delivery capacity calculated at 2.4 mL O2/min
  based on membrane area (50 cm^2) and mass transfer coefficient.
- Key results:
  - Delivery capacity: 2.4 mL O2/min
  - Required for 10^7 cells: ~1.8 mL O2/min
  - Safety margin: 33% above requirement
- Assumptions: 0.5 nmol/s per 10^6 cells, standard PDMS membrane

### Detected Conflicts
- None. Calculator assumption (0.5 nmol/s) falls within researcher range (0.4-0.8 nmol/s).

### Integration Notes
- Current bioreactor design has adequate safety margin (33%)
- If using primary hepatocytes (higher consumption), margin reduces to ~15%
- Consider sensitivity analysis for long-term culture scenarios
```

## Step 9: Synthesis Decision

With both outputs validated and aggregated, technical-pm can:

1. **Deliver directly**: If user only wanted research + calculation, present aggregated summary
2. **Continue to synthesis**: If user wants integrated document, invoke synthesizer with aggregated outputs

**In this case**: User wanted to "inform whether we need to redesign"

**Decision**: Present aggregated findings with recommendation:

```
## Parallel Analysis Complete

I've completed both the literature review and capacity calculation in parallel.

### Key Finding
Your current bioreactor design has a 33% safety margin for oxygen delivery,
which is adequate for standard conditions. However:

- If using primary hepatocytes (vs immortalized), margin drops to ~15%
- Long-term cultures (>7 days) have limited data and may have higher demands

### Recommendation
**No immediate redesign needed** for standard use case. Consider redesign if:
- Scaling to primary hepatocyte cultures
- Extended culture duration (>7 days)
- Higher cell density requirements

### Next Steps (your choice)
(A) Proceed with current design
(B) Run sensitivity analysis for primary hepatocytes
(C) Get full synthesis document combining both analyses
(D) Review individual outputs in detail
```

## Summary

| Metric | Value |
|--------|-------|
| Tasks parallelized | 2 (researcher, calculator) |
| Parallel duration | 50 minutes |
| Sequential estimate | 65 minutes |
| Time saved | 15 minutes (23%) |
| Quality gates passed | 2/2 |
| Contradictions detected | 0 |
| User decisions required | 1 (confirm parallel) |

## Lessons Learned

1. **First-time confirmation** valuable - user provided context that calculator could use assumptions
2. **Cross-check** found alignment, increasing confidence in parallel approach
3. **Checkpoints** enabled progress monitoring without interrupting tasks
4. **Template sentinels** confirmed no truncation in either task
