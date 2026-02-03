# Parallel Execution Reference

## Overview

This document defines the protocol for running researcher, calculator, and synthesizer as parallel Task agents with embedded templates. Use this for long-running, independent tasks to achieve 2-3x speedup.

## When to Use Parallel Execution

### Eligible Skills

Only these skills are eligible for parallel Task execution:

| Skill | Estimated Duration | Why Parallel? |
|-------|-------------------|---------------|
| researcher | 30-60 min | Long-running, produces self-contained output |
| calculator | 5-30 min | Long-running, produces self-contained output |
| synthesizer | 15-30 min | Can run parallel IF inputs are independent |

All other skills (devils-advocate, fact-checker, editor, archivist) remain sequential via Skill tool.

### Decision Criteria

Use parallel execution when:
- 2+ eligible skills are identified
- No data dependency between tasks (see dependency-detection.md)
- Combined estimated duration > 30 minutes
- Tasks operate on different topics/domains

Use sequential execution when:
- Tasks have explicit dependency ("based on", "using results from")
- Quality is more important than speed
- Tasks share more than 2 technical terms (ambiguous)
- User explicitly requests sequential (--sequential flag)

## Orchestrator Protocol

When launching parallel tasks, technical-pm acts as Process Orchestrator with these phases:

### Phase 1: Dispatch

1. **Generate batch_id**: Unique identifier for parallel execution group
   - Format: `batch-{timestamp}-{random4}`
   - Example: `batch-20260203-a1b2`

2. **Validate parallel-safe**: Run dependency detection for all task pairs
   - If any pair returns 'sequential': Fall back to sequential execution
   - If any pair returns 'ask_user': Confirm with user before parallelizing

3. **Assign task IDs**: Each task gets unique ID
   - Format: `{batch_id}-{skill_name}`
   - Example: `batch-20260203-a1b2-researcher`

4. **Create output directories**:
   ```bash
   mkdir -p scratchpad/{skill_name}/{batch_id}/
   ```

5. **Launch Tasks**: For each parallel task:
   - Load template from task-templates.md
   - Substitute variables: {OUTPUT_DIR}, {TASK_ID}, {GOAL}
   - Invoke: `Task(general-purpose, substituted_template)`

6. **Record batch in workflow state**:
   ```yaml
   parallel_batches:
     - batch_id: "batch-20260203-a1b2"
       workflow_id: "workflow-xyz789"
       tasks:
         - task_id: "batch-20260203-a1b2-researcher"
           skill: researcher
           status: in_progress
           output_dir: "scratchpad/researcher/batch-20260203-a1b2/"
         - task_id: "batch-20260203-a1b2-calculator"
           skill: calculator
           status: in_progress
           output_dir: "scratchpad/calculator/batch-20260203-a1b2/"
       started_at: "2026-02-03T14:00:00Z"
       status: in_progress
   ```

### Phase 2: Monitor

Wait for Task completion. Since Tasks run asynchronously:

1. **Check completion**: Each Task completes when output file exists
2. **Detect timeout**: If task exceeds timeout threshold (see below)
3. **Log progress**: "Parallel batch {batch_id}: N/M tasks complete"

### Phase 3: Collect

When all Tasks complete (or timeout):

1. **Gather outputs**: Read output file from each task's output_dir
2. **Validate each output**: Apply quality gate (see Quality Gates section)
3. **Aggregate status**: Count passed/failed/timeout

### Phase 4: Route

Based on aggregate status:

| Outcome | Action |
|---------|--------|
| All pass validation | Proceed to synthesis |
| Some fail validation | Apply error handling, present user options |
| All fail | Trigger catastrophic failure protocol |

## Resource Isolation Protocol

### Output Path Isolation

Each parallel task MUST write to unique paths to prevent collision.

**Path template**: `scratchpad/{skill_name}/{batch_id}/`

**Example paths**:
```
scratchpad/researcher/batch-20260203-a1b2/output.md
scratchpad/researcher/batch-20260203-a1b2/checkpoint.md
scratchpad/calculator/batch-20260203-a1b2/output.md
scratchpad/calculator/batch-20260203-a1b2/checkpoint.md
```

### Pre-Launch Collision Check

Before launching parallel batch, verify no path collisions:

```python
def validate_output_paths(parallel_tasks):
    paths = [task.output_dir for task in parallel_tasks]
    if len(paths) != len(set(paths)):
        raise ValueError(f"Output path collision detected: {paths}")
    for path in paths:
        if os.path.exists(path) and os.listdir(path):
            raise ValueError(f"Output directory not empty: {path}")
```

### Template Variable Substitution

Templates use placeholders, NOT hardcoded paths:

**In template**:
```
Write your output to: {OUTPUT_DIR}/output.md
Write checkpoints to: {OUTPUT_DIR}/checkpoint.md
```

**technical-pm substitutes before launch**:
```
{OUTPUT_DIR} -> scratchpad/researcher/batch-20260203-a1b2
{TASK_ID} -> batch-20260203-a1b2-researcher
{GOAL} -> [user's specific goal for this task]
```

**CRITICAL**: Never allow Task subagent to choose its own output path. Always substitute before launch.

## Timeout Configuration

### Per-Task Timeouts

| Skill | Default Timeout | Extended Timeout (retry) |
|-------|----------------|--------------------------|
| researcher | 60 min | 90 min |
| calculator | 30 min | 45 min |
| synthesizer | 45 min | 60 min |

### Batch Timeout

Maximum time for entire parallel batch: **2 hours**

Allows for retries within batch.

### Timeout Detection

Check after each Task completes:
1. Has any task exceeded its timeout? (no completion after timeout period)
2. If batch timeout exceeded: Treat remaining tasks as timeout

### On Timeout

When task timeout detected:

1. Mark task as TIMEOUT in workflow state
2. Check for partial output:
   - Look for `{OUTPUT_DIR}/checkpoint.md`
   - Look for `{OUTPUT_DIR}/output.md` (partial)
3. Present user options:

```
Task timeout: {task_id}

Task: {description}
Started: {start_time} ({elapsed} minutes ago)
Partial output: {yes/no, location if yes}

Options:
(A) Retry with extended timeout (+50%)
(B) Skip this task, continue with completed outputs
(C) Review partial output before deciding
(D) Abort parallel execution
```

## Quality Gates

### Gate 1: Pre-Launch Gate

**Location**: Before Task tool invocation

**Criteria**:
- [ ] Dependency detection confirms tasks are parallel-safe
- [ ] All required inputs available
- [ ] Output paths unique (no collision)
- [ ] Directories created
- [ ] Workflow state shows EXECUTING status

**On Failure**: Do not launch. Return to dependency detection.

### Gate 2: Task Completion Gate

**Location**: When each parallel Task completes

**Criteria by skill**:

**Researcher**:
- [ ] Output file exists at expected location
- [ ] Output >= 500 characters
- [ ] Contains at least 1 citation (regex: `\[\d+\]` or `(Author, \d{4})`)
- [ ] Has "Executive Summary" or "Key Findings" section
- [ ] No placeholder text ("TODO", "[INSERT HERE]")
- [ ] Template integrity sentinel echoed

**Calculator**:
- [ ] Output file exists at expected location
- [ ] Output >= 100 characters
- [ ] Contains numeric result with units
- [ ] Shows work or methodology ("Calculation", "Result", "=")
- [ ] Result is plausible (not NaN, not Infinity)
- [ ] Template integrity sentinel echoed

**Synthesizer**:
- [ ] Output file exists at expected location
- [ ] Output >= 300 characters
- [ ] Has "Executive Summary" section
- [ ] References multiple source documents
- [ ] No placeholder text
- [ ] Template integrity sentinel echoed

**On Failure**: Log specific criteria failed. Apply error handling.

### Gate 3: Batch Synthesis Gate

**Location**: Before synthesis phase

**Criteria**:
- [ ] At least 50% of tasks passed Gate 2 (or user override)
- [ ] No output file conflicts
- [ ] Combined context within capacity

**On Failure**: Present user options before proceeding.

## Error Handling

### Quality Gate Failure

When task output fails quality validation:

```
Quality validation failed for: {skill}

Criteria failed:
- [ ] {failed_criterion_1}: FAILED ({details})
- [ ] {failed_criterion_2}: FAILED ({details})

Output preview (first 200 chars):
"{preview_text}..."

Options:
(A) Accept anyway - proceed with this output (risky: downstream quality may suffer)
(B) Retry with explicit instruction: "{suggested_fix}"
(C) Review full output before deciding
(D) Skip this task's output, continue with others only
```

### Quality Threshold Levels

For flexibility, allow threshold adjustment:

| Level | Criteria | When to Use |
|-------|----------|-------------|
| HIGH (default) | All criteria must pass | Production, final deliverables |
| MEDIUM | Major criteria only (length, no placeholders) | Drafts, iteration |
| LOW | Output exists and >100 chars | Emergency, exploration |

### Catastrophic Failure Protocol

When ALL parallel tasks fail (no successful completions):

**Detection**:
- Batch complete with 0 successful tasks
- Or: All tasks timeout within 10 minutes of each other (correlated failure)

**Response**:

1. Do NOT auto-retry (likely systemic issue)
2. Check for ANY partial outputs (checkpoints)
3. Present comprehensive failure report:

```
PARALLEL EXECUTION FAILED: All {N} tasks failed

Failures:
1. researcher: {failure_reason}
2. calculator: {failure_reason}

Possible causes:
- System overload (correlated timeouts)
- Context window exhaustion
- Invalid templates (check template integrity)

Options:
(A) Retry SEQUENTIALLY (slower but more stable)
(B) Retry with reduced scope (specify narrowing)
(C) Abort workflow, preserve any partial outputs
(D) Debug: Show full error details
```

### Sequential Fallback

If user selects (A) after catastrophic failure:
- Switch to sequential execution
- Use Skill tool instead of Task tool
- Direct skill invocation with full skill context
- Trade parallelism for reliability

## Template Integrity Validation

### Sentinel Protocol

Every template includes at its END:

```
---
TEMPLATE_INTEGRITY_SENTINEL: {skill_name}_template_v1_complete
---
```

**In template instructions** (near start):

```
IMPORTANT: At the START of your response, confirm you received the complete
template by echoing: "TEMPLATE_RECEIVED: {skill_name}_template_v1_complete"
```

### Truncation Detection

After Task completes, check output:

1. Search for "TEMPLATE_RECEIVED:" in first 500 chars of output
2. If found: Extract sentinel value, compare to expected
3. If missing or different:
   - Flag: "Possible template truncation detected"
   - Options: (A) Review output anyway (B) Retry with shorter template (C) Abort

### Mitigation: Modular Templates

If truncation confirmed, restructure templates:
- Core instructions: 60 lines max (essential workflow)
- Quality checklist: 30 lines (at TOP, not bottom)
- Output format: 30 lines
- Personality: 20 lines (can be shortened)

## Parallel Cancellation Protocol

On Ctrl+C during parallel execution:

1. **Record status** per task: complete / in_progress / not_started
2. **Preserve completed outputs**: Do not delete successful task outputs
3. **Save workflow state** with parallel_status:

```yaml
parallel_execution:
  batch_id: "batch-20260203-a1b2"
  status: interrupted
  tasks:
    - task_id: "batch-20260203-a1b2-researcher"
      status: complete
      output: "scratchpad/researcher/batch-20260203-a1b2/output.md"
    - task_id: "batch-20260203-a1b2-calculator"
      status: interrupted
      partial: "scratchpad/calculator/batch-20260203-a1b2/checkpoint.md"
```

4. **Notify user**:

```
Parallel execution paused.

Completed:
- [x] researcher: scratchpad/researcher/batch-20260203-a1b2/output.md

Interrupted (partial may exist):
- [ ] calculator: scratchpad/calculator/batch-20260203-a1b2/checkpoint.md

To resume: Run technical-pm with resume=true
Options on resume:
(A) Continue incomplete tasks from checkpoint
(B) Restart incomplete tasks from beginning
(C) Proceed with completed outputs only
```

## Enhanced Dependency Detection

For parallel-eligible skills, add semantic dependency markers beyond keywords.

### Semantic Dependency Markers

1. **Value reference**: Does task B's description mention "values", "data", "results" from task A's domain?
2. **Calculation inputs**: Does calculator task reference any research topic?
3. **Synthesis requirements**: Does synthesizer explicitly list multiple inputs?

### First-Time Parallel Confirmation

When parallelizing for the FIRST time in a session, ALWAYS confirm:

```
I plan to run these tasks in PARALLEL:

Task A: Researcher - {description}
Task B: Calculator - {description}

These appear independent because:
- No explicit dependency markers ("based on", "using results")
- Limited shared terms ({count} shared terms)

However, parallel execution means:
- Calculator will NOT have access to literature values
- If calculation needs literature data, results may be incorrect

Confirm parallel execution? (Y/N/Help me decide)
```

### Post-Completion Cross-Check

Before synthesis, verify parallel outputs don't contradict:

1. Extract key numeric values from all parallel outputs
2. Identify shared terms/concepts mentioned in multiple outputs
3. If discrepancy found:

```
Potential contradiction detected:

researcher output: "{quote_1}"
calculator output: "{quote_2}"

These values differ by {magnitude}. Options:
(A) Proceed with researcher value (literature source)
(B) Proceed with calculator assumption
(C) Review both outputs to resolve
(D) Re-run calculator with researcher data
```

## Output Aggregation Protocol

### Step 1: Collect Outputs

```python
outputs = []
for task in completed_tasks:
    output = read_file(task.output_path)
    outputs.append({
        'task_id': task.task_id,
        'skill': task.skill,
        'content': output,
        'quality_passed': task.quality_passed
    })
```

### Step 2: Order for Synthesis

Order outputs for synthesis based on:
1. Natural skill order (researcher before calculator results typically)
2. Quality score (highest quality first)
3. Completion time (earlier first, as tie-breaker)

### Step 3: Create Aggregation Summary

For synthesis phase input:

```markdown
## Aggregated Parallel Outputs

**Batch ID**: {batch_id}
**Tasks Completed**: {completed}/{total}

### Output 1: Researcher
- File: {output_path}
- Summary: {50-word summary}
- Key findings: {bullet points}
- Known gaps: {from output}

### Output 2: Calculator
- File: {output_path}
- Summary: {50-word summary}
- Key results: {bullet points}
- Assumptions: {from output}

### Detected Conflicts
- {None or list conflicts}

### Integration Notes
- {Any notes for synthesizer}
```

## Examples

### Example 1: Successful Parallel Execution

**Goal**: "Research hepatocyte oxygenation literature AND calculate bioreactor oxygen delivery capacity"

**Dependency detection**:
- No explicit markers
- Shared terms: "oxygen" (1 term, < 2 threshold)
- Different domains: biology vs engineering
- Decision: PARALLEL

**Execution**:
1. Generate batch_id: batch-20260203-a1b2
2. Create directories:
   - scratchpad/researcher/batch-20260203-a1b2/
   - scratchpad/calculator/batch-20260203-a1b2/
3. Launch Tasks:
   - Task(general-purpose, researcher_template with {GOAL}="hepatocyte oxygenation literature")
   - Task(general-purpose, calculator_template with {GOAL}="bioreactor oxygen delivery capacity")
4. Wait for completion (estimated 60 min for researcher, 30 min for calculator)
5. Validate outputs:
   - Researcher: PASS (650 chars, 3 citations, has Key Findings)
   - Calculator: PASS (450 chars, numeric result with mL/min units)
6. Aggregate outputs
7. Proceed to synthesis (if needed)

### Example 2: Partial Failure Recovery

**Goal**: "Research membrane materials AND calculate manufacturing costs AND analyze regulatory requirements"

**Execution**:
1. Launch 3 parallel Tasks
2. Results:
   - Researcher: PASS
   - Calculator: TIMEOUT (no output after 45 min)
   - (Third task omitted for brevity - would be separate skill)
3. User chooses: (B) Skip calculator, continue with researcher output
4. Synthesis proceeds with partial results
5. Note in synthesis: "Manufacturing costs not calculated; may require follow-up"

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-03 | Initial parallel execution protocol |

---

**Last updated**: 2026-02-03
**Maintained by**: technical-pm skill
