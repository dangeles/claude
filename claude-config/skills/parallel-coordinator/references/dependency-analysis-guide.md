# Dependency Analysis Guide for Parallel Execution

## Purpose

This guide provides a systematic framework for analyzing task dependencies to determine whether tasks can be safely executed in parallel. Proper dependency analysis is critical for avoiding race conditions, data corruption, and incorrect results when coordinating parallel workflows.

## Core Dependency Types

### 1. Data Dependencies

**Definition**: Task B has a data dependency on Task A if Task B requires data produced by Task A.

**Categories**:

#### Read-After-Write (RAW) - True Dependency
- **Description**: Task B reads data that Task A writes
- **Example**: Task A generates a summary statistic; Task B creates a visualization using that statistic
- **Parallelization**: UNSAFE - Tasks must execute sequentially (A then B)

#### Write-After-Read (WAR) - Anti-Dependency
- **Description**: Task B writes to a location that Task A reads from
- **Example**: Task A reads configuration file; Task B updates the same configuration file
- **Parallelization**: UNSAFE if both operate on same resource - Use versioning or copies

#### Write-After-Write (WAW) - Output Dependency
- **Description**: Both Task A and Task B write to the same location
- **Example**: Task A writes results to output.txt; Task B also writes to output.txt
- **Parallelization**: UNSAFE - Results are non-deterministic, one will overwrite the other

#### Read-After-Read (RAR) - No Dependency
- **Description**: Both tasks only read from the same location
- **Example**: Task A reads config.json; Task B reads config.json
- **Parallelization**: SAFE - Multiple readers can access simultaneously

### 2. Control Dependencies

**Definition**: Task B has a control dependency on Task A if Task B's execution decision depends on Task A's outcome.

**Examples**:
- "If analysis shows error rate > 5%, then run debugging task"
- "Process dataset A first; if successful, process dataset B"
- "Run tests; if they pass, deploy"

**Parallelization**: UNSAFE - Task B cannot start until Task A completes and control decision is made

### 3. Resource Dependencies

**Definition**: Tasks depend on limited or exclusive resources.

**Categories**:

#### Exclusive Resources
- **Description**: Only one task can access the resource at a time
- **Examples**:
  - File locks
  - Database transactions requiring table locks
  - Hardware devices (specific GPU, instrument)
- **Parallelization**: UNSAFE without explicit coordination

#### Rate-Limited Resources
- **Description**: Resource has throughput or concurrency limits
- **Examples**:
  - API with rate limits (e.g., 100 requests/minute)
  - Database connection pool (e.g., max 10 concurrent connections)
  - Memory constraints (tasks combined exceed available RAM)
- **Parallelization**: CONDITIONAL - May need throttling or batching

#### Shared Read-Only Resources
- **Description**: Multiple tasks access the same resource without modification
- **Examples**:
  - Reading from same file
  - Querying same database table (no writes)
  - Using same API for independent queries
- **Parallelization**: SAFE - No conflicts from concurrent reads

### 4. Temporal Dependencies

**Definition**: Tasks have timing constraints or ordering requirements.

**Examples**:
- "Task A must complete before deadline X; Task B can start anytime"
- "Process data in chronological order"
- "Apply migrations in sequence"

**Parallelization**: Context-dependent - Analyze whether timing constraint creates true dependency

## Systematic Dependency Analysis Process

### Step 1: Enumerate Tasks

List all tasks explicitly with clear boundaries:

```
Task 1: [Clear description]
Task 2: [Clear description]
Task 3: [Clear description]
...
```

### Step 2: Identify Inputs and Outputs

For each task, document:

| Task | Inputs | Outputs | Resources Used |
|------|--------|---------|----------------|
| 1    | [List] | [List]  | [List]         |
| 2    | [List] | [List]  | [List]         |
| 3    | [List] | [List]  | [List]         |

### Step 3: Build Dependency Matrix

Create a matrix showing relationships between tasks:

|       | Task 1 | Task 2 | Task 3 |
|-------|--------|--------|--------|
| Task 1|   -    |   ?    |   ?    |
| Task 2|   ?    |   -    |   ?    |
| Task 3|   ?    |   ?    |   -    |

For each cell, mark:
- **D** (Dependent): Row task depends on column task
- **I** (Independent): No dependency
- **C** (Conflict): Tasks conflict on shared resource

### Step 4: Check for Circular Dependencies

Ensure no cycles exist:
- Task A depends on Task B
- Task B depends on Task C
- Task C depends on Task A (CIRCULAR - INVALID)

Circular dependencies indicate design problem and prevent both sequential and parallel execution.

### Step 5: Identify Independent Task Sets

From the dependency matrix, identify maximal sets of mutually independent tasks. These sets can execute in parallel.

**Example**:

If dependency matrix shows:
- Task 1: Independent of all others
- Task 2: Depends on Task 1
- Task 3: Independent of all others
- Task 4: Depends on Task 2

Then:
- **Parallel Set 1**: [Task 1, Task 3] can run concurrently
- **Sequential Phase**: After Set 1 completes
- **Parallel Set 2**: [Task 2, Task 4] cannot run together (Task 4 depends on Task 2)

## Practical Analysis Techniques

### Technique 1: Input-Output Intersection Test

```
For each pair of tasks (A, B):
  outputs_A = set(Task A outputs)
  inputs_B = set(Task B inputs)

  if outputs_A ∩ inputs_B != ∅:
    # Task B depends on Task A (RAW dependency)
    mark_dependent(B, A)

  if outputs_A ∩ outputs_B != ∅:
    # Tasks conflict on output (WAW dependency)
    mark_conflict(A, B)
```

### Technique 2: Resource Access Pattern Analysis

```
For each pair of tasks (A, B):
  resources_A = set(Task A resources)
  resources_B = set(Task B resources)
  shared = resources_A ∩ resources_B

  if shared == ∅:
    # No shared resources - safe to parallelize
    continue

  for resource in shared:
    if resource.is_mutable:
      if task_A_writes(resource) or task_B_writes(resource):
        # Shared mutable resource with writes - conflict
        mark_conflict(A, B)
      else:
        # Both read-only - safe
        continue
    else:
      # Immutable resource - safe
      continue
```

### Technique 3: Execution Trace Simulation

Mentally (or actually) simulate parallel execution:

1. Assume tasks start simultaneously
2. Walk through each task's steps
3. Identify any point where Task B would need data Task A hasn't produced yet
4. Identify any point where both tasks access same mutable resource
5. If any such point exists, tasks have dependency

## Common Dependency Patterns

### Pattern 1: Pipeline (Sequential Chain)

```
Task A → Task B → Task C → Task D
```

**Dependencies**: Each task depends on previous one
**Parallelization**: NONE - Must execute sequentially
**Alternative**: Look for opportunities to split pipeline into independent pipelines

### Pattern 2: Fan-Out (Independent Branches)

```
        → Task B
Task A → Task C
        → Task D
```

**Dependencies**: B, C, D all depend on A, but are independent of each other
**Parallelization**: After A completes, B, C, D can run in parallel
**Strategy**: Sequential then parallel

### Pattern 3: Fork-Join

```
Task A → Task B ↘
Task C → Task D → Task F
Task E ↗
```

**Dependencies**:
- B depends on A
- D depends on C
- F depends on B, D, E

**Parallelization**:
- Phase 1: A, C can run in parallel
- Phase 2: After Phase 1, B, D, E can run in parallel
- Phase 3: After Phase 2, F runs alone

**Strategy**: Multiple parallel phases

### Pattern 4: Fully Independent (Embarrassingly Parallel)

```
Task A
Task B
Task C
Task D
```

**Dependencies**: None
**Parallelization**: ALL tasks can run simultaneously
**Strategy**: Maximum parallelism

### Pattern 5: Shared Reader

```
        → Task B (reads X)
File X → Task C (reads X)
        → Task D (reads X)
```

**Dependencies**: All tasks read same file, but don't depend on each other
**Parallelization**: SAFE - All tasks can run in parallel
**Note**: Verify X is not modified during execution

## Red Flags: Signs of Hidden Dependencies

### Red Flag 1: Vague Task Descriptions

**Problem**: "Process data and generate report"
- Unclear what "process" means
- Unknown what "report" requires

**Solution**: Break down into specific steps:
- "Load data from CSV"
- "Calculate summary statistics"
- "Generate PDF report using statistics"

Now dependencies become clear: report generation depends on statistics.

### Red Flag 2: Shared Mutable State

**Problem**: Two tasks both "update configuration"
- Who updates first?
- What if updates conflict?
- Final state is non-deterministic

**Solution**:
- Make updates exclusive (sequential)
- Or use versioning (parallel with merge)
- Or partition configuration space (parallel on different sections)

### Red Flag 3: Implicit Ordering Assumptions

**Problem**: Task list appears ordered for a reason
- "Analyze data, clean data, visualize data"
- Order suggests potential dependencies

**Solution**: Verify if order is:
- **Necessary** (dependencies exist): Keep sequential
- **Arbitrary** (no dependencies): Can parallelize

### Red Flag 4: "Then" or "After" in Task Descriptions

**Problem**: "Fetch data, then analyze it"
- Explicit sequential language

**Solution**: If "then" implies dependency, tasks are sequential. If merely habitual phrasing, verify actual dependencies.

### Red Flag 5: Side Effects

**Problem**: Task has side effects not captured in output description
- "Analyze log file" (side effect: updates .last_read pointer)
- "Query database" (side effect: increments query counter)

**Solution**: Identify all side effects, treat them as outputs in dependency analysis

## Decision Tree for Parallelization

```
START
  │
  ├─ Do tasks share inputs?
  │    │
  │    ├─ YES → Are inputs read-only?
  │    │         │
  │    │         ├─ YES → Continue
  │    │         └─ NO → UNSAFE (abort parallelization)
  │    │
  │    └─ NO → Continue
  │
  ├─ Does any task's output serve as another's input?
  │    │
  │    ├─ YES → UNSAFE (tasks have data dependency)
  │    └─ NO → Continue
  │
  ├─ Do tasks write to same location?
  │    │
  │    ├─ YES → UNSAFE (output conflict)
  │    └─ NO → Continue
  │
  ├─ Does execution decision for one task depend on another's outcome?
  │    │
  │    ├─ YES → UNSAFE (control dependency)
  │    └─ NO → Continue
  │
  ├─ Do tasks require exclusive resources?
  │    │
  │    ├─ YES → UNSAFE (resource conflict)
  │    └─ NO → Continue
  │
  └─ SAFE to parallelize
```

## Verification Checklist

Before executing tasks in parallel, verify:

- [ ] Each task's inputs are available at start (no waiting for other tasks)
- [ ] No task reads data written by another task in this parallel set
- [ ] No two tasks write to the same file, database record, or output location
- [ ] All shared resources are read-only OR properly partitioned
- [ ] No task's execution is conditional on another task's outcome
- [ ] No circular dependencies exist
- [ ] Resource limits (memory, API rates, connections) can accommodate all parallel tasks
- [ ] Failure of one task won't corrupt state for other tasks
- [ ] Results can be aggregated after completion without ordering requirements

## Example Analysis

### Scenario

User requests: "Analyze file A.csv for outliers, calculate summary statistics for B.csv, and merge datasets C.csv and D.csv"

### Analysis

**Task Decomposition**:
1. Task 1: Analyze A.csv for outliers
2. Task 2: Calculate summary statistics for B.csv
3. Task 3: Merge C.csv and D.csv

**Inputs and Outputs**:

| Task | Inputs       | Outputs              | Resources Used    |
|------|--------------|----------------------|-------------------|
| 1    | A.csv        | outlier_report.txt   | A.csv (read)      |
| 2    | B.csv        | statistics.json      | B.csv (read)      |
| 3    | C.csv, D.csv | merged.csv           | C.csv, D.csv (read)|

**Dependency Matrix**:

|       | Task 1 | Task 2 | Task 3 |
|-------|--------|--------|--------|
| Task 1|   -    |   I    |   I    |
| Task 2|   I    |   -    |   I    |
| Task 3|   I    |   I    |   -    |

All cells show Independence (I).

**Resource Analysis**:
- No shared files (A, B, C, D all different)
- All resources read-only
- Different output files (no conflicts)

**Decision**: SAFE to parallelize all three tasks.

**Expected Speedup**: Near 3x (if tasks are similar duration and I/O-bound).

## Conclusion

Thorough dependency analysis is the foundation of safe and effective parallel execution. By systematically examining data flow, resource usage, control flow, and timing constraints, you can confidently identify truly independent tasks suitable for concurrent execution. When in doubt, conservative sequential execution is safer than risky parallel execution that could produce incorrect or non-deterministic results.
