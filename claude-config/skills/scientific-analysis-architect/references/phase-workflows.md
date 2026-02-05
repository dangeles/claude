# Phase Workflows

Detailed workflows for each phase of the scientific-analysis-architect skill.

## Phase 0: Initialization

**Duration**: 2-5 minutes
**Owner**: Orchestrator (main skill)
**Checkpoint**: Never (automatic)

### Workflow Steps

```
START Phase 0
    |
    v
[1] Check nbformat availability
    |
    +---> If missing: ERROR - abort with install instructions
    |
    v
[2] Check jupytext availability
    |
    +---> If missing: WARNING - proceed without Jupytext metadata
    |
    v
[3] Ask user for output directory
    |
    +---> Default: current working directory
    |
    v
[4] Validate output directory
    |
    +---> Check exists
    +---> Check writable (perform write test)
    +---> If fails: offer alternatives or ask for new path
    |
    v
[5] Create session directory
    |
    +---> Primary: {output_dir}/.scientific-analysis-session/
    +---> Fallback: /tmp/scientific-analysis-architect-session-{timestamp}/
    |
    v
[6] Initialize session-state.json
    |
    +---> Set status: "initialized"
    +---> Set current_phase: 0
    +---> Record config: output_directory
    |
    v
[7] Check for existing sessions
    |
    +---> If found and < 72 hours: offer resume
    |
    v
END Phase 0 -> Quality Gate 0
```

### Quality Gate 0 Criteria

- [ ] nbformat available (required)
- [ ] Session directory created and writable
- [ ] Output directory validated
- [ ] session-state.json initialized

---

## Phase 1: Birds-Eye Planning

**Duration**: ~12 minutes
**Owner**: research-architect (Sonnet 4.5)
**Timeout**: 15 minutes

### Workflow Steps

```
START Phase 1
    |
    v
[1] Load session state
    |
    v
[2] AskUserQuestion: Dataset and research goals
    |
    "Please describe your dataset and research goals.
     Include:
     - Data type (RNA-seq, proteomics, imaging, etc.)
     - Sample size and conditions
     - Main research questions"
    |
    v
[3] Analyze response for uncertainty signals
    |
    +---> Keywords: "not sure", "maybe", "could be", "uncertain"
    |
    +---> If HIGH uncertainty detected:
    |         |
    |         v
    |     [3a] Spawn brainstormers (parallel)
    |         |
    |         +---> Task: analysis-brainstormer
    |         |     "Suggest analysis types for: {dataset}"
    |         |
    |         +---> Task: method-brainstormer
    |         |     "Suggest methods for: {research_area}"
    |         |
    |         v
    |     [3b] Aggregate suggestions
    |         |
    |         v
    |     [3c] Present suggestions to user
    |         |
    |         AskUserQuestion: "Based on your description,
    |         here are analysis approaches to consider:
    |         {brainstormer_suggestions}
    |
    |         Which approaches interest you?"
    |
    v
[4] Generate research structure (3-7 chapters)
    |
    +---> Each chapter needs:
    |     - Title
    |     - Goal (atlas/hypothesis/mechanism)
    |     - List of analyses
    |     - Dependencies on other chapters
    |
    v
[5] Write research-structure.md
    |
    v
[6] Update session state
    |
    +---> Set current_phase: 1
    +---> Set completed_phases: [0, 1]
    +---> Record outputs.research_structure
    |
    v
END Phase 1 -> Quality Gate 1
```

### research-structure.md Format

```markdown
# Research Structure: {Project Title}

## Overview
{1-2 paragraph summary of research goals and approach}

## Dataset Description
{User-provided description, preserved verbatim}

## Chapters

### Chapter 1: {Chapter Title}
**Goal**: {atlas | hypothesis | mechanism}
**Analyses**:
1. {Analysis 1 name} - {brief description}
2. {Analysis 2 name} - {brief description}
**Dependencies**: None

### Chapter 2: {Chapter Title}
**Goal**: {atlas | hypothesis | mechanism}
**Analyses**:
1. {Analysis 1 name}
2. {Analysis 2 name}
**Dependencies**: Chapter 1 (requires normalized data)

...
```

### Quality Gate 1 Criteria

- [ ] research-structure.md exists
- [ ] Contains 3-7 chapters
- [ ] Each chapter has goal and at least 1 analysis
- [ ] Dependencies are valid (no circular, reference existing chapters)

---

## Phase 2: Subsection Planning

**Duration**: ~12 minutes for 3 chapters
**Owner**: analysis-planner (Sonnet 4.5)
**Timeout**: 20 minutes total

### Workflow Steps

```
START Phase 2
    |
    v
[1] Read research-structure.md
    |
    v
[2] For each chapter (sequential):
    |
    +---> [2a] Read chapter details
    |
    +---> [2b] Fan-out to expert panel (PARALLEL)
    |         |
    |         +---> Task: statistician-consultant
    |         |     "Review analyses for Chapter {N}: {analyses}"
    |         |     Timeout: 5 min
    |         |
    |         +---> Task: mathematician-consultant
    |         |     "Review algorithms for Chapter {N}: {analyses}"
    |         |     Timeout: 5 min
    |         |
    |         +---> Task: programmer-consultant
    |               "Review data requirements for Chapter {N}"
    |               Timeout: 5 min
    |
    +---> [2c] Fan-in: Wait for all consultants
    |         |
    |         +---> Handle failures (see Error Handling below)
    |
    +---> [2d] Aggregate recommendations
    |         |
    |         +---> If conflict detected:
    |               AskUserQuestion: "Consultants disagree on {topic}.
    |               Statistician: {opinion1}
    |               Mathematician: {opinion2}
    |               Which approach do you prefer?"
    |
    +---> [2e] Generate chapter{N}-notebook-plans.md
    |
    v
[3] Update session state
    |
    +---> Set current_phase: 2
    +---> Set completed_phases: [0, 1, 2]
    +---> Record outputs.chapter_plans
    |
    v
END Phase 2 -> Quality Gate 2
```

### Task Tool Invocation Pattern

```python
# Fan-out pattern for Phase 2 consultants
tasks = []

# Spawn all consultants in parallel
statistician_task = Task(
    prompt=f"Review statistical approaches for: {chapter_analyses}",
    agent="statistician-consultant",
    timeout=300000  # 5 minutes
)
tasks.append(statistician_task)

mathematician_task = Task(
    prompt=f"Review algorithms for: {chapter_analyses}",
    agent="mathematician-consultant",
    timeout=300000
)
tasks.append(mathematician_task)

programmer_task = Task(
    prompt=f"Review data requirements for: {chapter_analyses}",
    agent="programmer-consultant",
    timeout=300000
)
tasks.append(programmer_task)

# Wait for all tasks (fan-in handled by Task tool)
results = await asyncio.gather(*[t.run() for t in tasks], return_exceptions=True)
```

### Error Handling for Fan-Out

```
For each consultant result:
    |
    +---> If SUCCESS: Add to aggregation
    |
    +---> If FAILURE:
          |
          +---> Retry once (wait 30s)
          |
          +---> If retry fails:
                |
                +---> If statistician (critical):
                |     Escalate to user with Template 1
                |
                +---> If mathematician/programmer (optional):
                      Log warning, proceed without
```

### Quality Gate 2 Criteria

- [ ] All chapters have notebook plans
- [ ] No unresolved critical conflicts
- [ ] Each notebook has statistical approach defined
- [ ] Data flow is consistent (outputs match downstream inputs)

---

## Phase 3: Structure Review

**Duration**: ~5 minutes
**Owner**: structure-reviewer (Haiku)
**Timeout**: 10 minutes

### Workflow Steps

```
START Phase 3
    |
    v
[1] Read all planning artifacts
    |
    +---> research-structure.md
    +---> All chapter{N}-notebook-plans.md
    |
    v
[2] Check for issues:
    |
    +---> Missing dependencies
    +---> Redundant analyses
    +---> Logical flow problems
    +---> Incomplete specifications
    +---> Missing quality controls
    |
    v
[3] Categorize issues by severity
    |
    +---> Critical: Blocks execution
    +---> Major: Affects quality
    +---> Minor: Improvement opportunity
    |
    v
[4] Generate structure-review-report.md
    |
    v
[5] Present to user (USER APPROVAL GATE 1)
    |
    AskUserQuestion:
    "Structure Review Complete

    Summary:
    - {N} chapters planned
    - {M} notebooks total
    - {K} issues identified (X critical, Y major, Z minor)

    Critical Issues:
    {list}

    Approve / Request changes / Reject? [A/c/r]"
    |
    +---> If "A": Proceed to Phase 4
    |
    +---> If "c":
    |     Ask "What changes?"
    |     Route to appropriate phase
    |     Re-run affected phases
    |     Return to this gate
    |
    +---> If "r":
          Ask "Return to Phase 1 or 2? [1/2]"
          Jump to selected phase
    |
    v
[6] Update session state
    |
    +---> Record user_approvals.phase_3
    |
    v
END Phase 3 -> Phase 4
```

### Quality Gate 3 Criteria (User Approval)

- [ ] User reviewed structure summary
- [ ] User approved or requested specific changes
- [ ] All critical issues addressed (if any)

---

## Phase 4: Notebook Review

**Duration**: ~10 minutes
**Owner**: notebook-reviewer (Sonnet 4.5)
**Timeout**: 15 minutes total

### Workflow Steps

```
START Phase 4
    |
    v
[1] Read all chapter plans
    |
    v
[2] Fan-out: One reviewer per chapter (PARALLEL)
    |
    +---> For each chapter:
          |
          Task: notebook-reviewer
          Prompt: "Review notebook plans for Chapter {N}:
                   {chapter_plan_content}

                   Check:
                   - Pseudocode completeness
                   - Statistical correctness
                   - Data flow consistency
                   - Edge case coverage"
          Timeout: 5 min per chapter
    |
    v
[3] Fan-in: Aggregate review results
    |
    +---> Handle partial failures (proceed with available)
    |
    v
[4] Generate notebook-review-report.md
    |
    v
[5] Present to user (USER APPROVAL GATE 2)
    |
    AskUserQuestion:
    "Notebook Review Complete

    Per-Chapter Summary:
    - Chapter 1: {N} notebooks, {K} issues
    - Chapter 2: {N} notebooks, {K} issues
    ...

    Critical Issues:
    {list}

    Approve / Request changes / Reject? [A/c/r]"
    |
    +---> Handle responses (same as Phase 3)
    |
    v
[6] Update session state
    |
    +---> Record user_approvals.phase_4
    |
    v
END Phase 4 -> Phase 5
```

---

## Phase 5: Notebook Generation

**Duration**: ~7 minutes
**Owner**: notebook-generator (Sonnet 4.5)
**Timeout**: 15 minutes total

### Workflow Steps

```
START Phase 5
    |
    v
[1] Read all approved chapter plans
    |
    v
[2] Fan-out: One generator per chapter (PARALLEL)
    |
    +---> For each chapter:
          |
          Task: notebook-generator
          Prompt: "Generate .ipynb notebooks for Chapter {N}:
                   {chapter_plan_content}

                   Output to: {output_dir}/chapter{N}_{slug}/"
          Timeout: 7 min per chapter
    |
    v
[3] Fan-in: Collect generation results
    |
    +---> Track success/failure per chapter
    |
    v
[4] Validate generated notebooks
    |
    +---> For each .ipynb:
          python3 -c "
          import nbformat
          nb = nbformat.read('{path}', as_version=4)
          nbformat.validate(nb)
          "
    |
    v
[5] Copy to session directory (backup)
    |
    v
[6] Handle partial failures
    |
    +---> If some chapters failed:
          AskUserQuestion:
          "{M} of {N} chapters generated successfully.
          Failed: {list}

          (A) Proceed with available notebooks
          (B) Retry failed chapters
          (C) Abort"
    |
    v
[7] Update session state
    |
    +---> Record outputs.notebooks
    |
    v
END Phase 5 -> Quality Gate 5
```

### Notebook Generation Details

```python
# Generate valid .ipynb structure
import nbformat

notebook = nbformat.v4.new_notebook()

# Add metadata
notebook.metadata['kernelspec'] = {
    'display_name': 'Python 3',
    'language': 'python',
    'name': 'python3'
}

# Add title cell
notebook.cells.append(nbformat.v4.new_markdown_cell(
    f"# {notebook_title}\n\n## Goal\n{analysis_goal}"
))

# Add code cells with pseudocode
for cell in pseudocode_cells:
    code = f"# {cell['type'].upper()}\n"
    code += f"# {cell['pseudocode']}\n"
    code += "# TODO: Implement\n"
    notebook.cells.append(nbformat.v4.new_code_cell(code))

# Validate before saving
nbformat.validate(notebook)

# Write to file
with open(notebook_path, 'w') as f:
    nbformat.write(notebook, f)
```

### Quality Gate 5 Criteria

- [ ] All notebooks created (or partial with user approval)
- [ ] All notebooks pass nbformat validation
- [ ] Backup copies exist in session directory
- [ ] File naming follows convention

---

## Phase 6: Statistical Fact-Checking

**Duration**: ~10-20 minutes
**Owner**: statistical-fact-checker (Sonnet 4.5)
**Timeout**: 30 minutes

### Workflow Steps

```
START Phase 6
    |
    v
[1] Read all generated notebooks
    |
    v
[2] Analyze for statistical concerns
    |
    +---> Test mismatches
    +---> Multiple testing issues
    +---> Interpretation errors
    +---> Assumption violations
    +---> Effect size gaps
    |
    v
[3] Categorize concerns by severity
    |
    +---> Critical: Incorrect conclusions
    +---> Standard: Best practice violation
    +---> Minor: Improvement opportunity
    |
    v
[4] Count total concerns
    |
    +---> If 0 concerns:
    |     Show justification
    |     Proceed to completion
    |
    +---> If <= 5 concerns:
    |     Enter interview mode (one at a time)
    |
    +---> If > 5 concerns:
          Show summary first
          Offer batch options
          Then interview mode
    |
    v
[5] INTERVIEW MODE (see interview-protocol.md)
    |
    v
[6] Accumulate decisions
    |
    +---> accepted: list
    +---> rejected: list
    +---> skipped: list
    |
    v
[7] Generate reports
    |
    +---> statistical-review-report.md
    +---> corrections-manifest.json
    |
    v
[8] Apply corrections (if user approves)
    |
    AskUserQuestion:
    "Summary:
    - {X} accepted, {Y} rejected, {Z} skipped

    Apply all accepted corrections now? [yes/no]"
    |
    +---> If yes: Regenerate affected notebooks
    +---> If no: Save manifest for later
    |
    v
[9] Update session state
    |
    +---> Set status: "completed"
    |
    v
END Phase 6 -> Workflow Complete
```

### Quality Gate 6 Criteria (User Approval)

- [ ] All concerns reviewed (or explicitly skipped)
- [ ] User approved correction application decision
- [ ] Statistical review report generated
- [ ] Corrections manifest saved
