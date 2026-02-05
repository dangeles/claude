---
name: scientific-analysis-architect
description: Use when planning multi-chapter scientific research analyses with expert consultation. Produces Jupyter notebooks with pseudocode for RNA-seq, proteomics, or other data analysis workflows. Triggers on research planning, analysis architecture, scientific notebook generation, or multi-chapter analysis design requests.
version: 1.0.0
tags: [scientific-analysis, multi-agent, jupyter, research-planning, pseudocode]
---

# scientific-analysis-architect

Multi-phase workflow for planning scientific research analyses using Jupyter notebooks with pseudocode. Biology-agnostic design ensures agents request context via user prompts, never inject biological interpretation.

## When to Use

- Planning multi-chapter scientific data analysis (RNA-seq, proteomics, imaging)
- Need expert consultation (statistician, mathematician, programmer perspectives)
- Want Jupyter notebooks with pseudocode for implementation
- Research project requires 3-7 chapters of analysis

## When NOT to Use

- Need actual code implementation (use programming-pm after this skill)
- Need literature review (use lit-pm skill)
- Single analysis without chapter structure
- Already have detailed analysis plan

## Workflow Overview

```
User Request
     |
+----v--------------------+
| Phase 0: Initialization | ~2 min
| Session setup, validation|
| - Jupytext availability |
| - Output dir validation  |
+----+--------------------+
     |
+----v--------------------+
| Phase 1: Birds-Eye      | ~12 min
| Planning                |
| [research-architect]    |
+----+--------------------+
     |
research-structure.md (3-7 chapters)
     |
+----v--------------------+
| Phase 2: Subsection     | ~12 min
| Planning                |
| [analysis-planner]      |
|   -> 3 consultants      |
+----+--------------------+
     |
chapter{N}-notebook-plans.md
     |
+----v--------------------+
| Phase 3: Structure      | ~5 min
| Review                  |
| [structure-reviewer]    |
+----+--------------------+
     |
[USER APPROVAL GATE 1]
     |
+----v--------------------+
| Phase 4: Notebook       | ~10 min
| Review (parallel)       |
| [notebook-reviewer]     |
+----+--------------------+
     |
[USER APPROVAL GATE 2]
     |
+----v--------------------+
| Phase 5: Notebook       | ~7 min
| Generation (parallel)   |
| [notebook-generator]    |
| - nbformat validation   |
+----+--------------------+
     |
.ipynb files with pseudocode
     |
+----v--------------------+
| Phase 6: Statistical    | ~10-20 min
| Fact-Checking           |
| [statistical-fact-checker]
| INTERVIEW MODE          |
+----+--------------------+
     |
Corrected notebooks (final)
```

**Estimated Runtime**: 56-76 minutes for 3 chapters

## Phase 0: Initialization

**Owner**: Orchestrator
**Duration**: 2-5 minutes
**Checkpoint**: Never (automatic)

1. **Check dependencies**:
   ```bash
   # Verify jupytext is available
   python3 -c "import jupytext" 2>/dev/null || {
     echo "WARNING: jupytext not installed. Install with: pip install jupytext"
     echo "Notebooks will be created without Jupytext metadata."
   }

   # Verify nbformat is available (required)
   python3 -c "import nbformat" || {
     echo "ERROR: nbformat not installed. Required for notebook validation."
     echo "Install with: pip install nbformat"
     exit 1
   }
   ```

2. **Create session directory**:
   - Primary: `{output_directory}/.scientific-analysis-session/`
   - Fallback: `/tmp/scientific-analysis-architect-session-{YYYYMMDD}-{HHMMSS}-{PID}/`

3. **Validate output directory**:
   - Check exists and writable
   - Perform write test
   - If fails, offer alternatives

4. **Initialize session state**:
   - Create `session-state.json`
   - Set status: "initialized"

**Quality Gate 0**: Session directory created, output directory validated, nbformat available.

## Phase 1: Birds-Eye Planning

**Owner**: research-architect (Sonnet 4.5)
**Duration**: ~12 minutes
**Timeout**: 15 minutes

1. Ask user: "Please describe your dataset and research goals"
2. If uncertainty detected ("not sure", "maybe"):
   - Fan-out to analysis-brainstormer and method-brainstormer (Haiku)
   - Present brainstorming suggestions
3. Generate research-structure.md with 3-7 chapters

**Biology-Agnostic Behavior**:
- Agents ASK: "What biological questions are you trying to answer?"
- Agents DO NOT inject: "You should look at cell types"

**Output**: `{session_dir}/research-structure.md`

**Quality Gate 1**: Structure has 3-7 chapters, each with goal and analyses.

## Phase 2: Subsection Planning

**Owner**: analysis-planner (Sonnet 4.5)
**Duration**: ~12 minutes for 3 chapters
**Timeout**: 20 minutes total

For each chapter:
1. Fan-out to expert panel (parallel, all Haiku):
   - statistician-consultant: Statistical approach validation
   - mathematician-consultant: Algorithm requirements
   - programmer-consultant: Data requirements

2. Fan-in: Aggregate recommendations
3. If consultants disagree, present conflict to user
4. Generate chapter{N}-notebook-plans.md

**Consolidated Escalation** (if parallel failures):
- Wait for all retries before escalating
- Single prompt with all failures
- Statistician is critical; others are optional

**Output**: `{session_dir}/chapter{N}-notebook-plans.md` (one per chapter)

**Quality Gate 2**: All chapters have notebook plans, no unresolved conflicts.

## Phase 3: Structure Review

**Owner**: structure-reviewer (Haiku)
**Duration**: ~5 minutes
**Timeout**: 10 minutes

1. Review research-structure.md and all chapter plans
2. Check for missing dependencies, redundancies, logical issues
3. Generate structure-review-report.md

**Output**: `{session_dir}/structure-review-report.md`

**USER APPROVAL GATE 1**:
```
Structure Review Complete

Summary:
- {N} chapters planned
- {M} notebooks total
- {K} issues identified

Approve / Request changes / Reject? [A/c/r]
```

## Phase 4: Notebook Review

**Owner**: notebook-reviewer (Sonnet 4.5)
**Duration**: ~10 minutes
**Timeout**: 15 minutes total

1. Fan-out: One reviewer per chapter (parallel)
2. Check pseudocode completeness, statistical correctness, data flow
3. Fan-in: Aggregate review reports

**Output**: `{session_dir}/notebook-review-report.md`

**USER APPROVAL GATE 2**:
```
Notebook Review Complete

Per-Chapter Summary:
- Chapter 1: {N} notebooks, {K} issues
...

Approve / Request changes / Reject? [A/c/r]
```

## Phase 5: Notebook Generation

**Owner**: notebook-generator (Sonnet 4.5)
**Duration**: ~7 minutes
**Timeout**: 15 minutes total

1. Fan-out: One generator per chapter (parallel)
2. Create .ipynb files with pseudocode cells
3. Write to both session directory and output directory
4. Fan-in: Verify all notebooks created

**Partial Completion Handling**:
- If some chapters fail, offer to proceed with available
- Enable per-chapter regeneration later

**Output**:
- `{output_dir}/chapter{N}_{slug}/notebook{N}_{M}_{slug}.ipynb`
- `{session_dir}/notebooks/` (backup)

**Quality Gate 5**: All notebooks valid .ipynb format.

**Validation**:
```python
import nbformat
for notebook_path in generated_notebooks:
    with open(notebook_path) as f:
        nb = nbformat.read(f, as_version=4)
        nbformat.validate(nb)
        print(f"VALID: {notebook_path}")
```

## Phase 6: Statistical Fact-Checking

**Owner**: statistical-fact-checker (Sonnet 4.5)
**Duration**: ~10-20 minutes
**Timeout**: 30 minutes

**INTERVIEW MODE**:

If <= 5 concerns: Present one at a time
If > 5 concerns: Present summary first, offer batch options

**Concern Format**:
```
Statistical Concern {N} of {total}

Notebook: {path}
Cell: {number}

Issue: {description}

Current: {current_code}
Recommendation: {recommended_fix}

Accept? [yes/no/skip/explain]
```

**Batch Options** (after 5 concerns):
- Continue one-by-one
- Accept all remaining
- Reject all remaining
- Accept critical/standard, skip minor

**After Interview**:
```
Summary:
- {X} accepted, {Y} rejected, {Z} skipped

Apply corrections? [yes/no]
```

**Output**:
- `{session_dir}/statistical-review-report.md`
- `{session_dir}/corrections-manifest.json`
- Updated .ipynb files (if corrections applied)

## Session Management

### Session Directory Structure

```
{session_dir}/
├── session-state.json          # Resumable state
├── research-structure.md       # Phase 1 output
├── chapter1-notebook-plans.md  # Phase 2 output
├── chapter2-notebook-plans.md
├── structure-review-report.md  # Phase 3 output
├── notebook-review-report.md   # Phase 4 output
├── statistical-review-report.md # Phase 6 output
├── corrections-manifest.json   # Phase 6 corrections
├── notebooks/                  # Backup copies
│   ├── chapter1_data-atlas/
│   │   └── *.ipynb
│   └── chapter2_hypothesis-testing/
│       └── *.ipynb
└── logs/
    └── workflow.log
```

### Resume Protocol

On skill invocation:
1. Check for existing sessions (output_dir first, then /tmp)
2. If found and < 72 hours old:
   ```
   Found incomplete session from {timestamp}
   Project: {research_goals}
   Status: Phase {N}

   Resume? [yes/no]
   ```
3. If yes: Load state, continue from current phase
4. If no: Archive old session, start new

### Interrupt Handling

On Ctrl+C:
1. Save current state with status: "interrupted"
2. Print: "Session saved. Resume with: /scientific-analysis-architect"

## Error Handling

See [error-handling.md](references/error-handling.md) for complete specification.

### Timeout Configuration

| Phase | Timeout | Exceeded Action |
|-------|---------|-----------------|
| 0 | 5 min | Abort |
| 1 | 15 min | Escalate to user |
| 2 | 20 min | Proceed with available consultants |
| 3 | 10 min | Escalate to user |
| 4 | 15 min | Proceed with available reviews |
| 5 | 20 min | Proceed with partial, offer retry |
| 6 | 30 min | Pass with uncertainty note |

### Retry Protocol

- First failure: Wait 30s, retry automatically
- Second failure: Ask user (proceed without or abort)
- Maximum 2 retries per agent

### Circuit Breaker

- Open after 2 consecutive failures per agent
- Action: Escalate to user
- Reset: On successful execution

## Quality Gates Summary

| Gate | Phase | Owner | Pass Criteria |
|------|-------|-------|---------------|
| 0 | 0 | Orchestrator | Session created, output validated, nbformat available |
| 1 | 1 | research-architect | 3-7 chapters with goals |
| 2 | 2 | analysis-planner | All chapter plans, no critical conflicts |
| 3 | 3 | User | Approve structure |
| 4 | 4 | User | Approve notebook plans |
| 5 | 5 | notebook-generator | Valid .ipynb files (nbformat validated) |
| 6 | 6 | User | Interview complete, corrections applied |

## Dependencies

- **Tools**: Task, AskUserQuestion, Read, Write, Bash
- **Python Packages**: nbformat (required), jupytext (optional)
- **Complements**: lit-pm (literature), programming-pm (implementation)
- **Output**: Pseudocode notebooks for manual or programming-pm implementation

## References

- [agent-definitions.md](references/agent-definitions.md)
- [phase-workflows.md](references/phase-workflows.md)
- [interview-protocol.md](references/interview-protocol.md)
- [notebook-templates.md](references/notebook-templates.md)
- [session-schema.md](references/session-schema.md)
- [error-handling.md](references/error-handling.md)
- [quality-gates.md](references/quality-gates.md)

## Examples

- [rnaseq-analysis-plan.md](examples/rnaseq-analysis-plan.md)
- [statistical-interview-session.md](examples/statistical-interview-session.md)
