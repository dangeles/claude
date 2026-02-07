# Pre-Work Invocation Examples

## Pattern A: Orchestrator Pre-Work

A calling orchestrator invokes git-strategy-advisor via Task tool before starting work.

### Task Tool Prompt

```
Use git-strategy-advisor to determine git strategy for the following planned work.

mode: pre-work

Task: Add input validation module to the data pipeline. Will create 3-4 new
Python files in src/validation/ and update src/pipeline/main.py. Estimated
~150 lines of new code plus ~30 lines of modifications.

Current state: On main branch, remote configured, clean working tree.
```

### Expected Output

```yaml
git_strategy_recommendation:
  version: "1.0"
  timestamp: "2026-02-07T14:00:00Z"
  mode: "pre-work"
  confidence: "medium"

  analysis:
    files_changed: 5
    lines_changed: 180
    directories_spanned: 2
    work_type: "code"
    scope: "moderate"

  context:
    current_branch: "main"
    primary_branch: "main"
    is_on_main: true
    has_remote: true

  strategy:
    branch:
      action: "create-feature-branch"
      suggested_name: "feature/add-input-validation-pipeline"
      rationale: "~5 files, ~180 lines across 2 directories (moderate scope) on main -- feature branch isolates risk and enables PR review."
    push:
      action: "push-now"
      rationale: "Moderate scope with remote configured -- push to back up work and enable collaboration."
    pr:
      action: "create-pr"
      rationale: "Moderate scope warrants peer review via pull request."

  warnings: []

  summary: "Create feature branch `feature/add-input-validation-pipeline`, push to remote, and open a PR for review."
```

---

## Pattern C: Standalone Pre-Work

User asks for git advice before starting work.

### User Prompt

```
What git strategy should I use for adding a new validation module to the
data pipeline? I expect to create 3-4 new files and modify one existing file.
```

### Expected Human-Readable Output

```
Git Strategy Recommendation
---
Based on your planned changes (~5 files, ~180 lines across 2 directories):

  Branch: Create feature branch `feature/add-validation-module`
  Push:   Push to remote when ready
  PR:     Create a pull request

  Confidence: Medium (pre-work estimate)
  Rationale: Moderate code changes on main branch warrant isolation
             and peer review.

Full recommendation written to: /tmp/git-strategy-recommendation.yaml
```

---

## Bookend Pattern D (First Half): Pre-Work Setup

An orchestrator invokes git-strategy-advisor at the start of a task. The recommendation will be re-evaluated post-work (see `examples/post-work-invocation.md`).

### Task Tool Prompt

```
Use git-strategy-advisor to determine git strategy before starting work.

mode: pre-work

Task: Update README.md with new installation instructions.
Estimated: 1 file, ~20 lines changed.
```

### Expected Output

```yaml
git_strategy_recommendation:
  version: "1.0"
  timestamp: "2026-02-07T14:05:00Z"
  mode: "pre-work"
  confidence: "medium"

  analysis:
    files_changed: 1
    lines_changed: 20
    directories_spanned: 1
    work_type: "documentation"
    scope: "minor"

  context:
    current_branch: "main"
    primary_branch: "main"
    is_on_main: true
    has_remote: true

  strategy:
    branch:
      action: "create-feature-branch"
      suggested_name: "docs/update-installation-instructions"
      rationale: "1 file, 20 lines (minor scope) on main -- feature branch for documentation update on main."
    push:
      action: "push-after-review"
      rationale: "Minor scope with remote -- push after self-review."
    pr:
      action: "no-pr-needed"
      rationale: "Minor documentation change does not require PR review."

  warnings: []

  summary: "Create branch `docs/update-installation-instructions`, push after self-review. No PR needed for minor docs update."
```
