# Post-Work Invocation Examples

## Pattern B: Orchestrator Post-Work

A calling orchestrator invokes git-strategy-advisor via Task tool after completing work.

### Task Tool Prompt

```
Use git-strategy-advisor to determine git strategy for completed work.

mode: post-work

Git status:
 M src/auth/validator.py
 M src/auth/tokens.py
?? src/auth/refresh.py

Git diff stat:
 src/auth/validator.py | 45 +++++++++++++++++++++++++++++++++------------
 src/auth/tokens.py    | 12 ++++++------
 src/auth/refresh.py   | 38 ++++++++++++++++++++++++++++++++++++++

Current branch: main
Remote: origin (configured)
```

### Expected Output

```yaml
git_strategy_recommendation:
  version: "1.0"
  timestamp: "2026-02-07T15:30:00Z"
  mode: "post-work"
  confidence: "high"

  analysis:
    files_changed: 3
    lines_changed: 95
    directories_spanned: 1
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
      suggested_name: "feature/auth-token-refresh"
      rationale: "3 files, 95 lines in 1 directory (moderate scope by line count) on main -- feature branch isolates risk and enables PR review."
    push:
      action: "push-now"
      rationale: "Moderate scope with remote configured -- push to back up work and enable collaboration."
    pr:
      action: "create-pr"
      rationale: "Moderate scope warrants peer review via pull request."

  warnings: []

  summary: "Create feature branch `feature/auth-token-refresh`, push to remote, and open a PR for review."
```

---

## Clean Working Tree Example

Changes already committed to a feature branch but not pushed.

### Scenario

```
$ git status
On branch feature/add-caching
nothing to commit, working tree clean

$ git log origin/main..HEAD --oneline
a1b2c3d feat(cache): add Redis connection pool
d4e5f6g feat(cache): implement LRU eviction policy

$ git diff --stat origin/main..HEAD
 src/cache/pool.py     | 120 ++++++++++++++++++++++++++
 src/cache/eviction.py |  85 +++++++++++++++++++
 src/cache/config.yaml |  15 ++++
 tests/test_cache.py   |  95 +++++++++++++++++++++
 4 files changed, 315 insertions(+)
```

### Expected Output

```yaml
git_strategy_recommendation:
  version: "1.0"
  timestamp: "2026-02-07T15:45:00Z"
  mode: "post-work"
  confidence: "high"

  analysis:
    files_changed: 4
    lines_changed: 315
    directories_spanned: 2
    work_type: "code"
    scope: "major"

  context:
    current_branch: "feature/add-caching"
    primary_branch: "main"
    is_on_main: false
    has_remote: true

  strategy:
    branch:
      action: "use-current-branch"
      rationale: "4 files, 315 lines across 2 directories (major scope). Changes are already committed on feature/add-caching."
    push:
      action: "push-now"
      rationale: "Major scope with remote configured -- push to back up work immediately."
    pr:
      action: "create-pr"
      rationale: "Major scope warrants peer review. Consider breaking into smaller PRs for better review quality (research shows review effectiveness drops significantly above 200 lines)."

  warnings:
    - code: "UNPUSHED_COMMITS"
      message: "Analysis based on 2 unpushed commits on feature/add-caching."

  summary: "Push feature/add-caching to remote and open a PR for review. Consider splitting into smaller PRs given 315 lines changed."
```

---

## Bookend Pattern D (Second Half): Post-Work Re-evaluation

Continuing from the pre-work example in `examples/pre-work-invocation.md`, the actual work exceeded the estimate.

### Pre-Work Estimate (from earlier)

Estimated: 1 file, ~20 lines (minor documentation change).

### Actual Changes

```
$ git diff --numstat HEAD
25   5   README.md
48   0   docs/installation.md
32   0   docs/configuration.md
```

### Task Tool Prompt

```
Use git-strategy-advisor to re-evaluate git strategy after completing work.

mode: post-work

Previously estimated as minor documentation update (1 file, ~20 lines).
Actual scope grew during implementation.
```

### Expected Output

```yaml
git_strategy_recommendation:
  version: "1.0"
  timestamp: "2026-02-07T16:00:00Z"
  mode: "post-work"
  confidence: "high"

  analysis:
    files_changed: 3
    lines_changed: 110
    directories_spanned: 2
    work_type: "documentation"
    scope: "moderate"

  context:
    current_branch: "docs/update-installation-instructions"
    primary_branch: "main"
    is_on_main: false
    has_remote: true

  strategy:
    branch:
      action: "use-current-branch"
      rationale: "3 files, 110 lines across 2 directories (moderate scope, upgraded from pre-work minor estimate). Already on docs/update-installation-instructions."
    push:
      action: "push-now"
      rationale: "Moderate scope with remote configured -- push to back up work and enable collaboration."
    pr:
      action: "create-pr"
      rationale: "Scope upgraded from minor to moderate during implementation. Moderate scope warrants peer review via pull request."

  warnings:
    - code: "SCOPE_UPGRADE"
      message: "Post-work scope (moderate) exceeds pre-work estimate (minor). PR now recommended."

  summary: "Push docs/update-installation-instructions to remote and open a PR. Scope grew from minor to moderate during implementation."
```
