---
name: archive-workflow
last_updated: 2026-05-24
description: Use when organizing projects — detecting clutter, enforcing naming conventions, structuring directories, managing gitignore, or assessing project expandability.
---

# archive-workflow

A project organization skill with 1 PM orchestrator (library-pm) and 5 specialist agents, using a 4-wave workflow to manage clutter, naming, structure, and expandability across any project type.

## Overview

library-pm dispatches up to 4 READ-ONLY analyst agents across three waves, then hands off to a single WRITE executor agent for synthesis and execution.

**Architecture**: Hub-and-spoke multi-agent coordination
**Pattern**: CQRS (Command Query Responsibility Segregation) - analysts READ, integrator WRITES

## When to Use

- Organizing a new project with proper structure
- Cleaning up an existing project with accumulated clutter
- Enforcing consistent naming conventions
- Reviewing gitignore for missing patterns
- Assessing project scalability and modularity
- Major project reorganization before release

## When NOT to Use

- Code refactoring (organizing structure is in scope; rewriting code is not)
- Content creation (where documentation goes is in scope; writing it is not)
- CI/CD pipeline configuration
- Database schema design
- Secret management beyond gitignore context

## Delegation

You are library-pm, the orchestrator. Delegate specialist work when the subtask is
substantial and independent — a full clutter scan of an unfamiliar tree, a naming audit
against project-type conventions, a structure proposal plus migration plan, an
expandability review, or execution of an approved plan. Handle it directly when you could
finish it in a handful of tool calls, when the work is sequential, or when you need the
context in your own loop. See `../references/delegation-and-scope.md`.

Scale the pipeline to the request. A narrow ask — review the gitignore, rename one file,
check whether a directory name fits the convention — does not need all four analysts;
run the relevant check and skip the waves that would return nothing. Reserve the full
4-wave fan-out for whole-project reorganization, where the analyst tracks are genuinely
independent and each has real ground to cover.

You own: dispatching agents per wave, evaluating quality gates, presenting execution
plans for user approval, session state, and conflict resolution between analyst outputs
(applying the Conflict Resolution Rules).

Repo-modifying file operations go through archive-decision-integrator so they land in the
execution log.

## State Anchoring

Start each response with `[Wave N/4 - {wave_name}] {brief status}`. The authoritative
record of progress is `workflow-state.yaml` in the session directory; trust it over
memory.

## Pre-flight Checks (Phase 0)

Before any analysis or file operations:

### 1. Git Status Check
```bash
git status --porcelain
```
- If empty: PROCEED
- If non-empty: STOP and prompt user (stash, commit, or abort)

### 2. Detached HEAD Check
- If detached: STOP ("Cannot run archive-workflow in detached HEAD state")

### 3. Permissions Check
- Verify write permissions on project directory

### 4. Disk Space Check
- Warn if <100MB available in /tmp

## The 4-Wave Pipeline

### Wave 1: Clutter Analysis (Sequential)
**Agent**: archive-clutter-analyst (READ-ONLY)

Detects:
- Generated files (node_modules/, __pycache__, build/, dist/, .venv)
- Stale content (old branches, abandoned experiments, deprecated code)
- Organizational mess (duplicates, misplaced files, temp files)

**Output**: clutter-report.md

### Wave 2: Parallel Organization (Parallel)
**Agents**: archive-nomenclature-enforcer + archive-structure-organizer (both READ-ONLY)

**nomenclature-enforcer**:
- Audits file/directory naming against project-type conventions
- Detects existing patterns (adaptive mode)
- Output: naming-violations.md

**structure-organizer**:
- Analyzes current structure vs project-type template
- Prescriptive for new projects, adaptive for existing
- Output: structure-proposal.md

### Wave 3: Expandability Review (Sequential)
**Agent**: archive-expandability-reviewer (READ-ONLY)
**Input**: structure-proposal.md from Wave 2

Assesses:
- Scalability (can structure handle 10x files? 5x contributors?)
- Modularity (components decoupled? extension points?)
- Coupling issues

**Output**: expandability-assessment.md

### Wave 4: Synthesis & Execution (Sequential, User Approval Required)
**Agent**: archive-decision-integrator (READ + WRITE)

**Inputs**: the analyst reports produced in Waves 1-3
**Process**:
1. Merge outputs, apply conflict resolution rules
2. Generate execution-plan.md for user review
3. Get user approval (APPROVE ALL / APPROVE WITH EXCLUSIONS / REJECT)
4. Execute file operations
5. Generate documentation via editor skill
6. Produce execution-log.md and final-organization-report.md
7. Generate `.archive-metadata.yaml` in the project root, per the schema in
   `references/archive-metadata-schema.md`:
   a. Write to `.archive-metadata.yaml.tmp` first (atomic write pattern)
   b. Populate from Wave 1-3 analysis results:
      - project.type from Wave 1 detection heuristics
      - naming_conventions from Wave 2 nomenclature-enforcer results
      - structure from Wave 2 structure-organizer results
   c. Double-quote all string values
   d. Make all paths in full_reference absolute (expand ~ at generation time)
   e. Validate the temp file: parse with yaml.safe_load, verify required fields
   f. If validation passes: mv .archive-metadata.yaml.tmp .archive-metadata.yaml
   g. If validation fails: log ERROR, remove temp file, continue without metadata
8. Check .gitignore: if `git check-ignore -q .archive-metadata.yaml`, WARN "metadata file is in .gitignore"
9. Persist final-organization-report.md to docs/organization/ in the repo
   (create directory if needed)
10. Stage both files: git add .archive-metadata.yaml docs/organization/final-organization-report.md
    These are part of the existing Wave 4 commit (NOT a separate commit).

**Circular Dependency Prevention**: When invoking any specialist (editor, etc.) from
within archive-workflow, include in the Task tool handoff: `archival_context: "skip"`.
This prevents specialists from checking a stale/non-existent .archive-metadata.yaml
during archive-workflow's own execution.

## Quality Gates and Degradation

Each gate is an objective file check on the session directory, evaluated when the
preceding wave returns.

| Gate | Phase | Checks | On Failure |
|------|-------|--------|------------|
| QG1 | Phase 1 | Project type detected, session initialized | Escalate |
| QG2 | Wave 1 | clutter-report.md exists, has Summary section | Retry once, then proceed without |
| QG3 | Wave 2 | Both Wave 2 reports exist | Proceed with available; escalate if structure-organizer failed |
| QG4 | Wave 3 | expandability-assessment.md exists | Proceed with advisory flag |
| QG5 | Wave 4 | execution-log.md exists, no ERROR entries | Rollback + escalate |

If an agent stalls or returns nothing usable, retry it once, then continue with the
reports you have and note the gap in the execution plan.

## Execution Plan Approval

Before Wave 4 executes file operations, present to user:

**Category A** (Non-destructive): Renames, moves
- Execute unless user explicitly excludes

**Category B** (Clutter cleanup): Gitignore additions
- Execute unless user explicitly excludes

**Category C** (Deletions): File removals
- REQUIRES EXPLICIT APPROVAL per file

User options:
- APPROVE ALL
- APPROVE WITH EXCLUSIONS (specify files to skip)
- REJECT (abort workflow)

## Conflict Resolution Rules

Apply in order:

1. **Rename + Move**: Combined operation
   - `git mv old_path new_dir/new_name`

2. **Naming vs Structure Directory Name**: Nomenclature wins on naming questions
   - If structure proposes `/Data/`, but naming says kebab-case, use `/data/`

3. **Placement vs Naming**: Structure wins on file placement
   - File goes to best-fit location even if name doesn't perfectly match

4. **Expandability Concerns**: Adjust if critical issue flagged
   - Modify proposal before executing, log reasoning

5. **Clutter Priority**: Process clutter first
   - Add to .gitignore before organizing remaining files

6. **Unresolvable Naming Conflict**: escalate to user
   - When nomenclature and structure propose different names for same file
   - Present both options, document decision

7. **Multi-Analyst Existence Conflict**: escalate to user
   - When clutter says delete, but other analyst needs the file
   - NEVER auto-delete; require explicit user confirmation

## Rollback Procedure

### Case 1: Failure During Wave 4 Execution (Pre-Commit)

If failure occurs BEFORE final commit (during file operations):

**Step 1**: Stop execution immediately

**Step 2**: Check git status:
```bash
git status --porcelain
```

**Step 3**: Revert all unstaged changes:
```bash
git restore .
```

**Step 4**: Unstage all staged changes:
```bash
git restore --staged .
```

**Step 5**: Remove any newly created files (not in git):
```bash
# List new files
git status --porcelain | grep "^??" | cut -c4-

# Manually review and delete if appropriate
```

**Step 6**: Verify clean state:
```bash
git status
# Should show: "nothing to commit, working tree clean"
```

**Step 7**: Clean up partial metadata:
```bash
rm -f .archive-metadata.yaml.tmp
# If .archive-metadata.yaml was created in this session:
git rm -f .archive-metadata.yaml 2>/dev/null || rm -f .archive-metadata.yaml
```

**Step 8**: Clean session directory:
```bash
rm -rf /tmp/archive-workflow-session-{id}/
```

### Case 2: Failure After Commit (Post-Commit)

If failure occurs AFTER final commit (during testing or validation):

**Step 1**: Identify the commit before archive-workflow:
```bash
git log --oneline | head -5
# Find the commit SHA before the archive-workflow commit
```

**Step 2**: Hard reset to that commit:
```bash
git reset --hard <SHA-before-archive-workflow>
```

**Step 3**: Force push if already pushed to remote (ONLY if you're sure):
```bash
# WARNING: Destructive operation - confirm with user first
git push --force origin main
```

**Step 4**: Clean up partial metadata:
```bash
rm -f .archive-metadata.yaml.tmp
```

**Step 5**: Clean session directory:
```bash
rm -rf /tmp/archive-workflow-session-{id}/
```

## Session Directory Structure

```
/tmp/archive-workflow-session-{YYYYMMDD-HHMMSS-PID}/
├── workflow-state.yaml
├── project-type.md
├── clutter-report.md          (Wave 1 output)
├── naming-violations.md       (Wave 2 output)
├── structure-proposal.md      (Wave 2 output)
├── expandability-assessment.md (Wave 3 output)
├── execution-plan.md          (Pre-Wave 4)
├── execution-log.md           (Wave 4 output)
└── final-organization-report.md (Phase 6 output)
```

## Persistent Outputs (In Repo Root)

In addition to session-local files, archive-workflow generates persistent files in the target repo:

### .archive-metadata.yaml
- **Location**: Repo root (same level as CLAUDE.md)
- **Purpose**: Machine-readable archival guidelines for consumption by other workflows
- **Content**: Project type, naming conventions summary, structure summary, references to full docs
- **Lifecycle**: Created/overwritten on each archive-workflow run
- **Write pattern**: Atomic (write to .tmp, validate, rename)
- **Schema**: `references/archive-metadata-schema.md`
- **Git**: Committed as part of the archive-workflow Wave 4 commit

### docs/organization/final-organization-report.md
- **Location**: `docs/organization/` directory (created if needed)
- **Purpose**: Human-readable record of what was organized and why
- **Content**: Full organization report from Wave 4
- **Lifecycle**: Created/overwritten on each archive-workflow run
- **Git**: Committed as part of the archive-workflow Wave 4 commit

## Project Type Detection Heuristics

| Signal | Code | Research | Data | Weight |
|--------|------|----------|------|--------|
| package.json, pyproject.toml, Cargo.toml | +++ | - | - | HIGH |
| .ipynb files | + | +++ | + | MEDIUM |
| Large CSV/JSON/Parquet files | - | + | +++ | MEDIUM |
| .tex files, /papers/ directory | - | +++ | - | HIGH |
| src/, tests/ directories | +++ | + | - | HIGH |
| data/, raw/, processed/ | - | + | +++ | HIGH |

**Classification Logic**:
- If max_score > 2 * second_score: Clear winner
- Elif max_score > 1.5 * second_score: Primary + secondary
- Else: Mixed (prompt user to confirm)

## User Confirmation Points

1. **Before Wave 4**: Execution plan approval
2. **Before deletions**: Explicit confirmation per file (Category C)
3. **After completion**: Final report review

## Graceful Cancellation Handling

If user sends SIGINT (Ctrl+C) or requests cancellation:

1. Complete current atomic operation (e.g., single git mv)
2. Update execution-log.md with partial progress
3. Preserve session directory for resume
4. Report cancellation status to user
5. Provide rollback instructions if needed

## Agent Dispatch

Dispatch agents via Task tool. Every dispatch includes the session directory path and
`archival_context: "skip"` — the circular dependency prevention mechanism described
above. Without it, specialists would attempt to check .archive-metadata.yaml during the
very workflow that creates it. Launch the two Wave 2 agents in a single message so they
run concurrently.

| Agent | Reference to load | Reads | Writes |
|-------|-------------------|-------|--------|
| archive-clutter-analyst | `references/clutter-detection-rules.md` | {project_root} | {session_dir}/clutter-report.md (with a Summary section: total items, severity counts, clutter score) |
| archive-nomenclature-enforcer | `references/naming-conventions-{project_type}.md` | {project_root}, clutter-report.md (skip flagged files) | {session_dir}/naming-violations.md (detected patterns, violation severity) |
| archive-structure-organizer | `references/structure-template-{project_type}.md` | {project_root}, clutter-report.md | {session_dir}/structure-proposal.md (migration plan with impact analysis) |
| archive-expandability-reviewer | — | structure-proposal.md | {session_dir}/expandability-assessment.md (flagging issues that should block Wave 4) |
| archive-decision-integrator | Conflict Resolution Rules above | all analyst reports | execution-plan.md, execution-log.md, final-organization-report.md, .archive-metadata.yaml |

archive-decision-integrator categorizes plan operations as A (non-destructive), B
(cleanup), or C (deletions requiring approval), presents the plan for approval, then uses
`git mv` for history preservation and the atomic write pattern for
`.archive-metadata.yaml`.

## Handoffs

| Condition | Hand off to |
|-----------|-------------|
| Wave 1 start | archive-clutter-analyst |
| Wave 2 start | archive-nomenclature-enforcer + archive-structure-organizer (parallel) |
| Wave 3 start | archive-expandability-reviewer |
| Wave 4 start (after approval) | archive-decision-integrator |
| Documentation needed | editor skill |
| Workflow complete | User |
| Critical failure | User (with rollback instructions) |

## References

- references/archival-compliance-check.md (centralized compliance check for all consumers)
- references/archive-metadata-schema.md (.archive-metadata.yaml schema and versioning)
- references/naming-conventions-code.md
- references/naming-conventions-research.md
- references/naming-conventions-data.md
- references/naming-conventions-mixed.md
- references/structure-template-code.md
- references/structure-template-research.md
- references/structure-template-data.md
- references/structure-template-mixed.md
- references/gitignore-patterns.md
- references/clutter-detection-rules.md
- ../references/delegation-and-scope.md (shared orchestrator conventions)

## Examples

- examples/code-project-organization.md
- examples/research-project-organization.md
- examples/mixed-project-organization.md
