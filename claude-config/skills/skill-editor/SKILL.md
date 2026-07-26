---
name: skill-editor
last_updated: 2026-05-24
description: >
  This skill should be used when creating, modifying, or refactoring any Claude Code skill in
  this repository. It is the default entry point for all skill work — new skills from scratch,
  targeted edits, refactoring, and quality review alike. Trigger on "create a skill", "add a
  skill", "build a skill for", "make a skill that", "edit the X skill", "modify a skill",
  "refactor a skill", "improve the X skill", or any time a user wants to add or change a skill
  in this repository. Always prefer this over skill-creator or plugin-dev:skill-development for
  skills that belong in claude-config/.

# Handoff metadata (custom extension -- see workflow-coordinator/references/frontmatter-metadata-standard.md)
handoff:
  accepts_from:
    - "*"
  provides_to:
    - programming-pm
    - technical-pm
  schema_version: "3.0"
  schema_type: universal

categories:
  - skill-development
  - workflow-creation

input_requirements:
  - specification
  - skill-request

output_types:
  - skill
  - agent-configuration
  - documentation
---

# Skill Editor

Workflow for editing Claude Code skills in this repository: scope the change, review it
proportionally to its risk, then implement it through the repo's validate → sync → test →
commit sequence.

## When to Use This Skill

- Creating a new skill in the repository
- Updating, enhancing, or refactoring an existing skill
- Changes spanning multiple files, agents, or architectural decisions
- Skill quality review

## When NOT to Use This Skill

- **Non-skill changes**: modifying agents, settings, or other configuration
- **Exploratory work**: just browsing or understanding skills (use Read or Explore agent)
- **Need a full plugin (multiple skills, agents, commands)**: use `plugin-dev` —
  skill-editor handles individual skill files; plugin-dev covers the broader plugin
  architecture

Note on `skill-creator`: its grader, benchmark aggregator, and eval viewer are useful for
quantitative validation of a newly-created skill. They can be invoked optionally after the
skill is synced to `~/.claude/` — but skill-creator is not an alternative to this workflow
for repo skills.

## Delegation

You coordinate this workflow. Delegate specialist work when the subtask is substantial and
independent — refining an under-specified request through several rounds of user questions,
simulating failure modes for a new multi-agent skill, an adversarial read of a risky plan,
or executing a large multi-file change. Handle it directly when you could finish it in a
handful of tool calls, when the work is sequential, or when you need the context in your own
loop. See `../references/delegation-and-scope.md`.

Concretely: a typo fix, a description tweak, adding one example, or a single-file edit under
about 50 lines does not need any of the four specialist agents. Read the file, make the
change, and run the Phase 4 sequence. Reserve the panel for changes where an independent
perspective has something to find.

You own: mode and scope decisions, session state, gate evaluation, user communication, and
workflow routing.

## State Anchoring

For multi-phase runs, start each response with `[Phase N/4 - {phase_name}] {brief status}`,
e.g. `[Phase 3/4 - Decision] Synthesizing analysis into a plan`. The authoritative record
is `${SESSION_DIR}/session-state.json`; trust it over memory.

## Choosing How Much Review to Run

Match the review to the change. The four phases below are named stages, not a mandatory
gauntlet — run the ones that earn their cost.

| Change | Review | Phases |
|--------|--------|--------|
| Typo, comment, doc wording, one example, single file under ~50 lines | None | Edit directly, then Phase 4 |
| Substantive edit to one skill: new section, reworked workflow step, changed triggers | Light — refine scope with the user if the request is ambiguous, plan, implement | 1 (in-loop or via request-refiner), 3, 4 |
| New skill from scratch; multi-file or multi-agent change; changes to workflow phases, quality gates, or agent wiring; anything the user asks to be reviewed | Full panel | 1, 2, 3, 4 |

Scale signals for the middle and top rows are in
`references/complexity-detection-criteria.md`. When genuinely uncertain, ask the user rather
than defaulting to the heaviest path.

Record the chosen scope in session state.

## Pre-Workflow: Safety Checks

Run the session management script, which performs:
- Git safety checks (uncommitted changes, merge/rebase detection, detached HEAD)
- sync-config.py status verification
- Directory verification (must be repo root)
- Archival awareness detection
- Trap handler registration for graceful interrupt
- Session management commands (`--list-sessions`, `--cleanup`)
- Resume protocol with multi-session support (including legacy format migration)
- Session directory creation and state initialization

**Implementation**: see `references/session-management.sh` for the complete bash.

If checks fail: ask the user to resolve before proceeding.

### If User Cancels (Ctrl+C)

Session state is preserved in `${SESSION_DIR}/session-state.json`. On next invocation,
offer to resume from the last phase; if declined, the session remains in
/tmp/skill-editor-session/{session-id}. Re-sync if needed: `./sync-config.py push`.

## Phase 1: Refinement

**Objective**: turn an under-specified request into a clear specification.

**Agent**: `skill-editor-request-refiner` — use it when the request is genuinely ambiguous
or needs several rounds of clarifying questions. For a clear, narrow request, ask the one
or two questions you need with AskUserQuestion and write the spec yourself.

The specification covers: objective (one sentence), scope (IN/OUT), measurable success
criteria, files affected, and user approval. Written to
`${SESSION_DIR}/refined-specification.md`.

**Gate 1**: the user approves the specification. If not, refine further.

### Orchestrator Detection

If the target skill coordinates other skills or agents, orchestrator-specific patterns
apply during analysis and planning. Score the target with the signals in
`references/orchestrator-best-practices.md` (name keywords, coordination terms in the
description, Task tool delegation, named phases, quality gates, session state). A score of
4 or more suggests an orchestrator; 2-3 is ambiguous, so ask. For a new skill, ask directly
whether it will orchestrate other skills.

Record the result in session state:
```json
"orchestrator_detected": true/false,
"orchestrator_confidence": "high"/"medium"/"none",
"orchestrator_user_confirmed": true/false
```

## Phase 2: Analysis (full panel only)

**Objective**: analyze the proposed change from perspectives you don't already hold.

Two independent tracks. Launch them in a single message so they run concurrently — see
`../references/delegation-and-scope.md` on keeping fan-out proportional.

**Track A — `brainstorming-pm`** (multi-perspective swarm). Fill the template in
`references/swarm-challenge-templates.md` with skill name, change type, file and line
counts, orchestrator detection result, specification summary, and current skill summary.
Useful output contains at least a couple of specific recommendations and one alternative
approach; if it returns less than that, note it and proceed without swarm input. When
orchestrator_detected is true, add to the template: "This skill is an orchestrator --
evaluate against orchestrator patterns in references/orchestrator-checklist.md (6 core + 4
common patterns)".

**Track B — `skill-editor-edge-case-simulator`** (domain failure modes). Produces a
skill-specific failure mode matrix: YAML parsing failures, sync-config.py edge cases, Task
tool failures, git dirty state, and other Claude Code boundary conditions.

**Outputs**: `${SESSION_DIR}/swarm-synthesis.md` (absent if skipped) and
`${SESSION_DIR}/edge-cases.md`.

**Gate 2**: at least the edge-case report is present and substantive. Missing or thin swarm
output degrades gracefully — proceed with what you have and note the gap in the plan. If
the edge-case simulator returns nothing usable, retry once, then ask the user whether to
proceed on the specification alone.

## Phase 3: Decision

Synthesize what you have into an implementation plan. Read the swarm synthesis and
edge-case report if they exist, identify where they agree and where they conflict, and
resolve the conflicts — asking the user via AskUserQuestion for decisions that add agents
or change structure, and deciding minor ones (examples, docs) yourself.

The plan states exact file paths and specific changes, how the identified edge cases are
handled, the validation steps, and the rollback plan. Written to
`${SESSION_DIR}/implementation-plan.md`.

When orchestrator_detected is true, note absent orchestrator patterns from
`references/orchestrator-best-practices.md` and say whether the plan adopts them.

For a light-review change, this is a short plan written directly from the specification: the
objective, files to modify, validation steps, and rollback. If the target files include core
workflow or agent files, escalate to the full panel first.

### Adversarial Review

**Agent**: `skill-editor-adversarial-reviewer`. Worth its cost when the plan touches
workflow structure, agent wiring, multiple files, or anything destructive. It reads the plan
skeptically, challenges assumptions, checks the file paths and git workflow, verifies
alignment with the specification, and returns GO / CONDITIONAL / NO-GO to
`${SESSION_DIR}/adversarial-review.md`.

**Gate 2 (plan approval)**: the plan has exact file paths, the git workflow is safe, the
reviewer (if run) returned GO or a CONDITIONAL whose fixes are applied, and the user
approves. On NO-GO, revise and resubmit.

## Phase 4: Execution (Implement + Validate + Commit)

This sequence is destructive if done wrong — it overwrites `~/.claude/`. Follow the steps
as written. `skill-editor-executor` runs them for large or multi-file plans; for a small
edit, run them yourself.

#### Step 1: Pre-Implementation Safety
- `git status` — must be clean
- `./sync-config.py status` — must be synced
- `pwd` — must be repo root
- Stop if any check fails.

#### Step 2: Implement Changes
For each file in the plan:
- **Edit**: Read first, then Edit with exact string replacement
- **Create**: Write new file
- **Delete**: Remove file

#### Step 3: Quality Gate 3 - Pre-Sync Validation
- Validate YAML frontmatter for all modified skills
- Validate JSON for modified agents
- Dry-run sync: `./sync-config.py push --dry-run`

**If Gate 3 fails**: fix the issues and re-validate. Do not proceed until it passes — the
push-time gates exist because a malformed skill breaks the live config.

#### Step 4: Sync to ~/.claude/
- `./sync-config.py push`
- `./sync-config.py status` — verify no divergence

#### Step 5: Test Skill Invocation
- Verify skill file exists at `$HOME/.claude/skills/$SKILL_NAME/SKILL.md`
- Verify YAML parses
- Smoke test existing skills for regressions

If any of these objective checks fail: roll back via `git reset --hard HEAD`, re-sync, fix,
retry.

#### Step 6: Update Planning Journal
`./sync-config.py plan --title "[Brief description from refined spec]"`

Optionally run `claude-md-management:revise-claude-md` to incorporate workflow learnings
into CLAUDE.md documentation.

#### Optional: Git Strategy Advisory
Before committing, `git-strategy-advisor` can be invoked via Task tool in post-work mode for
scope-adaptive git recommendations. Advisory only — Step 7 takes precedence.

#### Step 7: Commit Changes

Stage specific files (never `-A` or `.`), commit with a HEREDOC multi-line message using
conventional commit format (`feat`/`fix`/`docs`) and including
`Co-Authored-By: Claude <noreply@anthropic.com>`, then mark the session completed.
Never bypass hooks and never force push.

#### Step 8: Report Completion

Summary of changes, validation and test results, commit SHA, planning journal entry path,
and session status.

## Escalation

From CONFIG_MANAGEMENT.md decision thresholds:

- **Ask the user (AskUserQuestion)** before adding a new agent to a workflow, changing the
  skill structure specification, modifying core workflow phases, changing an existing
  skill's core workflow, adding a new supporting skill, or changing a naming convention.
- **Decide yourself, then report**: adding an example, fixing a typo, updating reference
  material.

## Error Handling

**Retry and degradation** (event-based, evaluated when an agent returns):
- edge-case-simulator returns nothing usable: retry once, then ask the user
- brainstorming-pm returns nothing usable or stalls: proceed with the edge-case output only
- Global ceiling: if two or more components fail, stop and escalate rather than retrying
  further

**Rollback (Phase 4 failures)**:
1. Stop immediately
2. `git reset --hard HEAD` (revert uncommitted changes)
3. `./sync-config.py push` (re-sync from repo)
4. Document failure in planning journal
5. Report to user with options: retry, skip, or abort

**Interrupt (user cancels)**: check git status, `git reset --hard HEAD`, re-sync with
`./sync-config.py push`, leave session state in `${SESSION_DIR}/` for resume, and record
"Cancelled by user" in the planning journal.

## Integration with Existing Tools

### CONFIG_MANAGEMENT.md

This workflow extends the 7-step CONFIG_MANAGEMENT.md process:

- **Step 1 (Safety Check)**: Pre-workflow checks
- **Step 2 (Planning Entry)**: Phase 4, Step 6
- **Step 3 (Implement)**: Phase 4, Step 2
- **Step 4 (Quality Analysis)**: Phases 2-3
- **Step 5 (Preview/Sync)**: Phase 4, Steps 3-4
- **Step 6 (Test)**: Phase 4, Step 5
- **Step 7 (Commit)**: Phase 4, Step 7

### sync-config.py

- `./sync-config.py status` (pre-flight check)
- `./sync-config.py push --dry-run` (validation)
- `./sync-config.py push` (apply changes)
- `./sync-config.py plan` (create planning entry)

### Planning Journal

Entry created in Phase 4: title (brief description from the refined spec), objective,
files modified, validation and test results, outcome.

## Quality Gates

| Gate | Phase | Criteria | Failure Action |
|------|-------|----------|----------------|
| 1: Specification Approval | Phase 1 | Spec approved by user | Return to refinement |
| 2: Plan Approval | Phase 3 | Plan has exact paths; adversarial GO (if run); user approves | Revise plan |
| 3: Execution Verification | Phase 4 | YAML validates, sync succeeds, skill invokes, no regressions | Rollback |

Gate 3 is objective — validators, sync status, and invocation, not a prose self-review.

## Examples

**Small edit (no panel)**: "Fix the typo in researcher's description." Read the file, make
the edit, validate YAML, `push --dry-run`, `push`, verify, planning entry, commit.

**Substantive edit (light review)**: "Add parallel web search to the researcher skill."
Clarify scope (which phase, how many searches, how results merge), write a short plan naming
`claude-config/skills/researcher/SKILL.md` and the lines to change, implement, then Phase 4.

**New skill (full panel)**: "Create a skill for API documentation." Refine requirements
(which APIs, output format, tools) via request-refiner; run the swarm and edge-case
simulator in parallel; synthesize a plan covering file structure, workflow, and examples;
adversarial review; then execute Phase 4.

## References

See `skill-editor/references/` for:
- `swarm-challenge-templates.md`: challenge template for brainstorming-pm swarm delegation
- `session-management.sh`: git safety checks, session creation/resume, cleanup commands
- `complexity-detection-criteria.md`: change-scale signals for choosing review depth
- `anthropic-guidelines-summary.md`: Anthropic best practices
- `skill-structure-specification.md`: skill format and validation
- `quality-gates.md`: gate criteria and validation commands
- `orchestrator-checklist.md`: orchestrator pattern evaluation checklist
- `orchestrator-best-practices.md`: orchestrator pattern catalogue and templates

Shared conventions: `../references/delegation-and-scope.md`.
