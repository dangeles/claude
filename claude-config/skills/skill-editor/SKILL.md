---
name: skill-editor
description: Use when creating, modifying, or refactoring Claude Code skills that require structured multi-agent review and quality validation

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

Comprehensive multi-agent workflow system for editing Claude Code skills with structured phases, quality gates, and expert review.

## When to Use This Skill

Use this skill when:

1. **Creating new skills**: User wants to add a new skill to the repository
2. **Modifying existing skills**: User wants to update, enhance, or refactor a skill
3. **Complex skill changes**: Change involves multiple files, agents, or architectural decisions
4. **Quality assurance needed**: Change requires thorough review and validation

This skill provides:
- Structured 4-phase workflow
- Interactive requirements refinement
- Parallel expert analysis (4 simultaneous agents)
- Adversarial review before implementation
- Automated validation and testing
- Integration with sync-config.py and planning journal

## When NOT to Use This Skill

Do NOT use this skill when:

- **Simple documentation fixes**: Typo fixes, minor documentation updates (edit directly)
- **Non-skill changes**: Modifying agents, settings, or other configuration
- **Urgent hotfixes**: Emergency fixes that can't wait for full workflow
- **Exploratory work**: Just browsing or understanding skills (use Read or Explore agent)

## Delegation Mandate

You are an **orchestrator**. You coordinate specialist agents -- you do not perform specialist analysis, research, or implementation yourself.

**You ARE the coordinator who ensures** analysis, research, review, and implementation happen through delegation to specialist agents via Task tool.

**You are NOT** an analyst, researcher, reviewer, or implementor. You do not perform best-practices analysis, external research, edge-case simulation, knowledge-engineering analysis, adversarial review, or code implementation yourself.

**Orchestrator-owned tasks** (you DO perform these yourself):
- Session setup, directory creation, state file management
- Mode detection and user interaction for mode selection
- Quality gate evaluation (checking that agent outputs meet criteria)
- User communication (presenting options, gathering decisions)
- Workflow routing (determining which phase to execute next)
- Pre-flight validation (git checks, file existence)
- Orchestrator detection (determining if target skill is an orchestrator)

### When You Might Be Resisting Delegation

| Rationalization | Reality |
|----------------|---------|
| "This analysis is too simple to delegate" | Simple tasks still consume context window. Delegate. |
| "I can do it faster myself" | Speed is not the goal; context isolation and specialist quality are. |
| "The agent will just repeat what I already know" | The agent provides independent verification. Your knowledge may be incomplete. |
| "It's just a quick read of the file" | Reading specialist content to make specialist decisions IS specialist work. |

**Self-check before every action**: "Am I about to load specialist instructions into my context so I can do their work? If yes, use Task tool instead."

## State Anchoring

Start every response with a phase indicator:

```
[Phase N/4 - {phase_name}] {brief status}
```

Examples:
- `[Phase 1/4 - Refinement] Gathering requirements from user`
- `[Phase 2/4 - Analysis] 3/4 agents completed, waiting for external-researcher`
- `[Phase 3/4 - Decision] Synthesizing 4 analysis reports`
- `[Phase 4/4 - Execution] Implementing change 3/12`

**Protocol**:
1. Before starting any phase: Read `${SESSION_DIR}/session-state.json`. Confirm current_phase matches expectations.
2. After any user interaction: Answer the user, then re-anchor with phase indicator.
3. If phase indicator and state file disagree: Trust state file, not memory.

## Tool Selection

| Situation | Tool | Reason |
|-----------|------|--------|
| Phase 2 parallel analysis (4 agents) | Task tool | Context isolation, parallel execution |
| Phase 2.5 strategic review | Task tool | Separate specialist context |
| Phase 3 synthesis | Task tool | Independent decision-making context |
| Phase 3 adversarial review | Task tool | Independent, skeptical review |
| Phase 4 implementation | Task tool | Isolated execution environment |
| Loading reference docs for YOUR routing decisions | Read tool | Orchestrator decision support |
| Loading skill instructions to decide WHICH specialist to invoke | Read tool (brief scan) | Routing information, not specialist work |
| User interaction (questions, approvals, options) | AskUserQuestion | Structured user communication |
| File operations (create, modify files) | Write tool (via executor agent) | Delegated to executor specialist |
| Validation scripts, git operations | Bash tool | Infrastructure commands |

**Self-check**: "Am I about to load specialist instructions into my context so I can do their work? If yes, use Task tool instead."

## Workflow Overview

```
SIMPLE MODE (15-45 min)
├── Phase 1: Refinement (5-15 min)
├── Mode Selection: User confirms SIMPLE
├── [SKIP Phase 2: No parallel analysis]
├── [SKIP Phase 2.5: No strategic review]
├── Phase 3: Lightweight Decision (10-20 min)
│   └── Minimal synthesis from specification only
└── Phase 4: Execution (10-20 min)
    └── Gates 4 & 5 always run

STANDARD MODE (1.5-3 hours) [Current default]
├── Phase 1: Refinement (10-30 min)
├── Mode Selection: User confirms STANDARD
├── Phase 2: Parallel Analysis (30-60 min, 4 agents)
├── [Phase 2.5: Strategic Review - conditional, stricter triggers]
├── Phase 3: Decision & Review (45-90 min)
│   └── Full synthesis + adversarial review
└── Phase 4: Execution (60-120 min)
    └── Gates 4 & 5 always run

EXPERIMENTAL MODE (10-30 min) [User-requested]
├── Phase 1: Quick Refinement (5-10 min)
├── Mode Selection: User explicitly requests EXPERIMENTAL
├── [SKIP Phase 2]
├── [SKIP Phase 2.5]
├── Phase 3: Minimal Decision (5-10 min)
│   └── Direct implementation plan with experimental tags
└── Phase 4: Execution with rollback plan (5-15 min)
    └── Gates 4 & 5 always run + experimental tagging
```

## Workflow

### Pre-Workflow: Safety Checks

Before starting workflow, run the session management script which performs:
- Git safety checks (uncommitted changes, merge/rebase detection, detached HEAD)
- sync-config.py status verification
- Directory verification (must be repo root)
- Archival awareness detection
- Trap handler registration for graceful interrupt
- Session management commands (`--list-sessions`, `--cleanup`)
- Resume protocol with multi-session support (including legacy format migration)
- Session directory creation and state initialization

**Implementation**: See `references/session-management.sh` for complete bash.

If checks fail: Ask user to resolve before proceeding.

### If User Cancels (Ctrl+C)

Session state is preserved in `${SESSION_DIR}/session-state.json`.

On next invocation:
1. Offer to resume from last phase
2. If declined, session remains in /tmp/skill-editor-session/{session-id}
3. Re-sync if needed: `./sync-config.py push`

### Phase 1: Refinement (Interactive)

**Objective**: Transform user's request into detailed, unambiguous specification.

**Agent**: `skill-editor-request-refiner`

**Model**: Opus 4.6

**Process**:

1. Launch request-refiner agent via Task tool
2. Agent asks clarifying questions to understand:
   - What user wants to change
   - Why they want this change
   - What success looks like
   - What's in scope vs. out of scope
3. Agent reads existing skill (if modifying)
4. Agent establishes clear boundaries and success criteria
5. Agent presents refined specification to user

**Output File**: `${SESSION_DIR}/refined-specification.md` containing:
- Objective (one sentence)
- Scope (IN/OUT lists)
- Success criteria (measurable)
- Files affected
- User approval

**Quality Gate 1: Specification Approval**

User must approve:
- [ ] Specification matches intent
- [ ] Scope is appropriate
- [ ] Success criteria are clear
- [ ] Ready to proceed to analysis

**If Gate 1 fails**: Return to request-refiner for more refinement.

**If Gate 1 passes**: Update session state and proceed to Mode Selection.

**Post-Gate 1: Orchestrator Detection**

After specification approval, determine if the target skill is an orchestrator:

1. Read target SKILL.md (if editing an existing skill)
2. Score against detection criteria:

| Signal | Score | Check |
|--------|-------|-------|
| Name contains orchestrator keyword (pm, coordinator, orchestrator, pipeline, architect) | +1 | Check skill name |
| Description mentions coordination terms (coordinate, orchestrate, multi-agent, pipeline) | +1 | Check description field |
| Delegates to other skills via Task tool | +2 | Search for Task tool usage |
| Has named phases/stages with sequential progression | +1 | Search for Phase/Stage headers |
| Has quality gates between phases | +1 | Search for Gate references |
| Manages session state across phases | +1 | Search for state file management |

3. Apply thresholds:
   - Score >= 4: `"Detected as orchestrator (confidence: high). Apply orchestrator analysis? [Y/n]"`
   - Score 2-3: `"May be an orchestrator (confidence: medium). Apply orchestrator analysis? [y/N]"`
   - Score 0-1: Not an orchestrator. Skip orchestrator analysis.
   - Always append: `"If this IS an orchestrator, reply 'orchestrator' to enable pattern analysis."`

4. If creating a new skill: Ask directly: "Will this skill orchestrate other skills? [y/N]"

5. Record detection result in session state:
   ```json
   "orchestrator_detected": true/false,
   "orchestrator_confidence": "high"/"medium"/"none",
   "orchestrator_user_confirmed": true/false
   ```

Update session state to phase 1.5 with `agents_completed: ["request-refiner"]`.

---

### Mode Selection (After Phase 1)

**Objective**: Select workflow execution mode based on complexity detection and user preference.

**Trigger**: After Quality Gate 1 passes (specification approved)

Run the mode detection script which performs:
- Three-tier complexity detection (COMPLEX / SIMPLE / STANDARD / EXPERIMENTAL)
- Metrics extraction from specification (files changed, lines changed, scope keywords)
- Fail-safe default to STANDARD when uncertain
- Mode selection display with user prompt (A/B/C for SIMPLE/STANDARD/EXPERIMENTAL)
- User override confirmation for risky downgrades
- Experimental mode warning and acknowledgment
- Mode recording to session state (`mode-selection.json`)
- Mode-based branching to appropriate phase

**Implementation**: See `references/mode-detection.sh` for complete bash.

**Detection logic summary**:
- **COMPLEX** (triggers Phase 2.5): New skill, >4 files, >300 lines, explicit request, refactoring with moderate+ scope
- **SIMPLE**: Documentation/typo/comment with <=1 file and <=50 lines
- **STANDARD**: Default for everything else
- **EXPERIMENTAL**: User keyword triggers (experimental, quick, try, prototype)

---

### Phase 2: Parallel Analysis (4 Simultaneous Agents)

**Objective**: Analyze proposed change from multiple expert perspectives.

**Agents** (all run in parallel):
1. `skill-editor-best-practices-reviewer` (Opus 4.6) - Critical
2. `skill-editor-external-researcher` (Opus 4.6) - Supplementary
3. `skill-editor-edge-case-simulator` (Opus 4.6) - Critical
4. `skill-editor-knowledge-engineer` (Opus 4.6) - Critical [NEW]

**Process**:

Launch all 4 agents with wave-based execution to reduce resource contention:

**Wave 1 (T=0s)**: Launch critical analysis agents (best-practices-reviewer, edge-case-simulator)
**Wave 2 (T=30s)**: Launch structural analysis agent (knowledge-engineer)
**Wave 3 (T=60s)**: Launch supplementary research agent (external-researcher)

**Rationale for wave-based execution**: Staggering launches by 30-60 seconds reduces system resource contention and improves reliability for parallel agent execution.

**Important**: All 4 agents run in parallel (waves overlap). Wait for all to complete before proceeding to Phase 3.

**Orchestrator Analysis** (conditional -- only when orchestrator_detected is true in session state):

When the target skill is an orchestrator, Phase 2 agents perform additional analysis:

- **best-practices-reviewer**: Evaluates 6 REQUIRED patterns from `references/orchestrator-checklist.md`
- **knowledge-engineer**: Evaluates 4 RECOMMENDED patterns from `references/orchestrator-checklist.md`
- **Neither agent evaluates all 11 patterns.** Division of labor prevents cognitive overload.
- **external-researcher and edge-case-simulator**: No additional orchestrator-specific tasks.

**Agent Timeouts and Retry Logic**:

| Agent Type | Timeout | On Failure |
|-----------|---------|------------|
| Critical (best-practices, edge-case, knowledge) | 10 min | Auto-retry once (30s wait), then ask user |
| Supplementary (external-researcher) | 10 min | Proceed without |

Maximum 2 attempts per critical agent.

**Output Files** (must be created before proceeding to Phase 3):
- `${SESSION_DIR}/best-practices-review.md`
- `${SESSION_DIR}/external-research.md`
- `${SESSION_DIR}/edge-cases.md`
- `${SESSION_DIR}/knowledge-engineering-analysis.md`

**Quality Gate 2: Analysis Completion**

| Completed Agents | Critical Agents Status | Action |
|------------------|----------------------|--------|
| 4/4 | All critical complete | PASS - Proceed to Phase 3 |
| 3/4 | All critical (only external-researcher failed) | PASS - Proceed with note |
| 3/4 | 1 critical failed (first attempt) | RETRY - Retry failed critical agent once |
| 3/4 | 1 critical failed (after retry) | ASK USER - Proceed with placeholder or abort? |
| 2/4 or fewer | Multiple critical failed | FAIL - Retry all failed critical agents or abort |

**If Gate 2 passes**: Update session state to phase 3 and proceed.

---

### Phase 2.5: STRATEGIC REVIEW [CONDITIONAL]

**Purpose**: Strategic architectural assessment for complex changes using cross-domain pattern matching.

**When**: Conditionally executed based on complexity detection. Skipped for simple changes.

**Agent**: strategy-consultant (Opus 4.6)

**Process**:
1. Run complexity detection to determine if strategic review is needed
2. If triggered: Launch strategy-consultant agent with 30-minute timeout
3. Validate strategic review quality (>200 words, patterns identified, recommendations present)
4. Check for major refactoring opportunities (user decision: proceed/explore/abort)
5. Quality Gate 2.5: Verify complexity detection completed, review exists if needed, user decisions recorded

**Implementation**: See `references/phase-2-5-detection.sh` for complete bash including complexity detection, strategy consultant launch/timeout handling, quality validation, major refactoring detection, and Quality Gate 2.5.

**Note**: Phase 2.5 is optional and conditional. If skipped, strategic-review.md will not exist, and Phase 3 agents handle this gracefully.

---

### Phase 3 Variants by Mode

#### Phase 3: SIMPLE MODE (Lightweight Decision)

**Duration**: 10-20 minutes
**Trigger**: SELECTED_MODE = "SIMPLE"

Process:
1. Create minimal implementation plan from specification (objective, files to modify, validation steps, rollback plan)
2. Check if target files include core workflow/agent files — if so, offer upgrade to Standard Mode
3. Skip adversarial review unless core files affected

#### Phase 3: EXPERIMENTAL MODE (Minimal Decision)

**Duration**: 5-10 minutes
**Trigger**: SELECTED_MODE = "EXPERIMENTAL"

Process:
1. Create implementation plan with experimental flags and WARNING header
2. Optionally run adversarial review (user choice, adds ~15 min)
3. All output tagged as experimental/not production-ready

#### Phase 3 Mode Checkpoint

Before launching synthesis, offer mode change option. If currently in SIMPLE or EXPERIMENTAL, user can switch to STANDARD (which will run Phase 2, adding ~1.5 hours).

---

### Phase 3: Decision & Review (Synthesis + Adversarial)

**Objective**: Synthesize analyses, make decisions, create plan, get expert approval.

#### Part A: Decision Synthesis

**Agent**: `skill-editor-decision-synthesizer`

**Model**: Opus 4.6 (critical decision-making)

**Process**:

1. Read all 4 analysis reports + refined specification
2. Identify consensus and conflicts
3. Resolve conflicts or present options to user:
   - **Major decisions**: MUST ask user (new agents, structure changes)
   - **Medium decisions**: SHOULD ask user (workflow changes)
   - **Minor decisions**: Agent decides (examples, docs)
4. Create detailed implementation plan with:
   - Exact file paths
   - Specific changes (line numbers if possible)
   - Edge case handling
   - Git workflow
   - Validation steps
   - Rollback plan

**Orchestrator-Specific Synthesis** (when orchestrator_detected is true in session state):

5. Read orchestrator checklist results from both best-practices-review.md and knowledge-engineering-analysis.md
6. If REQUIRED patterns are ABSENT: Implementation plan MUST include adding those patterns
7. If RECOMMENDED patterns are ABSENT: Implementation plan SHOULD note them as suggested additions
8. Reference pattern templates from `orchestrator-best-practices.md` for copy-paste inclusion
9. Check Pattern Interactions section to avoid contradictions
10. For existing orchestrators: PARTIAL with a working variant is acceptable

**Output File**: `${SESSION_DIR}/implementation-plan.md`

#### Part B: Adversarial Review

**Agent**: `skill-editor-adversarial-reviewer`

**Model**: Opus 4.6 (expert review)

**Process**:

1. Read implementation plan with expert skepticism
2. Challenge assumptions and approach
3. Identify failure modes not caught by analysis
4. Verify exact file paths (run bash checks)
5. Verify git workflow safety
6. Check alignment with original specification
7. Provide go/no-go decision

**Output File**: `${SESSION_DIR}/adversarial-review.md` containing:
- Architecture assessment
- Failure mode analysis
- Integration risk assessment
- Exact file path verification
- Git workflow verification
- Final decision: GO / CONDITIONAL / NO-GO

**Quality Gate 3: Plan Approval**

Check:
- [ ] Implementation plan has exact file paths
- [ ] Git workflow is safe and correct
- [ ] Integration points identified
- [ ] No architectural concerns
- [ ] Adversarial reviewer approved (GO or CONDITIONAL with fixes applied)
- [ ] User approves plan

**If Gate 3 fails**:
- If CONDITIONAL: Fix issues, re-review
- If NO-GO: Return to decision-synthesizer, revise plan
- If user doesn't approve: Refine plan or return to Phase 1

**If Gate 3 passes**: Update session state to phase 4 and proceed.

---

### Phase 4: Execution (Implement + Validate + Commit)

**Objective**: Execute approved plan with validation at each step.

**Agent**: `skill-editor-executor`

**Model**: Opus 4.6

**Process**:

#### Step 1: Pre-Implementation Safety
- `git status` — must be clean
- `./sync-config.py status` — must be synced
- `pwd` — must be repo root
- Stop if any check fails.

#### Step 2: Implement Changes
For each file in implementation plan:
- **Edit**: Read first, then Edit with exact string replacement
- **Create**: Write new file
- **Delete**: Remove file

#### Experimental Mode Output Tagging
If workflow_mode = "EXPERIMENTAL", add experimental tags to output skill files (YAML frontmatter `experimental: true` + warning comment).

**Implementation**: See `references/experimental-tagging.sh` for complete bash.

#### Step 3: Quality Gate 4 - Pre-Sync Validation
- Validate YAML frontmatter for all modified skills
- Validate JSON for modified agents
- Dry-run sync: `./sync-config.py push --dry-run`

**If Gate 4 fails**: Fix issues, re-validate, do NOT proceed until pass.

#### Step 4: Sync to ~/.claude/
- `./sync-config.py push`
- `./sync-config.py status` — verify no divergence

#### Step 5: Test Skill Invocation
- Verify skill file exists at `$HOME/.claude/skills/$SKILL_NAME/SKILL.md`
- Verify YAML parses
- Smoke test existing skills (no regressions)

#### Step 6: Quality Gate 5 - Post-Execution Verification
- [ ] Original requirement met (from refined spec)
- [ ] Edge cases handled (from edge-case report)
- [ ] sync-config.py push successful
- [ ] Skill invokes without errors
- [ ] No regressions in existing skills
- [ ] Planning journal entry ready

**If Gate 5 fails**: Rollback via `git reset --hard HEAD`, re-sync, fix, retry.

#### Step 7: Update Planning Journal
`./sync-config.py plan --title "[Brief description from refined spec]"`

#### Optional: Git Strategy Advisory
Before committing, MAY invoke `git-strategy-advisor` via Task tool in post-work mode for scope-adaptive git recommendations. This is advisory only — Step 8 logic takes precedence.

#### Step 8: Commit Changes

Determine commit prefix based on mode:
- EXPERIMENTAL: `experimental` prefix + `[EXPERIMENTAL - requires full review]` suffix
- Others: `feat` prefix

Stage specific files (NEVER -A or .), commit with HEREDOC multi-line message including `Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>`, then mark session as completed.

**Git Safety Checklist**:
- [ ] Specific files staged (not -A or .)
- [ ] Conventional commit format (feat/fix/docs)
- [ ] Descriptive message
- [ ] Co-authored-by line
- [ ] No destructive operations
- [ ] No hook bypasses

#### Step 9: Report Completion

Generate completion report with:
- Summary of changes
- Validation results (Gates 4 & 5)
- Testing results
- Commit SHA
- Planning journal entry path
- Success criteria verification
- Session completion status

---

## Escalation Framework

Decision thresholds (from CONFIG_MANAGEMENT.md):

### Major Decisions → User Approval Required

- Add new agent to workflow
- Change skill structure specification
- Modify core workflow phases

**Action**: Use AskUserQuestion before proceeding

### Medium Decisions → User Approval Required

- Modify existing skill's core workflow
- Add new supporting skill
- Change skill naming convention

**Action**: Use AskUserQuestion with options

### Minor Decisions → Agent Decides

- Add example to existing skill
- Fix documentation typo
- Update reference material

**Action**: Proceed, notify user

## Error Handling

### Retry Protocol (Phase 2 Agent Failures)
- First failure: Wait 30s, retry automatically
- Second failure: User decision required (proceed with placeholder or abort)
- Maximum 2 attempts per critical agent
- Retried operations should be idempotent

### Graceful Degradation (Supplementary Agent Failures)
- external-researcher timeout: Proceed without research analysis
- knowledge-engineer timeout (after retry): Proceed with 3 analyses
- Decision-synthesizer notes missing perspectives in synthesis

### Circuit Breaker (Cascading Failures)
- If 2+ critical agents fail in Phase 2: Stop retrying, escalate to user
- User choices: retry all, proceed with available, or abort

### Rollback Protocol (Phase 4 Failures)
1. Stop immediately
2. `git reset --hard HEAD` (revert uncommitted changes)
3. `./sync-config.py push` (re-sync from repo)
4. Document failure in planning journal
5. Report to user with options: retry, skip, or abort

### Interrupt Handling (User Cancels)
1. Check git status
2. Rollback uncommitted changes: `git reset --hard HEAD`
3. Re-sync: `./sync-config.py push`
4. Session state preserved in `${SESSION_DIR}/` for potential resume
5. Document in planning journal: "Cancelled by user"

## Integration with Existing Tools

### CONFIG_MANAGEMENT.md

This workflow extends the 7-step CONFIG_MANAGEMENT.md process:

- **Step 1 (Safety Check)**: Pre-workflow checks
- **Step 2 (Planning Entry)**: Phase 4, Step 7
- **Step 3 (Implement)**: Phase 4, Step 2
- **Step 4 (Quality Analysis)**: Phases 2-3, Quality Gates
- **Step 5 (Preview/Sync)**: Phase 4, Steps 3-4
- **Step 6 (Test)**: Phase 4, Step 5
- **Step 7 (Commit)**: Phase 4, Step 8

### sync-config.py

Executor agent uses sync-config.py:
- `./sync-config.py status` (pre-flight check)
- `./sync-config.py push --dry-run` (validation)
- `./sync-config.py push` (apply changes)
- `./sync-config.py plan` (create planning entry)

### Planning Journal

Planning entry created in Phase 4, Step 7:
- Title: Brief description from refined spec
- Objective: From refined specification
- Changes: Files modified
- Testing: Validation and test results
- Outcome: Success/Partial/Failed

## Quality Gates Summary

| Gate | Phase | Owner | Criteria | Failure Action |
|------|-------|-------|----------|----------------|
| 1 | Phase 1 | request-refiner | Spec approved | Return to refinement |
| 2 | Phase 2 | decision-synthesizer | All analyses complete | Re-run agents |
| 3 | Phase 3 | adversarial-reviewer | Plan approved | Revise plan |
| 4 | Phase 4 | executor | Syntax validated | Fix issues |
| 5 | Phase 4 | executor | Implementation verified | Rollback |

## Examples

### Example 1: Add Parallel Execution to Researcher

**User Request**:
```
/skill-editor "Add parallel web search to researcher skill"
```

**Phase 1 Output**:
```markdown
Objective: Modify researcher skill to execute 3 WebSearch calls in parallel

Scope:
- IN: researcher/SKILL.md Phase 2 workflow
- OUT: No changes to agents or other phases

Success Criteria:
- 3 WebSearch calls execute simultaneously
- Results synthesized correctly
- No regressions
```

**Phase 2 Findings**:
- Best practices: Use Task tool for parallel calls
- Research: Community uses this pattern
- Edge cases: Handle timeout, network failure

**Phase 3 Plan**:
```markdown
Edit: claude-config/skills/researcher/SKILL.md
Lines 45-60: Replace sequential WebSearch with parallel

Implementation:
[3 Task tool calls in single message]
```

**Phase 4 Result**:
```
YAML validates, sync succeeds, skill invokes correctly
Commit: feat(researcher): Add parallel web search
```

### Example 2: Create New Skill

**User Request**:
```
/skill-editor "Create a new skill for API documentation"
```

**Process**:
- Phase 1: Refine requirements (which APIs? format? tools?)
- Phase 2: Analyze (best practices for doc skills, community patterns, edge cases)
- Phase 3: Plan (file structure, workflow steps, examples)
- Phase 4: Create files, validate, sync, test, commit

## Timeout Configuration

| Phase | Component | Timeout | Exceeded Action |
|-------|-----------|---------|-----------------|
| 1 | request-refiner | 30 min | Escalate to user |
| 2 | critical agents (Wave 1-2) | 10 min each | Auto-retry once, then user decision |
| 2 | supplementary agent (Wave 3) | 10 min | Proceed without |
| 2.5 | strategy-consultant | 30 min | User decision (proceed/retry/abort) |
| 3 | decision-synthesizer | 30 min | Escalate to user |
| 3 | adversarial-reviewer | 30 min | Escalate to user |
| 4 | executor | 60 min | Escalate to user |
| Global | entire workflow | 4 hours | Safety ceiling, force escalate |

## Notes

- **Parallel execution in Phase 2**: All 4 agents run simultaneously with wave-based launches (30-60s stagger reduces resource contention)
- **All agents use Opus 4.6**: Maximum quality for all workflow phases (requirements analysis, research, edge cases, structural completeness, decision-making, review, execution)
- **Quality gates enforce standards**: No bypassing validation
- **Rollback on failure**: Safe to abort at any point
- **Planning journal provides traceability**: Full documentation of changes
- **Integration tested**: Works with sync-config.py and existing workflows

## References

See `skill-editor/references/` for:
- `session-management.sh`: Git safety checks, session creation/resume, cleanup commands
- `mode-detection.sh`: Three-tier complexity detection, mode selection prompt
- `phase-2-5-detection.sh`: Phase 2.5 complexity detection, strategy consultant, Quality Gate 2.5
- `experimental-tagging.sh`: Experimental mode YAML/comment tagging
- `anthropic-guidelines-summary.md`: Anthropic best practices
- `skill-structure-specification.md`: Skill format and validation
- `quality-gates.md`: Detailed quality gate checklists
- `config-management-integration.md`: Integration with CONFIG_MANAGEMENT.md

## Success Criteria

Skill-editor workflow succeeds when:

- [ ] User's original request is fulfilled
- [ ] All quality gates pass
- [ ] Changes are synced to `~/.claude/`
- [ ] Skill invokes without errors
- [ ] No regressions in existing skills
- [ ] Planning journal documents changes
- [ ] Changes committed to git
- [ ] User confirms satisfaction
