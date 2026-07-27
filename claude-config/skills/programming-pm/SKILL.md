---
name: programming-pm
last_updated: 2026-06-09
description: Use when coordinating Python software development across multiple specialists with quality gates. NOT for ad-hoc bug fixes (use python-developer or copilot directly), non-Python work (use technical-pm), or single-developer feature work (use feature-dev plugin).

# v3.0 universal handoff metadata (see workflow-coordinator/references/frontmatter-metadata-standard.md)
handoff:
  accepts_from:
    - "*"
  provides_to:
    - python-developer
    - systems-architect
    - requirements-analyst
    - mathematician
    - statistician
    - copilot
    - skill-editor
  schema_version: "3.0"
  schema_type: universal
  # Legacy fields preserved for perspective-swarm v2.0 backward compatibility
  handoff_trigger: "--handoff {payload_path}"
  requires:
    - context.original_prompt
    - context.problem_type
  optional_consumes:
    - context.synthesis_summary
    - insights.convergent

categories:
  - implementation
  - architecture
  - project-management

input_requirements:
  - specification
  - requirements
  - architecture-decision

output_types:
  - implementation
  - code
  - tests
  - documentation
---

# Programming Project Manager

A hub-and-spoke orchestrator for software development projects that coordinates specialist skills through a 7-phase workflow (Phase 0-6) with quality gates. Central state lives with programming-pm; specialists do not communicate with each other directly.

## Orchestration and Delegation

Specialists own their domains: `systems-architect` for architecture, `mathematician` for algorithm design and complexity, `statistician` for statistical methods and Monte Carlo/MCMC validation, `python-developer` for implementation, `copilot` for review support, `notebook-writer` for notebooks, `requirements-analyst` for scoping.

Delegate specialist work when the subtask is substantial and independent — designing a multi-component architecture, implementing a component with its tests, designing a numerical algorithm, validating an MCMC implementation, reviewing a sizeable diff. Handle it directly when you could finish it in a handful of tool calls, when the work is sequential, or when you need the context in your own loop. See `../references/delegation-and-scope.md`.

**Orchestrator-owned work**: session setup, directory and state-file management, quality gate evaluation, user communication (summaries, approvals, status), workflow coordination and handoff creation, and running this skill's bash validation scripts.

If a specialist you need is unavailable, tell the user and apply the substitution in Specialist Availability rather than silently absorbing the work.

### Tool Selection

| Situation | Tool | Reason |
|-----------|------|--------|
| Specialist doing independent work | Task tool | Separate context, parallel execution |
| 2+ specialists working simultaneously | Task tool (multiple in one message) | Only way to parallelize |
| Loading domain knowledge for YOUR coordination decision | Skill tool | Shared context needed |

The bash blocks in this skill and its references are orchestration infrastructure you run yourself with the Bash tool — they are not specialist invocations.

## Dispatch Templates

Templates for each specialist invocation live in `references/dispatch-templates.md` (systems-architect Phase 3; python-developer Phase 4 and Phase 4 Step 0 pre-flight; copilot Phase 5; mathematician Phase 4; statistician Phase 4). Read that file when filling a dispatch payload.

**Content gate before dispatching Phase 1 requirements**: `requirements.problem_statement` is non-empty (> 50 characters) and `requirements.success_criteria` has at least 1 entry. If either fails, return to Phase 1 for clarification instead of dispatching.

**Paste verbatim**: fill templates with content copied verbatim from the source handoff files rather than summaries or paraphrases. If content is too large, pass the handoff file path and have the specialist read it directly.

## State Anchoring

Start every response with: "[Phase N/6 - {phase_name}] {brief status}"

Before starting a phase (Phase 1 onward), read `~/.claude/programming-pm-sessions/{workflow-id}/state.yaml` and confirm `current_phase` and `phases_completed` match expectations.

After a user interaction, answer the question, then re-anchor: "Returning to Phase N - {phase_name}. Next step: {action}."

While parallel agents are running, maintain a status board (Agent | Task | Running / Complete / Failed). When all agents have returned, proceed to quality gate evaluation.

## When to Use This Skill

Multi-component Python projects requiring architecture design and implementation; algorithm-heavy projects needing mathematician input for complexity analysis; statistical software requiring validation of Monte Carlo, MCMC, or bootstrap implementations; projects where implementation can be decomposed into independent parallel tasks; projects requiring formal quality gates (code review, testing, pre-mortem risk assessment).

## When NOT to Use This Skill

- **Simple scripts**: single-file Python scripts (<100 lines) — use copilot directly
- **Non-Python projects**: this skill is Python-first — use technical-pm
- **Bug fixes**: small changes to existing code — use python-developer or copilot
- **Research coordination**: literature reviews — use lit-pm
- **General feature work (non-bioinformatics)**: features not requiring the bioinformatics specialist team — use `feature-dev` for a lighter workflow (explore → clarify → architect → implement → review)
- **technical-pm instead**: coordinating research, writing, or analysis rather than code; tasks involving researcher, synthesizer, calculator rather than developers; flexible milestone tracking without quality gates; code is incidental rather than the primary deliverable

## Specialist Availability

The workflow needs `requirements-analyst` (Phase 1), `systems-architect` (Phase 3), and `copilot` (Phase 5) to run as designed. If one is missing, stop and give the user installation guidance.

Optional specialists degrade gracefully — inform the user, then:

- `edge-case-analyst` (Phase 2 pre-mortem): delegate a simplified pre-mortem to python-developer via Task tool. Alternatives: skip the pre-mortem, install edge-case-analyst, or have the user run it manually.
- `mathematician`: delegate algorithm design to python-developer; flag output "designed without specialist mathematician review."
- `statistician`: delegate statistical work to python-developer; flag "unvalidated -- no specialist statistician review."
- `notebook-writer`: delegate to python-developer with best-effort formatting; flag "created without notebook-writer specialized formatting."

Skill presence, when you need to check it, is `~/.claude/skills/{name}/SKILL.md` (some skills use lowercase `skill.md`).

## Tools

**Task** launches specialists for independent work; **Skill** loads domain knowledge into your own context for a coordination decision; **Read** covers existing codebase, patterns, and deliverables; **Write** covers deliverable documents, state files, and planning artifacts; **Bash** runs tests, linters, type checkers, git commands, and handoff validation scripts.

## Workflow State Persistence

Maintain workflow state in a YAML file for resume capability.

**State File**: `~/.claude/programming-pm-sessions/{workflow-id}/state.yaml`

```yaml
workflow:
  id: "prog-{project}-{date}"
  project_name: string
  created: ISO8601
  last_updated: ISO8601

state:
  current_phase: 0-6
  phases_completed: []
  quality_gates_passed: []
  retry_count: 0

session:
  session_dir: "~/.claude/programming-pm-sessions/{workflow-id}/"
  archival_guidelines_path: "{session_dir}/archival-guidelines-summary.md"
  guidelines_found: boolean
  guidelines_source: string  # Path to CLAUDE.md or "defaults"
  cleanup_on_complete: boolean  # Default true

team:
  composition: []
  active_tasks: []

artifacts:
  requirements: "/path/to/requirements.md"
  pre_mortem: "/path/to/pre-mortem.md"
  architecture: "/path/to/architecture.md"
  architecture_context: "/path/to/.architecture/context.md"  # Optional, generated in Phase 3
  implementation: []

exceptions:
  overrides: []
  accepted_risks: []
```

### State Recovery

On session resume: list sessions in `~/.claude/programming-pm-sessions/`, read the state file, check `last_updated` is within 7 days, display current phase and completed gates, then offer to continue from the current phase or restart.

## Workflow Phases

### Phase 0: Archival Guidelines Review

**Owner**: programming-pm (automatic) · **Checkpoint**: Never (always runs)

**Objective**: Initialize the workflow session directory and extract archival guidelines (preferring `.archive-metadata.yaml` over a CLAUDE.md fallback), with code-specific extraction focus.

**Entry**: Workflow invoked with a project goal.
**Exit**: `~/.claude/programming-pm-sessions/{workflow-id}/` created, archival-guidelines-summary.md written, `archival_context` block prepared for downstream handoffs.

**Key gate**: Session directory created, archival summary written. Session directory creation failure is fatal (ABORT — cannot proceed without session isolation).

**Full implementation detail** (session-dir creation, stale-session cleanup, primary/fallback archival extraction, output schema, downstream handoff block, failure handling, and the Phase 0→1 handoff validation script): see `references/phases/phase-0-session-setup.md`.

**Phase Transition**: Phase 0 complete -> Quality Gate 0 -> PROCEED to Phase 1: Requirements and Scoping

---

### Phase 1: Requirements and Scoping

**Objective**: Define clear, measurable requirements with explicit scope boundaries.
**Receives**: Session directory path and archival guidelines from Phase 0

**Steps**: Invoke `requirements-analyst` with the project goal and session context, review the resulting requirements document for completeness, and present it to the user for approval.

**Quality Gate 1: Requirements Approval**:
- Type: Human judgment (programming-pm review)
- Criteria:
  - [ ] Problem statement is specific (no vague terms like "better", "faster")
  - [ ] Success criteria are measurable (numbers, thresholds, or boolean conditions)
  - [ ] Scope boundaries (IN/OUT) explicitly defined
  - [ ] Dependencies identified
- Pass Condition: All criteria checked
- Fail Action: Return to requirements-analyst with feedback
- Override: User can accept partial requirements with documented gaps

**Handoff Validation** (Phase 1 → Phase 2): run `scripts/validate-handoff.py` on `${SESSION_DIR}/handoffs/phase1-requirements-handoff.yaml` with schema type `requirements_handoff`. Script and failure/override handling: `references/handoff-validation-gates.md`.

---

### Mode Selection (After Phase 1)

**Objective**: Select workflow execution mode based on project complexity.

**Trigger**: After Quality Gate 1 passes (requirements approved)

**Three execution modes**:
- **SIMPLE** (~1-2 hrs): Single component, no stats/math, <5 implementation tasks
- **STANDARD** (~4-6 hrs): Multi-component (2-5), optional stats/math, 5-15 tasks (default)
- **EXTENDED** (~8-12 hrs): >5 components OR both stats+math OR >15 tasks OR architectural complexity

**Tier triggers** (decision table):

| Tier | Select when |
|------|-------------|
| **EXTENDED** | Component count >5; OR requires BOTH stats AND math; OR task count >15; OR architectural-complexity keywords ("distributed system", "microservices", "event-driven", "real-time processing"); OR user requests "extended analysis"/"comprehensive review" |
| **SIMPLE** | Single component AND no stats/math; OR utility script with <5 tasks; OR single-component ETL pipeline |
| **STANDARD** (default) | Multiple components (2-5); single specialization (stats OR math, not both); 5-15 tasks; standard patterns ("web API", "CLI tool", "data analysis", "visualization") |

**Procedure**: run complexity detection via `detect_tier` from `references/mode-selection-criteria.md`, display the mode-selection prompt, confirm any user override (a risky SIMPLE override gets an explicit confirm), record the selection to `mode-selection.json` + `session-state.json`, export `PROGRAMMING_PM_MODE`, and branch behavior on the selected mode. Legacy sessions without `mode-selection.json` default to STANDARD.

**Full implementation detail** (all bash steps and the backwards-compatibility shim): see `references/phases/mode-selection.md`.

**Phase Transition**: Phase 1 complete -> Quality Gate 1 (user approval required) -> PROCEED to Phase 2: Pre-Mortem and Risk Assessment

---

### Phase 2: Pre-Mortem and Risk Assessment

**Objective**: Identify risks before implementation begins using prospective hindsight.

**Steps**: Invoke `edge-case-analyst` (or delegate a simplified pre-mortem to python-developer via Task tool), using the template in `assets/pre-mortem-template.md`. Document risks with likelihood, impact, and mitigation.

**Quality Gate 2: Pre-Mortem Completion**:
- Type: Automated (checklist validation)
- Criteria:
  - [ ] At least 3 risks identified
  - [ ] Each risk has likelihood rating (1-5) and impact rating (1-5)
  - [ ] Each risk has disposition: mitigate, accept, transfer, or avoid
  - [ ] Critical risks (score >= 15) have contingency plans
- Pass Condition: All risks have disposition
- Override: User can proceed with documented unmitigated risks

**Handoff Validation** (Phase 2 → Phase 3): run `scripts/validate-handoff.py` on `${SESSION_DIR}/handoffs/phase2-premortem-handoff.yaml` with schema type `premortem_handoff`. Script and failure/override handling: `references/handoff-validation-gates.md`.

**Phase Transition**: Phase 2 complete -> Quality Gate 2 -> PROCEED to Phase 3: Architecture Design

### Phase 3: Architecture Design

**Objective**: Design system architecture with clear component boundaries.

**Steps**:
1. **Dispatch systems-architect via Task tool** using the "Dispatch to systems-architect" template. Include the full requirements (Phase 1 handoff), all risks with architecture implications (Phase 2 handoff), and archival context (Phase 0).
2. **Receive and validate output**: architecture handoff YAML exists at the expected path and carries `producer: "systems-architect"`. Note whether the Architecture Context Document was generated or skipped.
3. **Review** against Quality Gate 3 criteria and **present to the user** for approval.

Scale the delegation to the project. Multi-component systems, non-obvious dependency graphs, scalability questions, and open technology choices go to systems-architect. For a SIMPLE-mode single-component project whose structure follows directly from the requirements, you can record the component interfaces, technology-choice rationale, and testing strategy yourself and note in session state that the architecture was authored inline.

#### Architecture Context Document

After architecture approval, systems-architect generates `.architecture/context.md`:

A persistent, version-controlled bird's-eye view of module structure, dependencies, and modification order for all implementation agents. **Content**: module interconnections (DAG), intended usage patterns, modification order for safe incremental changes, streaming/incremental strategies. **Lifecycle**: created in Phase 3, read by developers before implementation, updated when architectural changes occur (Phase 5 drift check).

See `systems-architect/references/architecture-context-template.md` for template details.

**Quality Gate 3: Architecture Approval**:
- Type: Human judgment (programming-pm + user review)
- Criteria:
  - [ ] All components identified with responsibilities
  - [ ] Data flow documented (inputs, outputs, transformations)
  - [ ] Technology choices justified (libraries, frameworks)
  - [ ] Component interfaces defined
  - [ ] Testing strategy outlined
  - [ ] Architecture Context Document generated (`.architecture/context.md` exists)
- Override: User can approve partial architecture for proof-of-concept

**Handoff Validation** (Phase 3 → Phase 4): run `scripts/validate-handoff.py` on `${SESSION_DIR}/handoffs/phase3-architecture-handoff.yaml` with schema type `architecture_handoff`. Incomplete architecture does not proceed to implementation. Script and failure/override handling: `references/handoff-validation-gates.md`.

**Phase Transition**: Phase 3 complete -> Quality Gate 3 (user approval required) -> PROCEED to Phase 4: Implementation

### Phase 4: Implementation

**Objective**: Implement architecture with specialist agents, dispatched via Task tool.

**Mode-based execution**:
- **SIMPLE**: Sequential execution (one specialist at a time)
- **STANDARD/EXTENDED**: Wave-based parallel execution

**Entry**: Architecture approved (Phase 3, Quality Gate 3) and `.architecture/context.md` exists.
**Exit**: All non-escalated tasks have deliverables; per-task code/math/stats handoffs validate against schema.

**Step overview** (full mechanics in the reference below):

| Step | Purpose |
|------|---------|
| **0. Implementability check** | Confirm the architecture handoff is implementable: interfaces typed, dependencies acyclic, tech choices compatible, order consistent. Read the handoff yourself when it is small; dispatch python-developer via Task tool when it is large or the gaps look substantive. FAIL → escalate to Phase 3 (max 2 cycles). SIMPLE mode with <=1 component skips this step. |
| **1. Task decomposition** | Parse architecture handoff into tasks; assign each to a specialist. Primary signal = `specialist_flags` (v1.4+ handoffs); fallback = keyword heuristic. Algorithm→mathematician, stats→statistician, notebook→notebook-writer, else python-developer. State each task's scope in its dispatch description. |
| **2. Wave-based execution** | SIMPLE = sequential. STANDARD/EXTENDED = three waves: Wave 1 launches math/stats (feed other tasks), Wave 2 launches independent implementation tasks, Wave 3 launches dependency-blocked tasks once their inputs land (bounded, then escalates). |
| **3. Completion tracking** | Record each specialist's deliverable as it returns; escalate tasks whose dependencies never resolve. |
| **4. Gate 4a completion check** | PASS at 100%; CONDITIONAL PASS at 75%+ only if all critical (math/stats) specialists completed; FAIL below 75%. |

**Key gate — Quality Gate 4b (Implementation Validation)**:
- Type: Automated (output validation)
- Criteria:
  - [ ] All specialist outputs exist and are >100 words
  - [ ] No critical blocking issues flagged
  - [ ] Handoffs validate against schema (validate-handoff.py)
  - [ ] Acceptance criteria met per task
- Pass Condition: All criteria checked OR 75%+ with critical specialists complete
- Fail Action: Retry incomplete tasks or escalate to user

**Full implementation detail** (Step 0 dispatch and FAIL escalation, Step 1 yq task-decomposition with specialist-flag normalization and keyword fallback, Step 2 wave-launch protocol and Wave-3 bounded retry, Step 4a completion decision table, and the Phase 4→5 handoff-validation script): see `references/phases/phase-4-implementation.md`.

**Phase Transition**: Phase 4 complete -> Quality Gate 4 -> PROCEED to Phase 5: Code Review and Testing

### Phase 5: Code Review and Testing

**Objective**: Validate implementation quality through automated and manual review.

**Steps**:
1. **Run automated checks**: linting (`ruff check .`), type checking (`mypy --strict src/` or `pyright src/`), tests (`pytest --cov`).
2. **Review the implementation.** Dispatch `copilot` via Task tool using the "Dispatch to copilot" template when a context-isolated adversarial pass adds something — several components changed, subtle numerics, security-sensitive input handling, or automated checks that passed but left you unsure. For small diffs, read the changed files and review them yourself. Record which route you took in session state.
3. **If copilot ran**: read `{SESSION_DIR}/deliverables/copilot-review.md`. If the file is absent, use copilot's Task tool return value as the review. A review without a `VERDICT:` line is incomplete — ask for the verdict.
4. **Act on findings by severity**:
   - CRITICAL: return to python-developer with the copilot feedback for fixes
   - MAJOR: present to the user, who decides between fixing and documented tech debt
   - MINOR/GOOD: proceed
   Re-review after CRITICAL fixes when the fix itself was substantial. If CRITICAL issues persist, present them to the user rather than looping further.
5. **(Optional) Supplemental review via pr-review-toolkit**: invoke `/pr-review-toolkit:review-pr` for specialized analysis copilot does not cover (silent failures, test coverage, type design, comment accuracy). Advisory — findings do not block Phase 6 but are presented to the user.
6. If deliverables include notebooks: invoke `notebook-writer` for reproducibility review.
7. Check for architecture drift and update the context document if needed.

#### Architecture Drift Check

If `.architecture/context.md` exists, check whether implementation introduced structural changes that require a context update:

**Drift detection heuristics** (narrow scope to reduce false positives):
- New files in `src/` or `modules/` directories
- Deleted module directories
- Changes to `__init__.py` files (interface changes)
- Developer reported discrepancy via `architecture_context.discrepancy_noted: true` in code handoff

**On drift**: invoke `systems-architect` for a targeted update of the affected sections of `.architecture/context.md` (not full regeneration), and commit the context update with the implementation changes. **On no drift**: proceed to Quality Gate 4.

This is a lightweight check. Fundamental architectural changes (a new module changing dependency-graph topology) are logged as "architectural drift requiring future Phase 3 review" rather than triggering heavyweight updates within Phase 5.

**Quality Gate 4: Code Review Approval**:
- Type: Human judgment (programming-pm + copilot review)
- Automated checks (all pass):
  - [ ] `ruff check .` returns 0 errors
  - [ ] Type checking passes: `mypy --strict src/` or `pyright src/` returns 0 errors (warnings acceptable)
  - [ ] Test coverage >= 80% for new code
- Review completed:
  - [ ] Review performed and recorded (copilot review document, or documented inline review)
  - [ ] No CRITICAL issues remain unaddressed
  - [ ] MAJOR issues addressed or documented as accepted tech debt
  - Override: programming-pm can approve with a "tech debt" tag if the deadline is critical; if copilot was invoked and unreachable, document "copilot unavailable" in session state and proceed with a documented inline review only, reporting this to the user
- Human review:
  - [ ] Code matches requirements specification
  - [ ] Edge cases from pre-mortem are handled
  - [ ] Documentation present (docstrings, type hints)
  - [ ] No obvious security issues
- Fail Action: Return to developer with specific feedback

**Quality Gate 5: Test Pass**:
- Type: Automated (test execution)
- Criteria:
  - [ ] All unit tests pass
  - [ ] All integration tests pass (if applicable)
  - [ ] Coverage >= 80% for new code
  - [ ] No regressions in existing tests
- Override: User can merge with failing tests for emergency (creates P0 issue)

**Handoff Validation** (Phase 5 → Phase 6): run `scripts/validate-handoff.py` on `${SESSION_DIR}/handoffs/phase5-review-handoff.yaml` with schema type `review_handoff`. An incomplete review handoff does not proceed to merge. Script and failure/override handling: `references/handoff-validation-gates.md`.

**Phase Transition**: Phase 5 complete -> Quality Gate 5 -> PROCEED to Phase 6: Version Control Integration

### Phase 6: Version Control Integration

**Objective**: Integrate changes with sync-config.py and version control.

**Entry**: Code review approved (Phase 5, Quality Gate 5); review handoff validates.
**Exit**: Changes committed (on a feature branch if starting from main/master), synced via sync-config.py, planning journal entry created, session directory cleaned up or preserved.

**Optional shortcuts** (use when git strategy is simple):
- **Commit Commands plugin**: `/commit-push-pr` handles commit + push + PR in one step, replacing manual Steps 2-3.
- **Git Strategy Advisory**: optionally invoke `git-strategy-advisor` (Task tool, `mode: post-work`) for scope-adaptive recommendations. Advisory only — Step 3's feature-branch logic takes precedence unconditionally; present the advisor's `summary` as an informational note. Skip silently on confidence "none"; caveat on "low".

**Step overview** (full mechanics in the reference below):

| Step | Purpose |
|------|---------|
| **1. Pre-merge validation** | Check sync-config.py availability; re-validate every handoff; `push --dry-run` to detect conflicts before committing. |
| **2. Quality Gate 6 validation** | Run `validate-gate.sh 6`; on FAIL, offer logged override. |
| **3. Commit and sync** | Create a feature branch if on main/master; stage specific files only (NEVER `git add .`); conventional commit; `sync-config.py push`. |
| **4. Planning journal entry** | `sync-config.py plan --title` documenting objective, specialists, files changed, testing, outcome. |
| **5. Session cleanup** | Mark session completed; delete on clean success, preserve if any override or validation error occurred. |

**Key gate — Quality Gate 6**:
- Type: Automated (VCS checks)
- Criteria:
  - [ ] All previous gates passed (or overrides documented)
  - [ ] No merge conflicts (verified by dry-run)
  - [ ] Review approved (from phase5-review-handoff.yaml)
  - [ ] Deliverable location documented
  - [ ] Files staged (if in git repo)
- Override: Repository admin can force merge (logged for audit)

**Post-Merge Verification**: After sync, prompt the user to verify the deliverable meets expectations. If issues are found, create a follow-up task (not a rollback unless critical).

**Full implementation detail** (the optional-shortcut/advisory notes in full plus every step's bash — pre-merge handoff revalidation loop, dry-run conflict handling, Gate 6 validation/override, branch + specific-file staging + conventional commit + sync, planning journal entry, and session-cleanup decision logic): see `references/phases/phase-6-version-control.md`.

## Gate Override Protocol

When a quality gate fails:

1. **Display failure details** with severity: CRITICAL (cannot override — security, runtime errors), HIGH (override requires explicit user approval), MEDIUM (override allowed with documentation), LOW (override allowed).
2. **Offer options**: [Fix] address issues and re-run the gate; [Override] proceed with documented risk acceptance; [Escalate] consult a specialist for a second opinion.
3. **If Override selected**: log the decision with timestamp, user, and rationale; mark the deliverable "GATE_OVERRIDE: {gate_name}"; continue the pipeline but flag it in the final PR description.

**Override cannot skip**: test failures indicating runtime errors, security vulnerabilities (P0), architecture compatibility failures.

## Exception Handling

**Specialist stalls or returns incomplete work**: diagnose from the specialist's output and progress file, then present options to the user — extend, narrow the scope, substitute a specialist (e.g. python-developer for mathematician), or escalate for guidance. Apply the choice and log it to `exceptions-log.md`. Per-phase and per-specialist duration expectations and the intervention schema: `references/timeout-config.md`.

**Circuit breaker**: after 3 consecutive failures of the same type, stop retrying, present the failure summary, and require an explicit user decision — retry with changes, skip the component, or abort the workflow.

## Role Conflict Resolution

Authority by domain: architecture decisions (Phase 3) rest with systems-architect; algorithm design with mathematician; statistical methods with statistician; implementation decisions (Phase 4) with python-developer within the architecture constraints.

When specialists return contradictory recommendations, classify the conflict. Minor conflicts (implementation detail) are python-developer's call. Major conflicts (an architecture change is required) go to the user with both positions and their rationales, the available options including any hybrid, and your own recommendation. Document the resolution in the architecture spec.

## Team Composition

| Skill | Role | Phase |
|-------|------|-------|
| programming-pm | Orchestrator | All |
| requirements-analyst | Requirements scoping | 1 |
| systems-architect | Architecture design | 3 |
| python-developer | Implementation | 4-5 |
| copilot | Code review support | 5 |

The table above is the **default team** — always included. For optional specialists (mathematician, statistician, notebook-writer) with their inclusion criteria, the team-selection decision tree, the RACI matrix, team-size guidelines, and `--include`/`--exclude`/`--minimal` override flags, see `references/team-composition.md`.

## Handoff Format

All handoffs between specialists use a standardized schema. See `references/handoff-schema.md`.

**Base handoff fields**:
```yaml
handoff:
  version: "1.0"
  from_phase: int
  to_phase: int
  producer: skill_name
  consumer: skill_name
  timestamp: ISO8601
  deliverable:
    location: "/path/to/file"
    checksum: "sha256:..."
  context:
    focus_areas: []
    known_gaps: []
  quality:
    status: "complete" | "partial"
    confidence: "high" | "medium" | "low"
```

## Supporting Resources

- `assets/pre-mortem-template.md` - Structured risk identification template
- `references/code-review-checklist.md` - Quality gate criteria for code review
- `references/git-workflow.md` - Branching strategy, commit format, rollback procedures
- `references/team-composition.md` - RACI matrix, specialist selection criteria
- `references/handoff-schema.md` - Interface contracts between specialists
- `references/handoff-validation-gates.md` - Handoff validation scripts and override handling per phase transition
- `references/timeout-config.md` - Per-phase and per-specialist duration expectations, intervention protocol
- `references/phases/phase-0-session-setup.md` - Phase 0 full detail (session setup, archival extraction)
- `references/phases/mode-selection.md` - Mode Selection full detail (complexity detection + branching)
- `references/phases/phase-4-implementation.md` - Phase 4 full detail (implementability check, decomposition, waves, gates)
- `references/phases/phase-6-version-control.md` - Phase 6 full detail (pre-merge, gate 6, commit/sync, journal, cleanup)
- `git-strategy-advisor` - Phase 6 git strategy consultation (optional, advisory)

## Integration with Existing Skills

This skill invokes but does not modify `requirements-analyst` (Phase 1), `systems-architect` (Phase 3), `copilot` (Phase 5), and `edge-case-analyst` (Phase 2 pre-mortem, optional).

It coordinates `python-developer` (Phase 4-5 implementation and fixes), `mathematician` (Phase 4 algorithm design), `statistician` (Phase 4 statistical validation), and `notebook-writer` (Phase 4 notebook creation, Phase 5 notebook review).
