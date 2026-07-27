# Dispatch Templates

Templates for invoking each specialist via the Task tool. Fill all `{placeholder}` values from session state and handoff files before dispatching.

## Contents

- Content gate for Phase 1 handoff
- Paste verbatim discipline
- Template: Dispatch to systems-architect (Phase 3)
- Template: Dispatch to python-developer (Phase 4)
- Template: Dispatch to python-developer (Phase 4 Step 0 Pre-flight)
- Template: Dispatch to copilot (Phase 5)
- Template: Dispatch to mathematician (Phase 4)
- Template: Dispatch to statistician (Phase 4)

## Content gate for Phase 1 handoff

Before pasting requirements into any template, check that `requirements.problem_statement`
is non-empty (> 50 characters) and `requirements.success_criteria` has at least 1 entry.
If either check fails, return to Phase 1 for requirements clarification instead of dispatching.

## Paste verbatim discipline

Fill templates with content copied verbatim from the source handoff files rather than summaries or paraphrases. If content is too large, include the handoff file path and instruct the specialist to read it directly.

---

## Template: Dispatch to systems-architect (Phase 3)

```
Use the systems-architect skill to design the system architecture.

## Requirements Summary
{paste verbatim from Phase 1 handoff: problem_statement, success_criteria, scope}

## Risk Assessment Summary
{paste verbatim from Phase 2 handoff: critical risks, architecture implications}

## Archival Context
{paste archival_context block from Phase 0}

## Expected Deliverables
1. Architecture document with components, data flow, technology choices
2. Architecture Context Document (.architecture/context.md) if applicable
3. Architecture handoff YAML at: {SESSION_DIR}/handoffs/phase3-architecture-handoff.yaml
4. Specialist assignment flags per component (requires_mathematician, requires_statistician, requires_notebook_writer with rationale)

## Required Handoff Field
Your architecture handoff YAML includes: `producer: "systems-architect"`

## Constraints
- Output must follow architecture_handoff schema (see handoff-schema.md)
- All components must have defined interfaces (inputs/outputs)
- Implementation order must be specified
- Each component includes specialist_flags with explicit boolean values and a rationale field
- SIMPLE mode minimum: at least 1 component with interfaces, 1 technology choice rationale, 1 testing strategy
```

---

## Template: Dispatch to python-developer (Phase 4)

```
Use the python-developer skill to implement the following component.

## Task Assignment
- Task ID: {task_id}
- Component: {component_name}
- Description: {component_description}
- Scope: {single function | single module | multi-file component}

## Architecture Specification
{paste verbatim relevant component spec from Phase 3 handoff}

## Architecture Context
- Context document path: {path to .architecture/context.md, or "none" if not generated}
- Component tier: {tier from context doc, or "unknown" if context doc not available}

## Dependencies
- Upstream: {list of completed tasks this depends on}
- Downstream: {list of tasks waiting on this}

## Acceptance Criteria
{list from task decomposition}

## Archival Context
{paste archival_context block}

## Deliverables
- Implementation files in: {output directory}
- Test files in: {test directory}
- Code handoff YAML at: {SESSION_DIR}/handoffs/phase4-code-handoff-{task_id}.yaml
```

---

## Template: Dispatch to python-developer (Phase 4 Step 0 Pre-flight)

```
Use the python-developer skill to validate the architecture handoff for implementability.
This is PRE-FLIGHT VALIDATION only — do not implement anything.

## Validation Task
Review the architecture handoff and produce a preflight_validation report.

## Architecture Handoff
{If <= 200 lines: paste verbatim. If > 200 lines: "Read architecture from: {SESSION_DIR}/handoffs/phase3-architecture-handoff.yaml"}

## Requirements Summary
{paste verbatim from Phase 1 handoff: problem_statement, success_criteria}

## Validation Checklist
1. All component interfaces fully specified (input types, output types present)
2. Component dependency chains are traceable without cycles
3. Technology choices are recognizable
4. Implementation order is consistent with dependencies

## Expected Output
Write validation report to: {SESSION_DIR}/deliverables/phase4-preflight.yaml

status: PASS | FAIL | PASS_WITH_WARNINGS
For each gap: component name, issue, severity (BLOCKING or WARNING)
Recommendations for resolving BLOCKING issues.
```

---

## Template: Dispatch to copilot (Phase 5)

```
Use the copilot skill to perform adversarial code review on the following implementation.

Note: Copilot performs static code review (Read-based analysis only). Automated checks
(linting, testing, coverage) were already run in Phase 5 Step 1, before this invocation.

## Review Scope
- Python files to review: {list all Python files changed in Phase 4, from code handoffs}
- Non-Python files (YAML, Docker, shell, etc.) are listed for awareness only; copilot review focuses on Python. Non-Python validation relies on automated checks.
- Total Python files: {count}

## Review Focus Areas
1. Correctness: Logic errors, off-by-one, wrong operators
2. Edge cases: Empty input, boundary values, missing data
3. Performance: Vectorization, memory efficiency, algorithmic complexity
4. Security: Input validation, injection, unsafe operations

## Requirements Context
{paste verbatim problem_statement and success_criteria from Phase 1 handoff}

## Pre-Mortem Risks to Verify
{paste verbatim critical and high risks from Phase 2 — verify these are handled in code}

## Architecture Context
{paste verbatim component descriptions from Phase 3 — verify implementation matches design}

## Expected Output
Write a review document to: {SESSION_DIR}/deliverables/copilot-review.md

The review includes:
- CRITICAL issues (must fix before merge) — with file:line, description, impact, fix
- MAJOR issues (should fix) — with file:line and suggestion
- MINOR issues (nice to have)
- GOOD practices observed
- Final line: `VERDICT: APPROVED` or `VERDICT: NEEDS REVISION`
```

---

## Template: Dispatch to mathematician (Phase 4)

```
Use the mathematician skill to design the algorithm for the following component.

## Problem Description
{paste verbatim algorithm requirements from architecture handoff}

## Performance Constraints
{paste verbatim performance requirements from architecture}

## Expected Deliverables
1. Algorithm specification with pseudocode
2. Complexity analysis (time and space)
3. Numerical stability assessment
4. Implementation guidance for python-developer
5. Math handoff YAML at: {SESSION_DIR}/handoffs/phase4-math-handoff-{task_id}.yaml
```

---

## Template: Dispatch to statistician (Phase 4)

```
Use the statistician skill to design the statistical approach for the following component.

## Statistical Problem
{paste verbatim statistical requirements from architecture handoff}

## Data Characteristics
{paste verbatim data description from requirements}

## Expected Deliverables
1. Method selection with rationale
2. Validation criteria and diagnostic checks
3. Implementation guidance for python-developer
4. Stats handoff YAML at: {SESSION_DIR}/handoffs/phase4-stats-handoff-{task_id}.yaml
```
