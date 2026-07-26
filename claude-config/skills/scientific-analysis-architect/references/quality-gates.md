# Quality Gates

Criteria for the quality gates in the scientific-analysis-architect workflow, plus the runnable validators for the two document-structure gates.

## Contents

- Gate 0: Initialization
- Gate 1: Research Structure
- Gate 2: Chapter Plans
- Gate 3: Structure Approval (user)
- Gate 4: Notebook Approval (user)
- Gate 5: Analysis Document Validation (validator)
- Gate 6: Statistical Review (user)
- Gate 7: Audience Document Validation (validator)
- Pass/Fail Actions Summary
- Minimum Thresholds
- Bypassing Gates

---

## Gate 0: Initialization

Phase 0, orchestrator, automated.

- Session directory created successfully and writable (test file created and deleted)
- Output directory exists and is writable
- `session-state.json` initialized with a valid schema

On failure: stop with a clear error message.

---

## Gate 1: Research Structure

Phase 1, research-architect, automated.

- `research-structure.md` exists and is valid markdown
- Contains 3-7 chapters inclusive
- Each chapter has a non-empty title, a goal from {atlas, hypothesis, mechanism}, and at least 1 analysis
- Dependencies reference valid chapter numbers, all earlier than the referencing chapter (no cycles, no forward references)

On failure: return to Phase 1 with the specific problem.

---

## Gate 2: Chapter Plans

Phase 2, analysis-planner, automated.

- Every `chapter{N}-notebook-plans.md` exists
- Each plan has at least 1 analysis, and every analysis has a statistical approach
- No unresolved critical conflicts recorded in the plan
- Data flow is consistent: every analysis input is either user data or an output produced upstream, and outputs are declared

On failure: return to Phase 2 with the specific problem.

---

## Gate 3: Structure Approval (USER)

Phase 3, human judgment.

```
Structure Review Complete

Summary:
- {N} chapters planned
- {M} analyses total
- {K} issues identified (X critical, Y major, Z minor)

Critical Issues:
{list or "None"}

Major Issues:
{list or "None"}

Approve / Request changes / Reject? [A/c/r]
```

Passes when the user approves ("A"), or when they request changes and those changes are implemented.

- **A (Approve)**: record approval, proceed to Phase 4
- **c (Changes)**: gather feedback, re-run affected phases, return to this gate
- **r (Reject)**: ask which phase to return to (1 or 2)

---

## Gate 4: Notebook Approval (USER)

Phase 4, human judgment.

```
Plan Review Complete

Per-Chapter Summary:
- Chapter 1: {N} analyses, {K} issues
- Chapter 2: {N} analyses, {K} issues
...

Critical Issues:
{list or "None"}

Approve / Request changes / Reject? [A/c/r]
```

Passes when the user approves, or when all critical issues are addressed after changes. Responses are handled as in Gate 3.

---

## Gate 5: Analysis Document Validation

Phase 5, orchestrator plus notebook-generator, automated.

- All expected analysis documents created
- Each has the required sections: Goal, Statistical Approach, Analysis Steps, Expected Outputs
- Each has at least one fenced code block, with balanced fences
- Master strategy overview exists with its required sections
- Backup copies exist in `{session_dir}/analyses/`

```python
import os
import re

REQUIRED_SECTIONS = {
    "goal": [r'^##\s+(Goal|Objective|Goals)\b'],
    "statistical_approach": [r'^##\s+(Statistical Approach|Statistical Method|Methods)\b'],
    "analysis_steps": [r'^##\s+(Analysis Steps|Steps|Workflow Steps)\b'],
    "expected_outputs": [r'^##\s+(Expected Outputs|Outputs|Results)\b']
}

def validate_analysis_document(path: str) -> tuple:
    """Validate markdown analysis document structure."""
    if not os.path.exists(path):
        return False, f"Document not found: {path}"

    with open(path) as f:
        content = f.read()

    if len(content.strip()) == 0:
        return False, f"Empty file: {path}"

    # Check required sections
    missing = []
    for section_name, patterns in REQUIRED_SECTIONS.items():
        found = any(
            re.search(p, content, re.MULTILINE | re.IGNORECASE)
            for p in patterns
        )
        if not found:
            missing.append(section_name)

    if missing:
        return False, f"Missing sections in {path}: {missing}"

    # Check for at least one fenced code block
    if '```' not in content:
        return False, f"No code blocks found in {path}"

    # Check balanced fences
    fence_count = content.count('```')
    if fence_count % 2 != 0:
        return False, f"Unbalanced code fences in {path} ({fence_count} backtick markers)"

    return True, f"Valid: {path}"

def validate_gate_5(session_state: dict) -> tuple:
    """Validate markdown analysis document generation (Gate 5)."""
    expected = session_state["outputs"].get("analyses", [])
    if not expected:
        return False, "No analysis documents generated"

    for doc_path in expected:
        valid, msg = validate_analysis_document(doc_path)
        if not valid:
            return False, msg

    # Check master overview
    overview_path = os.path.join(
        session_state["config"]["output_directory"],
        "analysis-strategy-overview.md"
    )
    if not os.path.exists(overview_path):
        return False, "Master strategy overview not found"

    # Check overview required sections
    OVERVIEW_SECTIONS = [
        r'^##\s+Project Objective',
        r'^##\s+Strategy at a Glance',
        r'^##\s+Chapter Summaries',
        r'^##\s+Data Flow',
        r'^##\s+Execution Order'
    ]
    with open(overview_path) as f:
        overview_content = f.read()
    for pattern in OVERVIEW_SECTIONS:
        if not re.search(pattern, overview_content, re.MULTILINE):
            return False, f"Master overview missing section: {pattern}"

    # Check backup directory
    backup_dir = os.path.join(
        session_state["config"]["session_directory"],
        "analyses"
    )
    if not os.path.exists(backup_dir):
        return False, "Backup directory not found"

    return True, f"Gate 5 passed: {len(expected)} valid analysis documents + master overview"
```

On partial failure: offer to proceed with what exists or retry the failed chapters. On total failure: return to Phase 5 with the error details.

---

## Gate 6: Statistical Review (USER)

Phase 6, human judgment via interview mode.

- Every concern has a decision recorded (accepted, rejected, or skipped)
- User confirmed the correction application decision
- If corrections were applied: analysis documents re-validated with `validate_analysis_document`, and the master overview refreshed and re-validated

On failure: resume the interview from the last concern.

---

## Gate 7: Audience Document Validation

Phase 7, orchestrator, automated.

- `{output_dir}/researcher-plan.md` exists and is non-empty (> 500 bytes)
- `{output_dir}/.research-architecture/architect-handoff.md` exists and is non-empty
- `{output_dir}/.research-architecture/engineering-translation.md` exists and is non-empty
- Researcher plan has required sections (case-insensitive): Research Overview, Research Questions, Expected Outcomes, Decision Points
- Architect handoff has required sections: Design Rationale, Current State, Open Questions, Continuation Guidance
- Engineering translation has required sections: System Overview, Pipeline Architecture, Data Specifications, Processing Stages, Resource Requirements, Dependencies
- No fenced code blocks (```) in `researcher-plan.md`
- Backup copies exist in `{session_dir}/audience-documents/`
- All three documents include provenance metadata (HTML comments with session ID)

```python
import re
import os

RESEARCHER_SECTIONS = {
    "research_overview": [r'^##\s+Research Overview'],
    "research_questions": [r'^##\s+Research Questions'],
    "expected_outcomes": [r'^##\s+Expected Outcomes'],
    "decision_points": [r'^##\s+Decision Points']
}

ARCHITECT_SECTIONS = {
    "design_rationale": [r'^##\s+Design Rationale'],
    "current_state": [r'^##\s+Current State'],
    "open_questions": [r'^##\s+Open Questions'],
    "continuation": [r'^##\s+Continuation Guidance']
}

ENGINEERING_SECTIONS = {
    "system_overview": [r'^##\s+System Overview'],
    "pipeline_architecture": [r'^##\s+Pipeline Architecture'],
    "data_specifications": [r'^##\s+Data Specifications'],
    "processing_stages": [r'^##\s+Processing Stages'],
    "resource_requirements": [r'^##\s+Resource Requirements'],
    "dependencies": [r'^##\s+Dependencies']
}

def validate_gate_7(output_dir: str, session_dir: str) -> tuple:
    """Validate audience document gate."""

    researcher_plan = os.path.join(output_dir, "researcher-plan.md")
    architect_handoff = os.path.join(output_dir, ".research-architecture", "architect-handoff.md")
    engineering_translation = os.path.join(output_dir, ".research-architecture", "engineering-translation.md")

    # Check existence and minimum size
    for path, name in [(researcher_plan, "researcher plan"),
                       (architect_handoff, "architect handoff"),
                       (engineering_translation, "engineering translation")]:
        if not os.path.exists(path):
            return False, f"{name} does not exist"
        if os.path.getsize(path) < 500:
            return False, f"{name} is too small (< 500 bytes)"

    # Read contents
    with open(researcher_plan) as f:
        researcher_content = f.read()
    with open(architect_handoff) as f:
        architect_content = f.read()
    with open(engineering_translation) as f:
        engineering_content = f.read()

    # Check sections (case-insensitive)
    for content, sections, name in [
        (researcher_content, RESEARCHER_SECTIONS, "researcher plan"),
        (architect_content, ARCHITECT_SECTIONS, "architect handoff"),
        (engineering_content, ENGINEERING_SECTIONS, "engineering translation")
    ]:
        for section_key, patterns in sections.items():
            found = False
            for pattern in patterns:
                if re.search(pattern, content, re.MULTILINE | re.IGNORECASE):
                    found = True
                    break
            if not found:
                return False, f"{name} missing section: {section_key}"

    # Check no code blocks in researcher plan
    if '```' in researcher_content:
        return False, "researcher plan contains code blocks (must be prose-only)"

    # Check backups
    backup_dir = os.path.join(session_dir, "audience-documents")
    for filename in ["researcher-plan.md", "architect-handoff.md", "engineering-translation.md"]:
        backup_path = os.path.join(backup_dir, filename)
        if not os.path.exists(backup_path):
            return False, f"backup copy missing: {filename}"

    # Check provenance metadata
    for content, name in [(researcher_content, "researcher plan"),
                          (architect_content, "architect handoff"),
                          (engineering_content, "engineering translation")]:
        if "Session:" not in content:
            return False, f"{name} missing provenance metadata"

    return True, "Gate 7 passed: All audience documents validated"
```

On failure: retry failed documents (up to 2 attempts each); if backups fail but documents exist, retry the backup only and warn if it still fails; on total failure, proceed without audience documents, warn the user, and mark Phase 7 "degraded".

---

## Pass/Fail Actions Summary

| Gate | On Pass | On Fail |
|------|---------|---------|
| 0 | Proceed to Phase 1 | Abort workflow |
| 1 | Proceed to Phase 2 | Retry Phase 1 |
| 2 | Proceed to Phase 3 | Retry Phase 2 |
| 3 | Proceed to Phase 4 | Handle user choice |
| 4 | Proceed to Phase 5 | Handle user choice |
| 5 | Proceed to Phase 6 | Partial proceed or retry |
| 6 | Proceed to Phase 7 | Resume interview |
| 7 | Complete workflow | Retry documents or proceed without |

## Minimum Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Minimum chapters | 3 | Ensures meaningful analysis structure |
| Maximum chapters | 7 | Prevents scope creep |
| Analysis documents per chapter | >= 1 | Each chapter must produce output |
| Statistical approach | Required | Core purpose of skill |
| Markdown structure validation | Required | Ensures usable output |

## Bypassing Gates

Gates 0, 1, 2, 5, and 7 are automated validation and are not bypassed. Gates 3, 4, and 6 are the user's: they can approve with known issues (3, 4) or skip remaining concerns (6). Log a bypass when it happens:

```json
{
  "gate": 3,
  "bypassed": true,
  "reason": "User approved with 2 unresolved minor issues",
  "timestamp": "2026-02-04T14:55:00Z",
  "known_issues": ["Issue 1", "Issue 2"]
}
```
