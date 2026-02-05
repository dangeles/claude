# Quality Gates

Definitions and criteria for all quality gates in the scientific-analysis-architect workflow.

## Gate Definitions

### Quality Gate 0: Initialization

**Phase**: 0
**Owner**: Orchestrator
**Type**: Automated

**Pass Criteria**:
- [ ] nbformat Python package available
- [ ] Session directory created successfully
- [ ] Session directory is writable (test file created and deleted)
- [ ] Output directory exists and is writable
- [ ] session-state.json initialized with valid schema

**Validation Code**:
```python
def validate_gate_0(session_state: dict) -> bool:
    """Validate initialization gate."""

    # Check nbformat
    try:
        import nbformat
    except ImportError:
        return False, "nbformat not installed"

    # Check session directory
    session_dir = session_state["config"]["session_directory"]
    if not os.path.exists(session_dir):
        return False, "Session directory does not exist"

    # Write test
    test_file = os.path.join(session_dir, ".write_test")
    try:
        with open(test_file, "w") as f:
            f.write("test")
        os.remove(test_file)
    except Exception as e:
        return False, f"Session directory not writable: {e}"

    # Check output directory
    output_dir = session_state["config"]["output_directory"]
    if not os.path.isdir(output_dir):
        return False, "Output directory does not exist"

    return True, "Gate 0 passed"
```

**Fail Action**: Abort workflow with clear error message

---

### Quality Gate 1: Research Structure

**Phase**: 1
**Owner**: research-architect
**Type**: Automated

**Pass Criteria**:
- [ ] research-structure.md file exists
- [ ] File is valid markdown
- [ ] Contains 3-7 chapters (inclusive)
- [ ] Each chapter has:
  - [ ] Title (non-empty string)
  - [ ] Goal (one of: atlas, hypothesis, mechanism)
  - [ ] At least 1 analysis listed
- [ ] Dependencies reference valid chapter numbers
- [ ] No circular dependencies

**Validation Code**:
```python
def validate_gate_1(session_dir: str) -> tuple:
    """Validate research structure gate."""

    structure_path = os.path.join(session_dir, "research-structure.md")

    # File exists
    if not os.path.exists(structure_path):
        return False, "research-structure.md not found"

    # Parse structure
    with open(structure_path) as f:
        content = f.read()

    chapters = parse_chapters(content)

    # Chapter count
    if not (3 <= len(chapters) <= 7):
        return False, f"Invalid chapter count: {len(chapters)} (expected 3-7)"

    # Validate each chapter
    valid_goals = {"atlas", "hypothesis", "mechanism"}
    for i, chapter in enumerate(chapters, 1):
        if not chapter.get("title"):
            return False, f"Chapter {i} missing title"
        if chapter.get("goal") not in valid_goals:
            return False, f"Chapter {i} invalid goal: {chapter.get('goal')}"
        if not chapter.get("analyses"):
            return False, f"Chapter {i} has no analyses"

    # Check dependencies
    chapter_nums = set(range(1, len(chapters) + 1))
    for chapter in chapters:
        for dep in chapter.get("dependencies", []):
            if dep not in chapter_nums:
                return False, f"Invalid dependency: Chapter {dep}"
            if dep >= chapter["number"]:
                return False, f"Circular/forward dependency: Chapter {chapter['number']} -> {dep}"

    return True, f"Gate 1 passed: {len(chapters)} chapters validated"
```

**Fail Action**: Return to Phase 1 with specific feedback

---

### Quality Gate 2: Chapter Plans

**Phase**: 2
**Owner**: analysis-planner
**Type**: Automated

**Pass Criteria**:
- [ ] All chapter{N}-notebook-plans.md files exist
- [ ] Each chapter plan has:
  - [ ] At least 1 notebook defined
  - [ ] Statistical approach for each notebook
  - [ ] No unresolved critical conflicts
- [ ] Data flow is consistent:
  - [ ] Outputs defined for each notebook
  - [ ] Inputs available from upstream notebooks or user data

**Validation Code**:
```python
def validate_gate_2(session_dir: str, num_chapters: int) -> tuple:
    """Validate chapter plans gate."""

    plans = []
    for i in range(1, num_chapters + 1):
        plan_path = os.path.join(session_dir, f"chapter{i}-notebook-plans.md")

        if not os.path.exists(plan_path):
            return False, f"chapter{i}-notebook-plans.md not found"

        with open(plan_path) as f:
            plan = parse_chapter_plan(f.read())

        # At least one notebook
        if not plan.get("notebooks"):
            return False, f"Chapter {i} has no notebooks"

        # Statistical approach required
        for nb in plan["notebooks"]:
            if not nb.get("statistical_approach"):
                return False, f"Chapter {i} notebook {nb['number']} missing statistical approach"

        # Check for unresolved conflicts
        if plan.get("unresolved_conflicts"):
            return False, f"Chapter {i} has unresolved conflicts: {plan['unresolved_conflicts']}"

        plans.append(plan)

    # Validate data flow
    available_outputs = set(["user_data"])
    for plan in plans:
        for nb in plan["notebooks"]:
            # Check inputs available
            for input_req in nb.get("data_requirements", {}).get("inputs", []):
                if input_req not in available_outputs:
                    return False, f"Input not available: {input_req}"

            # Add outputs to available
            for output in nb.get("data_requirements", {}).get("outputs", []):
                available_outputs.add(output)

    total_notebooks = sum(len(p["notebooks"]) for p in plans)
    return True, f"Gate 2 passed: {total_notebooks} notebooks across {num_chapters} chapters"
```

**Fail Action**: Return to Phase 2 with specific feedback

---

### Quality Gate 3: Structure Approval (USER)

**Phase**: 3
**Owner**: User
**Type**: Human judgment

**Presentation**:
```
Structure Review Complete

Summary:
- {N} chapters planned
- {M} notebooks total
- {K} issues identified (X critical, Y major, Z minor)

Critical Issues:
{list or "None"}

Major Issues:
{list or "None"}

Approve / Request changes / Reject? [A/c/r]
```

**Pass Criteria**:
- [ ] User explicitly approves (enters "A")
- [ ] OR user requests changes and changes are implemented

**Responses**:
- **A (Approve)**: Record approval, proceed to Phase 4
- **c (Changes)**: Gather feedback, re-run affected phases, return to gate
- **r (Reject)**: Ask which phase to return to (1 or 2)

**Fail Action**: Loop until approved or workflow aborted

---

### Quality Gate 4: Notebook Approval (USER)

**Phase**: 4
**Owner**: User
**Type**: Human judgment

**Presentation**:
```
Notebook Review Complete

Per-Chapter Summary:
- Chapter 1: {N} notebooks, {K} issues
- Chapter 2: {N} notebooks, {K} issues
...

Critical Issues:
{list or "None"}

Approve / Request changes / Reject? [A/c/r]
```

**Pass Criteria**:
- [ ] User explicitly approves (enters "A")
- [ ] OR all critical issues addressed after changes

**Responses**: Same as Gate 3

**Fail Action**: Loop until approved or workflow aborted

---

### Quality Gate 5: Notebook Validation

**Phase**: 5
**Owner**: notebook-generator
**Type**: Automated

**Pass Criteria**:
- [ ] All expected notebooks created
- [ ] Each notebook passes nbformat validation
- [ ] Each notebook has required metadata
- [ ] Backup copies exist in session directory

**Validation Code**:
```python
def validate_gate_5(session_state: dict) -> tuple:
    """Validate notebook generation gate."""

    import nbformat

    expected = session_state["outputs"].get("notebooks", [])
    if not expected:
        return False, "No notebooks generated"

    valid_count = 0
    for notebook_path in expected:
        # File exists
        if not os.path.exists(notebook_path):
            return False, f"Notebook not found: {notebook_path}"

        # nbformat validation
        try:
            with open(notebook_path) as f:
                nb = nbformat.read(f, as_version=4)
            nbformat.validate(nb)
        except Exception as e:
            return False, f"Invalid notebook {notebook_path}: {e}"

        # Metadata check
        if "scientific_analysis_architect" not in nb.metadata:
            return False, f"Missing provenance metadata: {notebook_path}"

        valid_count += 1

    # Check backups
    backup_dir = os.path.join(
        session_state["config"]["session_directory"],
        "notebooks"
    )
    if not os.path.exists(backup_dir):
        return False, "Backup directory not found"

    return True, f"Gate 5 passed: {valid_count} valid notebooks"
```

**Fail Action**:
- If partial failure: Offer to proceed with available or retry failed
- If total failure: Return to Phase 5 with error details

---

### Quality Gate 6: Statistical Review (USER)

**Phase**: 6
**Owner**: User
**Type**: Human judgment (Interview mode)

**Pass Criteria**:
- [ ] All concerns reviewed (accepted, rejected, or skipped)
- [ ] User confirmed correction application decision
- [ ] If corrections applied: notebooks re-validated

**Validation Code**:
```python
def validate_gate_6(session_state: dict) -> tuple:
    """Validate statistical review gate."""

    corrections = session_state.get("corrections", {})

    # All concerns have decisions
    total = (
        len(corrections.get("accepted", [])) +
        len(corrections.get("rejected", [])) +
        len(corrections.get("skipped", []))
    )

    expected = session_state.get("total_concerns", 0)
    if total < expected:
        return False, f"Not all concerns reviewed: {total}/{expected}"

    # If corrections applied, verify notebooks
    if corrections.get("applied"):
        for notebook_path in session_state["outputs"]["notebooks"]:
            try:
                with open(notebook_path) as f:
                    nb = nbformat.read(f, as_version=4)
                nbformat.validate(nb)
            except Exception as e:
                return False, f"Corrected notebook invalid: {notebook_path}: {e}"

    return True, "Gate 6 passed: Statistical review complete"
```

**Fail Action**: Resume interview from last concern

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
| 6 | Complete workflow | Resume interview |

## Minimum Thresholds

These thresholds cannot be bypassed:

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Minimum chapters | 3 | Ensures meaningful analysis structure |
| Maximum chapters | 7 | Prevents scope creep |
| Notebooks per chapter | >= 1 | Each chapter must produce output |
| Statistical approach | Required | Core purpose of skill |
| nbformat validation | Required | Ensures usable output |

## Quality Gate Flow Diagram

```
Start
  |
  v
[Gate 0] --FAIL--> Abort
  |
  PASS
  |
  v
Phase 1
  |
  v
[Gate 1] --FAIL--> Retry Phase 1
  |
  PASS
  |
  v
Phase 2
  |
  v
[Gate 2] --FAIL--> Retry Phase 2
  |
  PASS
  |
  v
Phase 3
  |
  v
[Gate 3: USER] --REJECT--> Return to Phase 1 or 2
  |           --CHANGES--> Implement changes, re-check
  |
  APPROVE
  |
  v
Phase 4
  |
  v
[Gate 4: USER] --REJECT--> Return to earlier phase
  |           --CHANGES--> Implement changes, re-check
  |
  APPROVE
  |
  v
Phase 5
  |
  v
[Gate 5] --PARTIAL--> Proceed with available / Retry failed
  |      --TOTAL FAIL--> Retry Phase 5
  |
  PASS
  |
  v
Phase 6
  |
  v
[Gate 6: USER] --INCOMPLETE--> Resume interview
  |
  COMPLETE
  |
  v
Workflow Complete
```

## Bypassing Gates

Gates 0, 1, 2, and 5 cannot be bypassed (automated validation).

Gates 3, 4, and 6 can be "bypassed" by user choice:
- Gate 3: User can approve with known issues
- Gate 4: User can approve with known issues
- Gate 6: User can skip remaining concerns

When a user bypasses a gate, a warning is logged:

```json
{
  "gate": 3,
  "bypassed": true,
  "reason": "User approved with 2 unresolved minor issues",
  "timestamp": "2026-02-04T14:55:00Z",
  "known_issues": ["Issue 1", "Issue 2"]
}
```
