# Handoff Schema Specification

Interface contracts for specialist communication in programming-pm orchestrated projects.

---

## Overview

All handoffs between specialists follow standardized schemas to ensure:
- Consistent communication format
- Clear deliverable expectations
- Validation at phase boundaries
- Traceability for debugging

**Principle**: All communication flows through programming-pm. Specialists do not communicate directly.

---

## Base Handoff Format

All handoffs include these base fields.

```yaml
handoff:
  # Metadata
  version: "1.2"
  from_phase: int  # 0-6
  to_phase: int    # 0-6
  producer: string  # skill name that created this
  consumer: string  # skill name that will receive this
  timestamp: ISO8601

  # Session context (from Phase 0)
  session:
    session_dir: string  # Path to /tmp/programming-pm-session-{...}/
    archival_guidelines_path: string  # Path to archival-guidelines-summary.md

  # Deliverable reference
  deliverable:
    location: "/absolute/path/to/file"
    type: "specification" | "code" | "tests" | "documentation"
    checksum: "sha256:..."  # optional, for verification

  # Context for consumer
  context:
    task_id: string
    description: string
    focus_areas: []  # What to pay attention to
    known_gaps: []   # What's incomplete or uncertain

  # Quality assessment
  quality:
    status: "complete" | "partial"
    confidence: "high" | "medium" | "low"
    notes: string  # Explanation of status/confidence
```

---

## Phase 0 -> Phase 1: Archival Setup -> Requirements

**Producer**: programming-pm (Phase 0)
**Consumer**: requirements-analyst (Phase 1)

```yaml
session_handoff:
  <<: *base_handoff

  session:
    session_dir: "/tmp/programming-pm-session-{timestamp}-{pid}/"
    archival_guidelines_path: "{session_dir}/archival-guidelines-summary.md"
    guidelines_found: boolean
    guidelines_source: string  # Path to CLAUDE.md or "defaults"

  archival_guidelines:
    code_directories:
      - name: "src/"
        purpose: "Source code"
      - name: "modules/"
        purpose: "Modular components"
      - name: "experiments/"
        purpose: "Experimental code"
      - name: "models/"
        purpose: "Mathematical/ML models"
    git_workflow:
      commit_after_edit: boolean
      stage_specific_files: boolean
      no_destructive_ops: boolean
      conventional_commits: boolean
    testing_conventions:
      coverage_target: float  # If specified in CLAUDE.md
      test_directory: string
    documentation_conventions:
      docstrings_required: boolean
      type_hints_required: boolean
      readme_updates: boolean
    code_style:
      linter: string  # e.g., "ruff"
      formatter: string  # If specified
```

**Session Context Propagation**:

All downstream agents receive session context via handoffs:
- `requirements-analyst`: Uses session dir for intermediate files
- `systems-architect`: Uses archival guidelines for component naming
- `senior-developer`: Uses git workflow, code style, testing conventions
- `junior-developer`: Same as senior-developer (enforced by review)
- `mathematician`: Uses code directories for model output location
- `statistician`: Uses code directories for analysis output location
- `copilot`: Uses archival guidelines for code review criteria

---

## Phase 1 -> Phase 2: Requirements -> Pre-Mortem

**Producer**: requirements-analyst
**Consumer**: programming-pm (for pre-mortem facilitation)

```yaml
requirements_handoff:
  <<: *base_handoff

  requirements:
    problem_statement: string
    success_criteria: []
    scope:
      in_scope: []
      out_of_scope: []
    constraints: []
    dependencies: []

  stakeholders:
    primary: string
    consulted: []

  approval:
    approved_by: string
    approved_date: ISO8601
    conditions: []  # Any conditional approvals
```

---

## Phase 2 -> Phase 3: Pre-Mortem -> Architecture

**Producer**: programming-pm (pre-mortem results)
**Consumer**: systems-architect

```yaml
premortem_handoff:
  <<: *base_handoff

  risks:
    identified: int  # count
    critical: []     # risk IDs with score >= 15
    high: []         # risk IDs with score 10-14
    mitigated: []    # risks with mitigation plans
    accepted: []     # risks accepted without mitigation

  risk_summary:
    - id: string
      description: string
      score: int
      disposition: "mitigate" | "accept" | "transfer" | "avoid"
      mitigation: string  # if mitigate

  architecture_implications:
    - risk_id: string
      implication: string  # How this affects architecture
```

---

## Phase 3 -> Phase 4: Architecture -> Implementation

**Producer**: systems-architect
**Consumer**: programming-pm (for task decomposition)

```yaml
architecture_handoff:
  <<: *base_handoff

  components:
    - name: string
      responsibility: string
      interfaces:
        inputs: []
        outputs: []
      dependencies: []
      estimated_effort: string  # T-shirt size or hours

  data_flow:
    description: string
    diagram_location: string  # optional

  technology_choices:
    - category: string  # e.g., "database", "web framework"
      choice: string
      rationale: string

  testing_strategy:
    unit: string
    integration: string
    coverage_target: float

  implementation_order:
    - component: string
      priority: int
      dependencies: []

  architecture_context:  # Optional (v1.2+)
    path: string  # Absolute path to .architecture/context.md if generated
    generated: boolean  # Whether context document was created in Phase 3
    version: string  # Template version (e.g., "1.0")
```

---

## mathematician -> developer

**Producer**: mathematician
**Consumer**: senior-developer

```yaml
math_handoff:
  <<: *base_handoff

  algorithm:
    name: string
    description: string
    pseudocode: |
      ...multiline...

  complexity_analysis:
    time:
      best_case: string  # e.g., "O(n)"
      average_case: string
      worst_case: string
    space:
      auxiliary: string
      total: string
    analysis_notes: string

  numerical_stability:
    stable: boolean
    conditions: string
    precision_requirements: string
    failure_modes:
      - condition: string
        symptom: string
        mitigation: string

  implementation_guidance:
    recommended_approach: string
    libraries:
      - name: string
        usage: string
    pitfalls: []

  verification_criteria:
    invariants: []  # Properties that must hold
    test_cases:
      - name: string
        input: string
        expected: string
    edge_cases:
      - name: string
        input: string
        expected: string
        note: string
```

---

## statistician -> developer

**Producer**: statistician
**Consumer**: senior-developer

```yaml
stats_handoff:
  <<: *base_handoff

  method:
    name: string
    description: string
    rationale: string

  assumptions:
    data_requirements: []
    distributional: []
    violations_impact:
      - assumption: string
        impact: string
        mitigation: string

  implementation_guidance:
    library: string
    function: string
    parameters: {}
    code_example: |
      ...multiline...

  power_analysis:  # if applicable
    effect_size: float
    alpha: float
    power: float
    required_n: int
    calculation_method: string

  validation_criteria:
    diagnostic_checks:
      - name: string
        method: string
        threshold: string
    sensitivity_analyses: []

  # For MCMC specifically
  mcmc_config:  # if applicable
    n_chains: int
    warmup: int
    samples: int
    thinning: int
    convergence_criteria:
      ess_threshold: int
      rhat_threshold: float

  interpretation_guide:
    result_format: string
    significant_threshold: float
    interpretation_template: string
```

---

## developer -> code_review

**Producer**: senior-developer or junior-developer
**Consumer**: senior-developer (review) and programming-pm

```yaml
code_handoff:
  <<: *base_handoff

  task:
    id: string
    description: string
    assigned_to: string

  changes:
    files_changed:
      - path: string
        type: "added" | "modified" | "deleted"
        changes: string  # brief description
    lines_added: int
    lines_removed: int

  summary: string  # min 100 chars

  test_coverage:
    new_lines: int
    covered_lines: int
    coverage_percent: float

  self_review_checklist:
    tests_pass: boolean
    ruff_clean: boolean
    mypy_clean: boolean
    documentation_updated: boolean
    type_hints_present: boolean

  open_questions: []
  known_limitations: []

  # If revision
  revision:
    number: int
    changes_made: []
    previous_feedback_addressed: []

  # Architecture context (v1.2+)
  architecture_context:  # Optional
    read: boolean  # Whether developer read .architecture/context.md before implementation
    component_tier: int  # 0=foundation, 1=core, 2=application (from context doc)
    discrepancy_noted: boolean  # True if developer found context-code mismatch
    discrepancy_details: string  # Description of mismatch (required if discrepancy_noted=true)
    stale: boolean  # True if staleness warning was shown during pre-flight
```

---

## code_review -> merge

**Producer**: senior-developer (reviewer)
**Consumer**: programming-pm (for gate decision)

```yaml
review_handoff:
  <<: *base_handoff

  review:
    reviewer: string
    reviewed_date: ISO8601
    status: "approved" | "changes_requested"

  automated_checks:
    ruff: "pass" | "fail"
    mypy: "pass" | "fail"
    tests: "pass" | "fail"
    coverage: float

  manual_review:
    code_quality: "pass" | "issues"
    documentation: "pass" | "issues"
    testing: "pass" | "issues"
    architecture: "pass" | "issues"

  required_changes: []  # if changes_requested

  suggestions: []  # non-blocking

  approval:
    approved: boolean
    conditions: []  # e.g., "Address suggestions in follow-up"
```

---

## Handoff Validation

Before Phase N starts, validate handoff from Phase N-1.

### Validation Steps

1. **Schema validation**: All required fields present
2. **Type validation**: Fields have correct types
3. **Cross-reference check**: References resolve correctly
4. **Consistency check**: No contradictions with previous handoffs

### Validation Script

```python
import yaml
from pathlib import Path

def validate_handoff(handoff_path: Path, schema: dict) -> list[str]:
    """Validate handoff against schema.

    Returns list of validation errors (empty if valid).
    """
    errors = []

    with open(handoff_path) as f:
        handoff = yaml.safe_load(f)

    # Check required fields
    for field in schema.get('required', []):
        if field not in handoff:
            errors.append(f"Missing required field: {field}")

    # Check types
    for field, expected_type in schema.get('types', {}).items():
        if field in handoff and not isinstance(handoff[field], expected_type):
            errors.append(f"Invalid type for {field}: expected {expected_type}")

    return errors
```

### On Validation Failure

1. Log validation errors
2. Present to user with options:
   - Fix issues and re-validate
   - Override with documented gaps
3. Do not proceed until resolved or overridden

---

## Handoff File Naming

```
/tmp/programming-pm-session-{timestamp}-{pid}/   # Session directory (Phase 0)
  ├── archival-guidelines-summary.md             # Phase 0 output
  └── handoffs/
      ├── phase0-session-handoff.yaml
      ├── phase1-requirements-handoff.yaml
      ├── phase2-premortem-handoff.yaml
      ├── phase3-architecture-handoff.yaml
      ├── phase4-math-handoff-TASK-001.yaml
      ├── phase4-stats-handoff-TASK-002.yaml
      ├── phase4-code-handoff-TASK-001.yaml
      ├── phase5-review-handoff-TASK-001.yaml
      └── phase6-merge-handoff.yaml
```

---

## Versioning

Handoff schema version is included in each handoff.

Current version: **1.2**

### Version History

| Version | Changes |
|---------|---------|
| 1.0 | Initial schema |
| 1.1 | Added Phase 0 (Archival Setup), session context in all handoffs, phase range 0-6 |
| 1.2 | Added architecture_context fields to architecture_handoff and code_handoff schemas |

### Compatibility

- Minor version changes (1.x) are backwards compatible
- Major version changes (x.0) may break compatibility
- programming-pm checks version compatibility at validation

### Session Cleanup

On successful Phase 6 completion, programming-pm deletes the session directory:
```bash
rm -rf /tmp/programming-pm-session-{timestamp}-{pid}/
```

On workflow failure or abort, the session directory is retained for debugging. The path is logged to the user for inspection.
