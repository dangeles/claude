# Orchestrator Handoff Schema

Structured YAML handoff contracts for communication between cfd-bioreactor orchestrator
and its specialist agents (cfd-mathematician, cfd-reviewer, brainstorming-pm). Adapted
from the programming-pm handoff schema (version 1.2) with CFD-specific extensions.

**Usage**: Agents copy the exact YAML templates below into their output. The orchestrator
validates required fields at each quality gate. All communication flows through the
orchestrator -- agents never communicate directly.

**Schema Version**: 1.0

---

## 1. Base Handoff Schema

All handoff types include these base fields. CFD-specific types extend this with
additional fields defined in Section 2.

```yaml
handoff:
  # Metadata
  version: "1.0"
  from_phase: int          # 0-5 (matches orchestrator phase number)
  to_phase: int            # 0-5
  producer: string         # "cfd-mathematician" | "cfd-reviewer" | "brainstorming-pm" | "cfd-bioreactor"
  consumer: "cfd-bioreactor"  # Always the orchestrator (hub-and-spoke)
  timestamp: ISO8601       # e.g. "2026-02-21T14:30:00Z"

  # Deliverable reference
  deliverable:
    location: "/absolute/path/to/file"
    type: "specification" | "review" | "synthesis" | "code" | "validation"
    checksum: "sha256:..."  # Optional, for file-based deliverables

  # Context for consumer
  context:
    task_id: string         # Unique identifier for this task invocation
    description: string     # Brief human-readable summary of what was produced
    focus_areas: []         # What the consumer should pay attention to
    known_gaps: []          # What is incomplete or uncertain

  # Quality assessment
  quality:
    status: "complete" | "partial"
    confidence: "high" | "medium" | "low"
    notes: string           # Explanation of status/confidence
```

---

## 2. CFD-Specific Handoff Types

### 2a. math-analysis Handoff

Produced by cfd-mathematician. Contains the mathematical specification for a
simulation phase (mesh, flow, or transport).

**Required fields** (orchestrator will reject if missing):
- `variational_form`
- `function_spaces`
- `convergence_order`
- `dimensionless_numbers`

**Optional fields** (orchestrator proceeds with defaults if missing):
- `stability_conditions`
- `solver_strategy`
- `known_risks`

```yaml
handoff:
  version: "1.0"
  from_phase: 1            # or 2, 3
  to_phase: 1              # Same phase (intra-phase handoff)
  producer: "cfd-mathematician"
  consumer: "cfd-bioreactor"
  timestamp: "2026-02-21T14:30:00Z"
  deliverable:
    location: "{session_dir}/handoffs/phase{N}-math-analysis.yaml"
    type: "specification"
  context:
    task_id: "phase{N}-math-{attempt}"
    description: "Mathematical analysis for phase {N}"
    focus_areas:
      - "variational formulation correctness"
      - "function space stability"
    known_gaps: []
  quality:
    status: "complete"
    confidence: "high"
    notes: ""

  # CFD math-analysis extension
  math_analysis:
    variational_form: |
      Find (u, p) in V x Q such that:
      a((u,p), (v,q)) = L((v,q)) for all (v,q) in V x Q
      where a(...) = mu * inner(grad(u), grad(v)) * dx - p * div(v) * dx + q * div(u) * dx
      and L(...) = inner(f, v) * dx
    function_spaces:
      velocity: "P2 (Lagrange, degree 2)"
      pressure: "P1 (Lagrange, degree 1)"
      justification: "Taylor-Hood pair satisfies inf-sup condition (Babuska-Brezzi)"
    convergence_order:
      velocity_L2: "O(h^3)"
      velocity_H1: "O(h^2)"
      pressure_L2: "O(h^2)"
      basis: "Cea's lemma with a priori interpolation estimates"
    dimensionless_numbers:
      Re: 0.5
      Pe: null               # null if not applicable to this phase
      Da: null
      implications: "Re < 1: Stokes regime; direct solver adequate"
    stability_conditions:     # Optional
      - "inf-sup (Babuska-Brezzi) satisfied by Taylor-Hood P2/P1"
      - "Lax-Milgram applicable to Stokes bilinear form (coercive on kernel of B)"
    solver_strategy:          # Optional
      method: "direct (MUMPS)"
      justification: "DOF count < 50K; direct solver faster than iterative"
      continuation: false
    known_risks: []           # Optional
```

### 2b. engineering-review Handoff

Produced by cfd-reviewer. Contains severity-rated challenges and an approval decision.

**Required fields**:
- `challenges` (list, may be empty only if APPROVED)
- `approval_status`

**Required if REJECTED**:
- `blocking_issues` (list of specific issues that must be resolved)

**Severity levels**: CRITICAL (must fix before proceeding), WARNING (should fix, document if not), NOTE (awareness only)

**Approval statuses**: APPROVED, APPROVED_WITH_WARNINGS, REJECTED

```yaml
handoff:
  version: "1.0"
  from_phase: 1
  to_phase: 1
  producer: "cfd-reviewer"
  consumer: "cfd-bioreactor"
  timestamp: "2026-02-21T14:35:00Z"
  deliverable:
    location: "{session_dir}/handoffs/phase{N}-engineering-review.yaml"
    type: "review"
  context:
    task_id: "phase{N}-review-{attempt}"
    description: "Engineering review for phase {N}"
    focus_areas:
      - "physical plausibility"
      - "numerical stability"
    known_gaps: []
  quality:
    status: "complete"
    confidence: "high"
    notes: ""

  # CFD engineering-review extension
  engineering_review:
    challenges:
      - severity: "WARNING"
        description: "Boundary layer near membrane requires refinement"
        impact: "Under-resolved BL may miss concentration gradients"
        suggested_fix: "Add graded refinement with 5 layers, growth ratio 1.2"
      - severity: "NOTE"
        description: "STEP file units not verified"
        impact: "Mesh may be in mm instead of m"
        suggested_fix: "Add unit check: measure bounding box, compare to expected dimensions"
    approval_status: "APPROVED_WITH_WARNINGS"
    blocking_issues: []      # Required if approval_status is REJECTED
```

### 2c. swarm-synthesis Handoff

Produced by brainstorming-pm (via its standard synthesis output). The orchestrator
extracts these fields from the brainstorming-pm synthesis.

**Required fields**:
- `convergent_insights`
- `divergent_alternatives`
- `confidence_score`

```yaml
handoff:
  version: "1.0"
  from_phase: 1
  to_phase: 1
  producer: "brainstorming-pm"
  consumer: "cfd-bioreactor"
  timestamp: "2026-02-21T14:25:00Z"
  deliverable:
    location: "{session_dir}/handoffs/phase{N}-swarm-synthesis.yaml"
    type: "synthesis"
  context:
    task_id: "phase{N}-swarm"
    description: "Multi-perspective analysis for phase {N} decision"
    focus_areas: []
    known_gaps: []
  quality:
    status: "complete"
    confidence: "medium"
    notes: ""

  # CFD swarm-synthesis extension
  swarm_synthesis:
    convergent_insights:
      - "All perspectives agree: graded mesh refinement near membrane is essential"
      - "Consensus on structured hex elements for channel region"
    divergent_alternatives:
      - "Innovator suggests adaptive mesh refinement (AMR) instead of a priori grading"
      - "Pragmatist recommends uniform mesh with post-hoc convergence study"
    confidence_score: 4      # 1-5 scale; below 3 triggers user notification
```

### 2d. mesh-plan Handoff

Produced by the orchestrator after Phase 1 quality gate passes. Aggregates
mathematician + reviewer + swarm inputs into the approved mesh plan.

**Required fields**:
- `element_type`
- `element_order`
- `refinement_zones`
- `quality_thresholds`
- `memory_estimate_mb`

```yaml
handoff:
  version: "1.0"
  from_phase: 1
  to_phase: 2
  producer: "cfd-bioreactor"
  consumer: "cfd-bioreactor"
  timestamp: "2026-02-21T14:40:00Z"
  deliverable:
    location: "{session_dir}/handoffs/phase1-mesh-plan.yaml"
    type: "specification"
  context:
    task_id: "phase1-mesh-plan"
    description: "Approved mesh plan for code generation"
    focus_areas:
      - "element sizes near boundaries"
      - "memory feasibility"
    known_gaps: []
  quality:
    status: "complete"
    confidence: "high"
    notes: "Reviewer APPROVED_WITH_WARNINGS; 1 WARNING documented"

  # CFD mesh-plan extension
  mesh_plan:
    element_type: "triangle"       # "triangle" | "tetrahedron" | "quadrilateral" | "hexahedron"
    element_order: 1               # Geometric element order (1 = linear, 2 = quadratic)
    refinement_zones:
      - region: "membrane_surface"
        target_size: 0.0005        # meters
        growth_ratio: 1.2
        layers: 5
      - region: "channel_bulk"
        target_size: 0.002
    quality_thresholds:
      min_jacobian: 0.1
      max_aspect_ratio: 10.0
      min_angle_deg: 15.0
    memory_estimate_mb: 450
    estimated_elements: 25000
    estimated_dofs: 150000         # After function space application
```

### 2e. flow-result Handoff

Produced by the orchestrator after Phase 2 execution and validation.

**Required fields**:
- `solver_used`
- `convergence_achieved`
- `mass_conservation_error`
- `validation_metrics`

```yaml
handoff:
  version: "1.0"
  from_phase: 2
  to_phase: 3
  producer: "cfd-bioreactor"
  consumer: "cfd-bioreactor"
  timestamp: "2026-02-21T15:00:00Z"
  deliverable:
    location: "{session_dir}/handoffs/phase2-flow-result.yaml"
    type: "validation"
  context:
    task_id: "phase2-flow-result"
    description: "Flow solver execution results and validation"
    focus_areas:
      - "mass conservation"
      - "velocity field quality for transport"
    known_gaps: []
  quality:
    status: "complete"
    confidence: "high"
    notes: ""

  # CFD flow-result extension
  flow_result:
    solver_used: "MUMPS (direct)"
    convergence_achieved: true
    newton_iterations: null         # null for Stokes (linear); integer for NS
    mass_conservation_error: 1.2e-12
    validation_metrics:
      poiseuille_L2_error: 0.003    # null if not a Poiseuille case
      max_velocity: 0.015           # m/s
      pressure_drop: 12.5           # Pa
    output_files:
      velocity_xdmf: "{session_dir}/scripts/velocity.xdmf"
      pressure_xdmf: "{session_dir}/scripts/pressure.xdmf"
```

### 2f. transport-result Handoff

Produced by the orchestrator after Phase 3 execution and validation.

**Required fields**:
- `stabilization_method`
- `regularization_params`
- `species_conservation_error`
- `min_concentration`

```yaml
handoff:
  version: "1.0"
  from_phase: 3
  to_phase: 4
  producer: "cfd-bioreactor"
  consumer: "cfd-bioreactor"
  timestamp: "2026-02-21T15:20:00Z"
  deliverable:
    location: "{session_dir}/handoffs/phase3-transport-result.yaml"
    type: "validation"
  context:
    task_id: "phase3-transport-result"
    description: "Transport solver execution results and validation"
    focus_areas:
      - "species conservation"
      - "negative concentration check"
    known_gaps: []
  quality:
    status: "complete"
    confidence: "high"
    notes: ""

  # CFD transport-result extension
  transport_result:
    stabilization_method: "SUPG with conditional xi formula"
    regularization_params:
      method: "sqrt regularization"
      epsilon: 1.0e-10              # eps = 1e-10 * c_inlet
      formula: "c_pos = (c + sqrt(c^2 + eps^2)) / 2"
    species_conservation_error: 3.5e-8
    min_concentration: -2.1e-11     # Slightly negative due to numerical error
    negative_concentration_warning: true  # true if min_concentration < 0
    newton_iterations: 8
    output_files:
      concentration_xdmf: "{session_dir}/scripts/oxygen.xdmf"
```

---

## 3. Handoff Validation Rules

The orchestrator validates handoffs at each quality gate using these rules.

### Required vs Optional Fields

| Handoff Type | Required Fields | Optional Fields |
|---|---|---|
| math-analysis | variational_form, function_spaces, convergence_order, dimensionless_numbers | stability_conditions, solver_strategy, known_risks |
| engineering-review | challenges, approval_status | blocking_issues (required if REJECTED) |
| swarm-synthesis | convergent_insights, divergent_alternatives, confidence_score | (none) |
| mesh-plan | element_type, element_order, refinement_zones, quality_thresholds, memory_estimate_mb | estimated_elements, estimated_dofs |
| flow-result | solver_used, convergence_achieved, mass_conservation_error, validation_metrics | newton_iterations, output_files |
| transport-result | stabilization_method, regularization_params, species_conservation_error, min_concentration | negative_concentration_warning, newton_iterations, output_files |

### Missing Field Behavior

**Missing required field**:
1. Re-invoke the producing agent with explicit feedback: "Your handoff is missing required field: {field_name}. Please include it."
2. Maximum 1 retry. If still missing after retry: log warning and proceed with a conservative default (see below).

**Missing optional field**:
- Proceed without the field. Log an informational note.

**Conservative defaults for missing required fields**:
- Missing `approval_status` on engineering-review: treat as REJECTED (conservative)
- Missing `convergence_achieved` on flow-result: treat as false
- Missing `min_concentration` on transport-result: treat as negative (triggers warning)
- Missing `confidence_score` on swarm-synthesis: treat as 1 (lowest confidence)

### Lenient Parsing

The orchestrator checks common synonyms before reporting a missing field:

| Canonical Name | Accepted Synonyms |
|---|---|
| `approval_status` | `review_status`, `status`, `decision` |
| `challenges` | `issues`, `findings`, `concerns` |
| `variational_form` | `weak_form`, `variational_formulation` |
| `convergence_order` | `convergence_rate`, `error_order` |
| `blocking_issues` | `blockers`, `critical_issues` |

If a synonym is found, the orchestrator uses it but logs a note: "Agent used non-canonical field name '{synonym}' for '{canonical}'. Consider updating agent output."

---

## 4. Canonical Field Names

Agents MUST use these exact field names in their handoff YAML. The lenient parsing
in Section 3 is a fallback, not a license to use non-canonical names.

| Canonical Name | Type | Used By |
|---|---|---|
| `variational_form` | string (multiline) | math-analysis |
| `function_spaces` | object | math-analysis |
| `convergence_order` | object | math-analysis |
| `dimensionless_numbers` | object | math-analysis |
| `stability_conditions` | list of strings | math-analysis |
| `solver_strategy` | object | math-analysis |
| `known_risks` | list of strings | math-analysis |
| `challenges` | list of objects | engineering-review |
| `approval_status` | enum string | engineering-review |
| `blocking_issues` | list of strings | engineering-review |
| `convergent_insights` | list of strings | swarm-synthesis |
| `divergent_alternatives` | list of strings | swarm-synthesis |
| `confidence_score` | integer (1-5) | swarm-synthesis |
| `element_type` | enum string | mesh-plan |
| `element_order` | integer | mesh-plan |
| `refinement_zones` | list of objects | mesh-plan |
| `quality_thresholds` | object | mesh-plan |
| `memory_estimate_mb` | number | mesh-plan |
| `solver_used` | string | flow-result |
| `convergence_achieved` | boolean | flow-result |
| `mass_conservation_error` | number | flow-result |
| `validation_metrics` | object | flow-result |
| `stabilization_method` | string | transport-result |
| `regularization_params` | object | transport-result |
| `species_conservation_error` | number | transport-result |
| `min_concentration` | number | transport-result |
| `negative_concentration_warning` | boolean | transport-result |

---

## 5. Error History Extension

For self-correction loops (Steps 2.5b, 3.5b in the orchestrator), handoffs include
an `error_history` field that tracks previous failed attempts. This prevents agents
from recommending fixes that have already been tried and failed.

```yaml
  # Appended to any handoff during self-correction
  error_history:
    - attempt_number: 1
      error_type: "solver_divergence"
      error_message: "Newton solver did not converge in 50 iterations"
      fix_attempted: "Reduced relaxation parameter to 0.5"
      fix_outcome: "Still diverged at iteration 38"
    - attempt_number: 2
      error_type: "solver_divergence"
      error_message: "Newton solver did not converge in 50 iterations"
      fix_attempted: "Switched to Picard iteration with 100 max iterations"
      fix_outcome: "Converged at iteration 72"
```

**Error type taxonomy**:
- `import_error`: Missing module or wrong import path
- `mesh_error`: gmsh failure, invalid geometry, quality threshold violation
- `solver_divergence`: Newton/Picard did not converge within iteration limit
- `assertion_failure`: Version check, conservation check, or quality check failed
- `numerical_instability`: NaN/Inf in solution, negative concentrations beyond tolerance

**Usage by agents**:
- cfd-reviewer in error diagnosis mode reads `error_history` to avoid suggesting previously-failed fixes
- Orchestrator appends each new attempt to the list before re-invoking agents
- Maximum 3 entries (after 3 failed attempts, escalate to user)
