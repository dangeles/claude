# Swarm Framing Templates

Pre-written challenge templates for FULL mode multi-perspective analysis at the
three decision points in the cfd-bioreactor orchestrator workflow. Each template
injects concrete numerical parameters (Re, Pe, Da, element counts, memory budgets)
to force the perspective agents to engage with specific CFD trade-offs rather than
producing generic advice.

**Usage**: The orchestrator reads the appropriate template below, fills in the
placeholder variables `{...}` with values from the problem specification and
preceding phase results, then passes the filled template to each perspective
agent via parallel Task tool invocations (see SKILL.md Section 8b).

**Swarm invocation**: FULL mode only. LITE and DIRECT modes skip all perspective steps.

---

## 1. Swarm 1: Pre-Mesh Challenge

Invoked at Phase 1, Step 1.1, before mathematical analysis and adversarial review.

### Challenge Template

```
CHALLENGE: Mesh Strategy for Bioreactor Geometry

Geometry: {geometry_description}
Dimensions: {dimensions} ({dimension_type})
Physical groups required: {physical_groups}
Available memory: {memory_budget_mb} MB
Target element count range: {element_count_min}-{element_count_max}

Key considerations:
- Boundary layer resolution needed at: {bl_surfaces}
- Expected flow regime: Re ~ {re_estimate}
- CAD complexity: {cad_complexity}

Choose between:
(A) Uniform mesh with h = {uniform_h}, ~{uniform_elements}K elements
(B) Graded mesh: h_min = {graded_h_min} near {bl_surfaces}, h_max = {graded_h_max} in bulk, ~{graded_elements}K elements
(C) Structured boundary layer mesh: {bl_layers} layers with growth ratio {growth_ratio}, unstructured bulk

Evaluate trade-offs: accuracy per DOF, memory usage, mesh generation robustness,
solver conditioning, and downstream impact on flow/transport quality.
```

### Context Injection Checklist

Before filling the template, the orchestrator must have:
- [ ] Geometry dimensions from user request or STEP file bounding box
- [ ] Physical group list (inlet, outlet, walls, membranes)
- [ ] Memory budget (from environment check or user constraint; default: 70% of system RAM)
- [ ] Element count range (estimated from geometry volume and target h)
- [ ] Boundary layer surfaces identified (membranes, narrow channels)
- [ ] Reynolds number estimate (from user-provided flow rate and geometry)
- [ ] CAD complexity assessment (simple parametric vs. imported STEP with many surfaces)

### Quality Threshold

The swarm synthesis for Swarm 1 must contain:
- At least 2 specific numerical recommendations (e.g., "h_min = 0.0005 m near membrane", "target 25K elements")
- At least 1 concrete alternative approach with numerical justification

If below threshold: log "Swarm 1 output below quality threshold; proceeding with mathematician analysis only."

---

## 2. Swarm 2: Pre-Flow Challenge

Invoked at Phase 2, Step 2.1, after mesh plan is approved (Phase 1 complete).

### Challenge Template

```
CHALLENGE: Flow Solver Strategy for Bioreactor Simulation

Geometry: {geometry_description}
Flow regime: Re = {re_value} ({re_regime})
Boundary conditions:
  Inlet: {inlet_bc_type} ({inlet_bc_value})
  Outlet: {outlet_bc_type} ({outlet_bc_value})
  Walls: {wall_bc_type}
  Membranes: {membrane_bc_type}
Mesh: {element_count}K elements, {dof_count}K DOFs ({element_type}, order {element_order})

Evaluate:
(A) Direct solver (MUMPS): robust for < 50K DOFs, high memory
(B) Iterative solver (GMRES + ILU): scales to large problems, may need tuning
(C) Newton iteration (for Re > 1): full nonlinear, fast convergence if good initial guess
(D) Picard iteration: simpler, more robust but slower convergence

Additional considerations:
- If Re > 10: need continuation strategy (Stokes initial guess, Re ramping)
- Pressure reference point: {pressure_reference}
- Function space: Taylor-Hood P2/P1 (default) vs alternatives

Recommend a solver strategy with convergence criteria and fallback plan.
```

### Context Injection Checklist

Before filling the template, the orchestrator must have:
- [ ] Reynolds number (computed from fluid properties and geometry)
- [ ] Boundary condition types and values for all boundaries
- [ ] Mesh size from Phase 1 mesh-plan handoff (element_count, estimated_dofs)
- [ ] Element type and order from Phase 1 mesh-plan handoff
- [ ] Pressure reference point specification (which boundary, what value)
- [ ] Solver options available in the environment (MUMPS presence from pre-flight)

### Quality Threshold

The swarm synthesis for Swarm 2 must contain:
- A specific solver strategy recommendation (not just "consider iterative methods")
- Numerical convergence criteria (e.g., "relative tolerance 1e-8", "max 50 Newton iterations")

If below threshold: log "Swarm 2 output below quality threshold; proceeding with mathematician analysis only."

---

## 3. Swarm 3: Pre-Transport Challenge

Invoked at Phase 3, Step 3.1, after flow solution is validated (Phase 2 complete).

### Challenge Template

```
CHALLENGE: Transport Stabilization and Reaction Model for Species Transport

Species: {species_name} ({species_description})
Transport parameters:
  Diffusivity: D = {diffusivity} m^2/s
  Peclet number: Pe = {pe_value} ({pe_regime})
  Damkohler number: Da = {da_value} ({da_regime})
Reaction model: {reaction_model}
  Michaelis-Menten parameters: V_max = {vmax}, K_m = {km}
Membrane boundary condition: {membrane_bc_type}
  Permeability: P = {membrane_permeability} m/s
  External concentration: c_ext = {external_concentration}
Flow field: max velocity = {max_velocity} m/s from Phase 2

Compare:
(A) Standard SUPG with conditional xi formula (xi = coth(Pe_h/2) - 2/Pe_h, capped)
(B) Residual-based stabilization (GLS or VMS)
(C) DG (discontinuous Galerkin) transport with upwind flux

For reaction model:
(A) Direct Michaelis-Menten with sqrt regularization: c_pos = (c + sqrt(c^2 + eps^2)) / 2
(B) Linearized reaction at each Newton step
(C) Operator splitting: solve transport and reaction separately

Address stabilization parameter overflow risk if Pe > 710 (cosh/sinh overflow).
Address negative concentration risk and mitigation strategies.
```

### Context Injection Checklist

Before filling the template, the orchestrator must have:
- [ ] Peclet number (computed from max velocity, characteristic length, diffusivity)
- [ ] Damkohler number (computed from V_max, K_m, c_inlet, characteristic time)
- [ ] Michaelis-Menten parameters (V_max, K_m from user request or parameter tables)
- [ ] Membrane permeability and external concentration (if Robin BC)
- [ ] Maximum velocity from Phase 2 flow-result handoff
- [ ] Element size from Phase 1 mesh-plan handoff (for element Peclet number)
- [ ] Diffusivity from species parameters

### Quality Threshold

The swarm synthesis for Swarm 3 must contain:
- A specific stabilization parameter recommendation (not just "use SUPG")
- Explicit statement on overflow risk mitigation for the given Pe
- At least 1 concrete alternative to the default approach

If below threshold: log "Swarm 3 output below quality threshold; proceeding with mathematician analysis only."

---

## 6. CFD-Specific Perspective Definitions

The decomposed pipeline uses domain-specific perspectives instead of generic
brainstorming archetypes. Each perspective is a Task-spawned sub-agent with a
focused analysis lens.

### 6a. Perspective Prompts

**Numerical Analyst:**
```
You are an expert in numerical analysis for finite element methods. Your focus
is mathematical stability, convergence properties, and error estimation.
Evaluate the challenge through the lens of: Will the numerical method converge?
What are the expected error bounds? Are stability conditions satisfied?
Cite specific theorems (Lax-Milgram, Babuska-Brezzi, Cea's lemma) when applicable.
```

**Mesh Engineer:**
```
You are an expert in computational mesh generation for CFD. Your focus is mesh
quality, element sizing, memory efficiency, and generation robustness.
Evaluate the challenge through the lens of: Is the mesh adequate for the physics?
Are boundary layers resolved? Is the memory budget feasible? What are the mesh
quality metrics (scaled Jacobian, aspect ratio)?
Provide specific numbers for element sizes, growth ratios, and element counts.
```

**Physical Modeler:**
```
You are an expert in bioreactor fluid dynamics and mass transport. Your focus
is physical accuracy, model assumptions, and parameter plausibility.
Evaluate the challenge through the lens of: Are the physics correctly captured?
Are model assumptions valid for this operating regime? Are parameter values
physically reasonable? What simplifications might cause problems?
Reference dimensionless numbers (Re, Pe, Da) and their regime implications.
```

**Computational Pragmatist:**
```
You are an expert in practical CFD computation. Your focus is solver efficiency,
runtime estimation, memory management, and robustness.
Evaluate the challenge through the lens of: What solver strategy minimizes
wall-clock time? What are the memory requirements? What happens if the solver
diverges? What is the fallback plan?
Provide specific solver recommendations with expected iteration counts and runtimes.
```

**Validation Strategist:**
```
You are an expert in CFD verification and validation. Your focus is benchmark
comparisons, convergence testing, and uncertainty quantification.
Evaluate the challenge through the lens of: How will we know the result is correct?
What benchmarks should we compare against? What convergence behavior should we
expect? What are the sources of uncertainty?
Reference specific analytical solutions and expected error magnitudes.
```

### 6b. Phase-Perspective Mapping

| Phase | Required Perspectives (min 3) | Optional (if budget allows) |
|-------|-------------------------------|----------------------------|
| Phase 1 (Mesh) | Mesh Engineer, Numerical Analyst, Computational Pragmatist | Physical Modeler, Validation Strategist |
| Phase 2 (Flow) | Numerical Analyst, Physical Modeler, Computational Pragmatist | Mesh Engineer, Validation Strategist |
| Phase 3 (Transport) | Numerical Analyst, Physical Modeler, Validation Strategist | Computational Pragmatist, Mesh Engineer |

### 6c. Perspective Reference Material

Each perspective agent receives relevant reference excerpts in its Task prompt.
The orchestrator loads these sections and includes them inline (perspectives
agents CAN use Read tool, but inline delivery is more reliable):

| Perspective | Reference Sections to Include |
|-------------|------------------------------|
| Numerical Analyst | physics-models.md Sections 1-3 (equations), Section 6 (dimensionless numbers) |
| Mesh Engineer | mesh-generation-guide.md Sections 4-6 (refinement, memory, quality) |
| Physical Modeler | physics-models.md Section 5 (parameter tables), validation-benchmarks.md Benchmark 1 |
| Computational Pragmatist | troubleshooting-guide.md Stage 4 (solver errors), fenicsx-patterns.md Section 8 (solver config) |
| Validation Strategist | validation-benchmarks.md Sections 1-4 (benchmarks, convergence protocol) |

---

## 5. Swarm Output Handling

### Current Protocol (v2.0 Decomposed Pipeline)

The synthesis agent produces a structured swarm-synthesis handoff YAML directly
(see SKILL.md Section 8b and orchestrator-handoff-schema.md Section 2c). No
extraction from prose is needed.

The orchestrator validates the synthesis handoff:
1. Parse YAML and check for required fields: `convergent_insights`,
   `divergent_alternatives`, `confidence_score`
2. Apply minimum quality criteria (see below)
3. If valid: pass to mathematician as supplementary context
4. If invalid or missing: proceed without swarm input

### Minimum Quality Criteria

The swarm synthesis is usable if it meets ALL of:
1. Contains at least 2 specific numerical recommendations in `convergent_insights`
   (element count, refinement ratio, solver parameter, convergence criterion,
   stabilization value)
2. Contains at least 1 concrete alternative approach in `divergent_alternatives`
   (not just "consider alternatives" -- must name a specific method or parameter change)

### Below-Threshold Behavior

If the synthesis fails the quality criteria:
- Log: "Swarm output for Phase {N} below quality threshold (found {count}
  numerical recommendations, needed 2)."
- Proceed WITHOUT swarm input. The mathematician and reviewer operate independently.
- The orchestrator notes in downstream handoffs: "Swarm synthesis omitted due
  to insufficient specificity."
- This is NOT an error condition. The swarm is advisory; the mathematician and
  reviewer provide the primary analysis.

### Legacy Protocol (Reference Only)

The following extraction protocol was used in the original brainstorming-pm
integration design and is retained for reference:

1. Read the brainstorming-pm synthesis output (stage-3-synthesis.md)
2. Identify convergent insights: statements where 3+ perspectives agree
3. Identify divergent alternatives: unique suggestions from individual perspectives
4. Assess confidence and map to 1-5 integer score
5. Construct swarm-synthesis handoff YAML

This protocol is no longer used because the decomposed pipeline (Section 8b
of SKILL.md) produces structured YAML output directly via the synthesis agent.
