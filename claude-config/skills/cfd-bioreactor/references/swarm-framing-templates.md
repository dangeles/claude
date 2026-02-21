# Swarm Framing Templates

Pre-written challenge templates for brainstorming-pm swarm invocations at the three
decision points in the cfd-bioreactor orchestrator workflow. Each template injects
concrete numerical parameters (Re, Pe, Da, element counts, memory budgets) to force
the swarm perspectives to engage with specific CFD trade-offs rather than producing
generic advice.

**Usage**: The orchestrator reads the appropriate template below, fills in the
placeholder variables `{...}` with values from the problem specification and
preceding phase results, then passes the filled template as the challenge to
brainstorming-pm via the Task tool.

**Swarm invocation**: FULL mode only. LITE and DIRECT modes skip all swarm steps.

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

## 4. Standard Archetypes and CFD-Specific Framing

The swarm uses the standard 5 archetypes from brainstorming-pm / perspective-swarm.
The CFD-specific framing comes from the challenge templates above, not from custom
archetypes. Each archetype naturally maps to a useful CFD perspective:

| Archetype | Natural CFD Perspective | Typical Contribution |
|---|---|---|
| Optimist | Best-case performance advocate | "With proper refinement, we can achieve 2nd-order convergence" |
| Critic | Failure mode detector | "Pe > 710 will cause SUPG parameter overflow with cosh/sinh" |
| Analyst | Quantitative trade-off evaluator | "Direct solver needs 2.1 GB RAM for 80K DOFs; iterative needs 0.4 GB" |
| Innovator | Alternative approach proposer | "Consider MINI element instead of Taylor-Hood for this low-Re case" |
| Pragmatist | Implementation feasibility assessor | "Newton with Stokes initial guess converges in 5 iterations for Re < 5" |

**No custom archetypes in v2.0.** The concrete numerical framing in the challenge
templates provides sufficient CFD domain context. Custom archetypes (Theoretician,
Experimentalist, Computationalist, Skeptic, Practitioner) are a Could-Have for
future iterations.

---

## 5. Swarm Output Handling

The orchestrator extracts structured information from the brainstorming-pm synthesis
output and converts it into a swarm-synthesis handoff (see orchestrator-handoff-schema.md
Section 2c).

### Extraction Protocol

1. **Read** the brainstorming-pm synthesis output (the final convergence analysis).
2. **Identify convergent insights**: statements where 3+ perspectives agree. Extract as a list of concrete recommendations.
3. **Identify divergent alternatives**: unique suggestions from individual perspectives that were not adopted by the majority. Extract as a list of alternative approaches.
4. **Assess confidence**: Map the brainstorming-pm confidence level to a 1-5 integer score:
   - 5 = "high confidence, strong consensus" (4-5 perspectives agree on key points)
   - 4 = "good confidence, majority consensus" (3 perspectives agree)
   - 3 = "moderate confidence, split opinions" (2-3 perspectives agree, significant dissent)
   - 2 = "low confidence, no clear consensus" (perspectives mostly disagree)
   - 1 = "very low confidence, contradictory" (perspectives contradict each other)

### Minimum Quality Criteria

The swarm synthesis is usable if it meets ALL of:
1. Contains at least 2 specific numerical recommendations (element count, refinement ratio, solver parameter, convergence criterion, stabilization value)
2. Contains at least 1 concrete alternative approach (not just "consider alternatives" -- must name a specific method or parameter change)

### Below-Threshold Behavior

If the synthesis fails the quality criteria:
- Log: "Swarm output for Phase {N} below quality threshold (found {count} numerical recommendations, needed 2)."
- Proceed WITHOUT swarm input. The mathematician and reviewer operate independently.
- The orchestrator notes in the quality section of downstream handoffs: "Swarm synthesis omitted due to insufficient specificity."
- This is NOT an error condition. The swarm is advisory; the mathematician and reviewer provide the primary analysis.
