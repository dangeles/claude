---
name: cfd-bioreactor
description: >
  Use when user needs to simulate fluid flow through bioprocess cartridges,
  bioreactor geometries, or membrane devices using FEniCSx. Orchestrates
  cfd-mathematician, cfd-reviewer, and brainstorming-pm agents through a
  5-phase workflow with quality gates, self-correction loops, and three-tier
  mode selection (DIRECT/LITE/FULL).

# v2.0 orchestrator handoff metadata
handoff:
  accepts_from:
    - "*"
  provides_to:
    - cfd-mathematician
    - cfd-reviewer
    - brainstorming-pm
  schema_version: "1.0"
  schema_type: cfd-orchestrator

categories:
  - scientific-simulation
  - orchestrator
---

# CFD Bioreactor Simulation Orchestrator

Generate complete, runnable FEniCSx Python scripts for research-grade CFD simulation
of bioprocess cartridge designs. This skill covers geometry import (STEP/IGES), mesh
generation, incompressible Navier-Stokes flow solving, species transport with O2
consumption (Michaelis-Menten), membrane permeation (Fick's law), and interactive 3D
visualization via PyVista.

All generated code targets FEniCSx v0.10 and uses the two-phase segregated solve
architecture: Phase A solves fluid flow, Phase B uses the velocity field for species
transport. Each phase is independently validatable.

---

## 1. Delegation Mandate

You are an **orchestrator**. You coordinate specialists -- you do not perform specialist
work yourself.

**Delegate to specialists via Task tool:**
- **Mathematical analysis** --> `cfd-mathematician` (variational formulations, stability,
  convergence, dimensionless analysis)
- **Engineering review** --> `cfd-reviewer` (adversarial challenge of plans, error
  diagnosis, physical plausibility)
- **Brainstorming** --> `brainstorming-pm` (multi-perspective discussion at decision
  points; FULL mode only)

**Orchestrator-owned tasks** (do NOT delegate these):
- Session setup, state management, phase transitions
- Quality gate evaluation and user communication
- Pre-flight validation
- **Code generation**: Code generation is a mechanical translation from validated
  mathematical specifications into FEniCSx Python using parameterized patterns from
  reference files. It is template-filling, not creative specialist work. The
  mathematician produces the specification; you fill in the code template.

**Graceful degradation**: If a required agent skill is unavailable (SKILL.md not found
at deployed path), degrade gracefully:
- Missing `cfd-mathematician`: Use reference file defaults for mathematical decisions
- Missing `cfd-reviewer`: Proceed with APPROVED_WITH_WARNINGS, log "NO ADVERSARIAL REVIEW"
- Missing `brainstorming-pm`: Skip swarm steps (equivalent to LITE mode)

---

## 2. When to Use This Skill

Use this skill when the user needs to:

- Simulate fluid flow through bioprocess cartridges, bioreactors, or membrane devices
- Generate FEniCSx Python scripts for Navier-Stokes + species transport
- Model O2 consumption with Michaelis-Menten kinetics
- Model membrane permeation with Fick's law Robin boundary conditions
- Import CAD geometry (STEP/IGES) for mesh generation via gmsh
- Validate CFD results against analytical solutions (Poiseuille flow, 1D diffusion-reaction)
- Create interactive 3D visualizations of flow and concentration fields
- Run mesh convergence studies for publication-quality results

**Trigger phrases**: "simulate flow through", "bioreactor CFD", "Navier-Stokes",
"oxygen transport", "FEniCSx", "mesh from STEP file", "species transport",
"Michaelis-Menten", "membrane permeation", "flow visualization"

---

## 3. When NOT to Use This Skill

Do NOT use this skill for:

- **Compressible flow or Mach effects**: Not relevant for bioreactors (Ma << 0.01)
- **Turbulence modeling**: Laminar flow assumed (Re < 100). For turbulent flows, suggest OpenFOAM
- **Structural mechanics or fluid-structure interaction**: Use a dedicated FEM/FSI tool
- **HPC job submission or cloud orchestration**: This skill runs locally or provides scripts
- **Thermal effects**: Isothermal assumption throughout. Heat transfer is out of scope
- **Commercial solver workflows**: COMSOL, ANSYS, OpenFOAM have their own interfaces
- **Quick order-of-magnitude estimates**: Use the `calculator` skill instead. A back-of-envelope
  calculation takes 5 minutes; a full CFD simulation takes 30 minutes to hours
- **Transient simulations with adaptive time-stepping**: Steady-state only in v2.0

**If in doubt**: Start with the `calculator` skill for a quick feasibility estimate,
then come here for the full CFD simulation.

---

## 4. State Anchoring

Prefix every response with the current phase and status:

```
[Phase N/5 - {phase_name}] {brief status}
```

**Examples**:
- `[Phase 0/5 - Pre-Flight] Environment validated. Selecting mode.`
- `[Phase 2/5 - Flow Planning] Mathematician analysis received. Invoking reviewer.`
- `[Phase 2/5 - Flow Execution] Self-correction loop: attempt 2 of 3.`

**Rules**:
- Before starting any phase: Read state file, confirm `current_phase`
- After any user interaction: Answer the question, then re-anchor with phase status
- During FULL mode with parallel agents: Show status board

Status board format (FULL mode):
```
Phase {N} Agent Status:
  Swarm:        [DONE | RUNNING | SKIPPED]
  Mathematician: [DONE | RUNNING | PENDING]
  Reviewer:      [DONE | RUNNING | PENDING]
```

---

## 5. Mode Selection -- Three-Tier System

### 5a. DIRECT Mode (Tier 1 Only)

Skip ALL agents. Use the direct-generation approach for well-understood validation
benchmarks.

- **Triggered by**: Tier 1 Poiseuille validation OR user override to DIRECT
- **Workflow**: Phase 0 --> direct code gen for mesh + flow --> validation --> Phase 5
- **Time estimate**: ~10 minutes
- **Rationale**: Well-understood validation benchmark; agents add no value

In DIRECT mode, suppress all agent-related status messages. The user sees a simple
linear workflow identical to v1.0 behavior.

### 5b. LITE Mode (Tier 2, Default for 2D)

Skip swarm discussions. Use cfd-mathematician + cfd-reviewer sequentially at each
decision point.

- **Triggered by**: Tier 2 2D problems without hidden complexity flags
- **Workflow**: Phase 0 --> Phases 1-3 (mathematician + reviewer, no swarm) --> Phase 4-5
- **Time estimate**: ~20-25 minutes

### 5c. FULL Mode (Tier 3-4, Default for 3D)

Full swarm + mathematician + reviewer at every decision point.

- **Triggered by**: Tier 3-4 OR user override OR hidden complexity detected in Tier 2
- **Workflow**: Phase 0 --> Phases 1-3 (swarm + mathematician + reviewer) --> Phase 4-5
- **Time estimate**: ~40-60 minutes

### 5d. Complexity Heuristics for Mode Upgrade

During Phase 0 or during any reviewer analysis, upgrade the mode if:

| Indicator | Threshold | Recommendation |
|---|---|---|
| Peclet number | Pe > 100 | Recommend FULL even for 2D |
| Membrane interfaces | > 1 membrane surface | Recommend FULL |
| Reynolds number | Re > 10 | Recommend FULL (Newton continuation needed) |
| Channel geometry | Width < 10x boundary layer thickness estimate | Recommend FULL |

Display the recommendation with rationale. The user can override.

### 5e. Mode Display Template

Present mode selection to user at end of Phase 0:

```
Problem: [description] (Tier [N])
Available modes: DIRECT (10 min) | LITE (20 min) | FULL (40 min)
Recommended: [MODE] ([reason])
Override: You can select a different mode.
```

---

## 6. Tool Selection Table

| Situation | Tool | Notes |
|---|---|---|
| Generate Python simulation script | Write | ALWAYS complete, runnable scripts. Never fragments. |
| Execute 2D simulation (< 5 min) | Bash (timeout=300000) | Set 5-minute timeout. Fast enough for inline. |
| Execute 3D simulation (> 5 min) | Write only | Save script. User runs manually. Too slow for Bash. |
| Check simulation output / logs | Read | Examine error messages, convergence reports. |
| Load reference patterns for code gen | Read | Load fenicsx-patterns.md, physics-models.md, etc. |
| Quick analytical feasibility estimate | Skill(calculator) | Before full CFD, estimate Re, Pe, Da numbers. |
| Look up cell-type parameters | Read | Load physics-models.md Section 5 (parameter tables). |
| Debug FEniCSx import errors | Read | Load troubleshooting-guide.md Stage 0. |
| Delegate mathematical analysis | Task(cfd-mathematician) | Invoke with problem spec + reference loading instructions. |
| Delegate engineering review | Task(cfd-reviewer) | Invoke with plan + mathematician output + error history. |
| Delegate brainstorming | Task(brainstorming-pm) | Invoke with filled challenge template. FULL mode only. |

---

## 7. Pre-Flight Validation Protocol -- Phase 0

**MANDATORY**: Before executing ANY workflow step, run pre-flight validation.

### Step 0.1: Environment Validation

1. Read `references/environment-setup.md` Section 4 (Pre-Flight Validation Script)
2. Generate `preflight_check.py` and execute it via Bash
3. Interpret results per decision table:

| Pre-Flight Result | Action |
|---|---|
| All PASS | Proceed to workflow |
| FAIL on dolfinx | Provide conda/Docker install instructions from environment-setup.md. STOP. |
| FAIL on gmsh OCC | Can still do parametric geometry. Warn user that STEP import is disabled. |
| WARN on PyVista | Proceed; export to VTK instead of interactive plots. |
| WARN on MUMPS | Proceed; use iterative solver (GMRES+ILU) instead of direct solver. |
| FAIL on petsc4py | Cannot solve. Provide install instructions. STOP. |
| WARN on RAM < 8 GB | Proceed for 2D only. Warn that 3D will require coarse meshes. |

**Critical Rule**: If pre-flight fails on dolfinx or petsc4py, do NOT attempt to generate
simulation scripts. Help the user set up the environment first.

### Step 0.2: Agent Availability Checks

Verify specialist agents are deployed:

1. Check `cfd-mathematician` SKILL.md exists at deployed path
2. Check `cfd-reviewer` SKILL.md exists at deployed path
3. Check `brainstorming-pm` SKILL.md exists (required only for FULL mode)
4. If agents missing: "New agent skills not synced. Run: `./sync-config.py push`"
5. If agents remain unavailable: degrade gracefully (see Section 1)

### Step 0.3: Reference File Checks

Verify all reference files exist:
- `references/environment-setup.md`
- `references/physics-models.md`
- `references/mesh-generation-guide.md`
- `references/fenicsx-patterns.md`
- `references/validation-benchmarks.md`
- `references/troubleshooting-guide.md`
- `references/orchestrator-handoff-schema.md`
- `references/agent-loading-guide.md`
- `references/swarm-framing-templates.md`

### Step 0.4: Session Initialization

1. Create session directory: `/tmp/cfd-bioreactor-session-{YYYYMMDD-HHMMSS}-{PID}/`
2. Create subdirectories: `handoffs/`, `scripts/`, `logs/`
3. Initialize state file (see Section 20 for schema)

### Step 0.5: Problem Parsing and Mode Selection

1. Parse user request: extract geometry type, fluid properties, species parameters,
   boundary conditions, target tier
2. Compute initial dimensionless estimates (Re, Pe, Da) if parameters available
3. Select mode based on tier + complexity heuristics (Section 5)
4. Present mode display to user (Section 5e)

**Quality Gate 0**: Pre-flight passes (dolfinx and petsc4py available). Session directory
created. State file initialized. Mode selected.

---

## 8. Phase 1 -- Mesh Planning

### Step 1.1: Swarm Discussion (FULL Mode Only)

1. Read `references/swarm-framing-templates.md` Swarm 1 template
2. Fill in placeholders with problem parameters (geometry, physical groups, memory budget)
3. Invoke `brainstorming-pm` via Task tool with the filled challenge template
4. Check swarm output against quality threshold:
   - At least 2 specific numerical recommendations
   - At least 1 concrete alternative approach
5. If below threshold: log "Swarm 1 output below quality threshold; proceeding with
   mathematician analysis only"

### Step 1.2: Mathematical Analysis

Invoke `cfd-mathematician` via Task tool:
- **Input**: Geometry description, target error tolerance, element type options,
  swarm synthesis (if FULL mode)
- **Loading instructions**: Include per agent-loading-guide.md (physics-models.md
  Sections 1-3, 6; validation-benchmarks.md Benchmark 1)
- **Task**: Recommend element type/order, estimate required mesh density, analyze
  boundary layer resolution needs
- **Timeout**: 8 minutes
- **On timeout**: Use reference file defaults (linear triangles/tetrahedra, h estimated
  from geometry scale)

### Step 1.3: Adversarial Review

Invoke `cfd-reviewer` via Task tool:
- **Input**: Mesh plan from Steps 1.1-1.2, mathematician handoff output
- **Loading instructions**: Include per agent-loading-guide.md
- **Task**: Challenge mesh quality thresholds, refinement adequacy, memory feasibility,
  STEP unit assumptions
- **Timeout**: 8 minutes
- **On timeout**: Proceed with APPROVED_WITH_WARNINGS, log "NO ADVERSARIAL REVIEW"

### Quality Gate 1

**Mechanical checks**:
- Handoff YAML exists with required fields
- Memory estimate < 70% available RAM
- Element count > 0
- Mesh quality threshold (min scaled Jacobian > 0.1) is specified

**Agent assessment**:
- Reviewer `approval_status` is APPROVED or APPROVED_WITH_WARNINGS
- If APPROVED_WITH_WARNINGS: no CRITICAL-severity challenges
- If REJECTED: pass reviewer `blocking_issues` as HARD CONSTRAINTS to mathematician
  (retry Step 1.2)
- Max **1 retry** of mathematician after rejection. If still REJECTED: escalate to user
  (see Section 16)

**Handoff**: Write Phase 1 mesh-plan handoff to `{session_dir}/handoffs/phase1-mesh-plan.yaml`

---

## 9. Phase 2 -- Flow Solver Planning, Code Generation, and Execution

### Step 2.1: Swarm Discussion (FULL Mode Only)

Same pattern as Phase 1 with Swarm 2 template from `references/swarm-framing-templates.md`.

### Step 2.2: Mathematical Analysis

Invoke `cfd-mathematician` via Task tool:
- **Input**: Flow parameters (Re, fluid properties), mesh plan from Phase 1,
  swarm synthesis (if FULL mode)
- **Task**: Write variational formulation (weak form), verify Taylor-Hood inf-sup
  stability, estimate convergence order, recommend continuation strategy if Re > 10
- **Timeout**: 8 minutes

### Step 2.3: Adversarial Review

Invoke `cfd-reviewer` via Task tool:
- **Input**: Flow plan from Steps 2.1-2.2, mathematician handoff output
- **Task**: Challenge BC physical consistency, solver convergence criteria, check for
  missing pressure reference point, verify Newton vs. Picard choice, check memory estimates
- **Timeout**: 8 minutes

### Quality Gate 2a: Flow Plan Approved

Same pattern as Quality Gate 1:
- Mechanical checks + reviewer approval
- Max 1 retry on rejection, then escalate to user

### Step 2.4: Code Generation (Orchestrator-Owned)

1. Read reference file sections per agent-loading-guide.md Phase 1 + Phase 2 maps
2. Generate mesh script from Phase 1 validated plan:
   - Import STEP geometry via `gmsh.model.occ.importShapes()` (with error handling)
     OR construct parametric geometry via gmsh built-in/OCC kernel
   - Call `gmsh.model.occ.synchronize()` after ALL OCC operations
   - Define physical groups (inlet, outlet, walls, membrane, cell_region)
   - Apply mesh refinement near walls and membrane interfaces
   - Estimate memory requirements before meshing
   - Check mesh quality after generation
   - Convert to DOLFINx via `dolfinx.io.gmsh.model_to_mesh()`
3. Generate flow solver script from Phase 2 validated plan:
   - Define Taylor-Hood P2/P1 function space
   - Apply boundary conditions from physical groups
   - For Re < 1: assemble Stokes variational form (linear solve)
   - For Re >= 1: Navier-Stokes with Newton continuation:
     a. Solve Stokes for initial guess
     b. Ramp Re through intermediate values if Re > 10
     c. Use previous solution as initial guess for each step
   - Select solver: MUMPS for < 50K DOFs, GMRES+ILU otherwise
   - Include solver progress monitor for 3D
   - Save velocity + pressure fields to XDMF checkpoint
4. Apply ALL code generation rules from Section 13 (Code Generation Protocol)
5. Run syntax check: `python3 -c "import ast; ast.parse(open('script.py').read())"`

### Step 2.5: Execution and Validation

1. Execute mesh + flow scripts (Bash for 2D, Write for 3D)
2. Validate:
   - Solver converged (residual < tolerance)
   - Mass conservation: |integral u.n ds_inlet + integral u.n ds_outlet| < 1e-6
   - No NaN or Inf in solution arrays
   - No unphysical flow patterns (reversed flow at inlet, etc.)
   - If Poiseuille comparison available: L2 error < 1% (P2 elements give ~1e-12 for
     quadratic solutions; see `references/validation-benchmarks.md` Benchmark 1)

### Step 2.5b: Self-Correction Loop (If Execution or Validation Fails)

1. **Collect** error output (traceback, log messages, solver convergence history)
2. **Classify** error type: `import_error`, `mesh_error`, `solver_divergence`,
   `assertion_failure`, `numerical_instability`
3. **Simple errors** (import path, syntax, assertion): Fix directly and retry
   (max 2 direct retries)
4. **Complex errors** (solver divergence, numerical instability): Invoke `cfd-reviewer`
   in error diagnosis mode:
   - Pass: error output, error_history from previous attempts, mathematical specification
   - Reviewer classifies root cause and recommends specific fix
   - Check error_history to avoid recommending previously-failed fixes
5. **Regenerate code** with reviewer feedback applied
6. **Track** error_history (append each attempt: error_type, error_message, fix_attempted,
   fix_outcome)
7. **After 3 total retries**: Escalate to user with full error summary:
   ```
   EXECUTION FAILED after 3 attempts.
   Error type: [classification]
   Attempts:
   1. [fix attempted] -> [outcome]
   2. [fix attempted] -> [outcome]
   3. [fix attempted] -> [outcome]
   Recommendation: [suggested next step]
   ```

**Quality Gate 2b**: Flow validation passes (mass conservation, convergence, no NaN/Inf).

**Handoff**: Write Phase 2 flow-result handoff to `{session_dir}/handoffs/phase2-flow-result.yaml`

---

## 10. Phase 3 -- Transport Planning, Code Generation, and Execution

Mirrors Phase 2 structure with transport-specific content.

### Step 3.1: Swarm Discussion (FULL Mode Only)

Same pattern as Phase 1-2 with Swarm 3 template from `references/swarm-framing-templates.md`.

### Step 3.2: Mathematical Analysis

Invoke `cfd-mathematician` via Task tool:
- **Input**: Transport parameters (Pe, Da, species properties), velocity field summary
  from Phase 2, swarm synthesis (if FULL mode)
- **Task**: Write transport variational form with SUPG, analyze stabilization parameter
  (Pe-dependent xi formula), verify regularized MM convergence properties, estimate
  transport convergence order
- **Timeout**: 8 minutes

### Step 3.3: Adversarial Review

Invoke `cfd-reviewer` via Task tool:
- **Input**: Transport plan from Steps 3.1-3.2, mathematician handoff output
- **Task**: Challenge SUPG parameter for overflow risk (Pe > 710 with cosh/sinh), check
  MM regularization epsilon adequacy, verify Robin BC formulation for membrane, check
  Newton solver convergence expectations
- **Timeout**: 8 minutes

### Quality Gate 3a: Transport Plan Approved

Same pattern as Quality Gates 1 and 2a:
- Mechanical checks + reviewer approval
- Max 1 retry on rejection, then escalate to user

### Step 3.4: Code Generation (Orchestrator-Owned)

1. Read reference file sections per agent-loading-guide.md Phase 3 map
2. Generate transport script from Phase 3 validated plan:
   - Define P1 function space for concentration
   - Estimate Peclet number: Pe = |u| * h / (2D)
   - Implement SUPG stabilization (see Section 13 -- ALWAYS rules)
   - Implement regularized Michaelis-Menten (see Section 13 -- ALWAYS rules)
   - Implement membrane permeation as Robin BC (Fick's law)
   - Use Newton solver (nonlinear due to Michaelis-Menten)
   - Import NonlinearProblem from `dolfinx.fem.petsc`, NewtonSolver from `dolfinx.nls.petsc`
   - Include post-solve checks: negative concentration warning, species conservation
3. Apply ALL code generation rules from Section 13
4. Run syntax check

### Step 3.5: Execution and Validation

1. Execute transport script (Bash for 2D, Write for 3D)
2. Validate:
   - Newton solver converged
   - Species conservation: |inlet_flux + outlet_flux + membrane_flux + total_reaction| < 1%
     of total reaction
   - min(c) check: warn if significantly negative (suggests insufficient mesh resolution)

### Step 3.5b: Self-Correction Loop

Same pattern as Step 2.5b. Invoke cfd-reviewer in error diagnosis mode for complex errors.
Max 3 total retries, then escalate to user.

**Quality Gate 3b**: Transport validation passes (species conservation, Newton convergence,
no large negative concentration regions).

**Handoff**: Write Phase 3 transport-result handoff to
`{session_dir}/handoffs/phase3-transport-result.yaml`

---

## 11. Phase 4 -- Post-Processing and Visualization

**Orchestrator-owned** (no agent involvement).

### Step 4.1: Visualization Scripts

1. Read `references/fenicsx-patterns.md` Section 12 (PyVista patterns)
2. Generate visualization script with:
   - Velocity magnitude contour plot
   - O2 concentration field with depletion zones highlighted
   - Streamlines from inlet (3D)
   - Cross-section slices at multiple positions (3D)
   - Centerline line probe (O2 vs. position)
   - Screenshots saved to PNG files
   - VTK files saved for ParaView
   - Headless fallback: `if not pyvista.system_supports_plotting(): pyvista.OFF_SCREEN = True`
3. Execute 2D visualizations; save 3D scripts for manual execution

### Step 4.2: Mesh Convergence Study (Optional -- for Publication Quality)

1. Read `references/validation-benchmarks.md` Section 4 (convergence protocol)
2. Generate convergence study script:
   - 4 mesh refinement levels: h, h/2, h/4, h/8 (or 3 for 3D)
   - Solve on each level
   - Compute QoI: outlet average O2, max velocity, min O2 in cell region
   - Richardson extrapolation for grid-independent value
   - Convergence rate: p = log(e_{i-1}/e_i) / log(2)
   - Non-monotone convergence detection and warning

**Quality Gate 4a** (convergence study):
- Monotone convergence (error decreases with refinement)
- Convergence rate matches expected order (O(h^2) for P1 concentration, O(h^3) for P2 velocity)
- QoI changes < 1% between last two levels

---

## 12. Phase 5 -- Summary and Deliverables

1. Collect all handoff documents from `{session_dir}/handoffs/`
2. Report to user:
   - What was simulated (geometry, physics, parameters)
   - Key agent decisions (mathematician recommendations, reviewer concerns, resolutions)
   - Validation results (mass conservation, species conservation, convergence)
   - Files produced (scripts, XDMF checkpoints, VTK files, PNG images)
3. Offer to preserve or clean up session directory
4. Update state file with `status: "complete"`

---

## 13. Code Generation Protocol

When generating FEniCSx Python scripts, ALWAYS follow these rules. These rules are
embedded inline because they guard against code-level errors that produce silently
wrong results. They must always be in the orchestrator's context.

### ALWAYS Include

1. **Version assertion header** at the top of every script:
   ```python
   import importlib.metadata as _meta
   _ver = tuple(int(x) for x in _meta.version("fenics-dolfinx").split(".")[:2])
   assert _ver >= (0, 10), f"Requires FEniCSx >= 0.10. See references/environment-setup.md."
   ```

2. **Reproducibility metadata header**: date, dolfinx version, basix version, gmsh version,
   numpy version, Python version, OS, all physical parameters

3. **Complete, self-contained scripts**: Every script must run standalone from the command
   line. Never generate code fragments that require manual assembly.

4. **Inline comments** explaining each step: what it does, why it is needed, what the
   expected output is

5. **Regularized Michaelis-Menten**: Always use `c_pos = (c + sqrt(c^2 + eps^2)) / 2`
   with `eps = 1e-10 * c_inlet`. Never use the raw `c / (Km + c)` form.

6. **SUPG stabilization** for all species transport (check Pe first, but enable by default):
   Use the numerically stable form `xi = conditional(gt(Pe, 1.0), 1.0 - 1.0/Pe, Pe/3.0)`.
   Never use cosh/sinh (overflows for Pe > 710).

7. **Post-solve quality checks**: mass conservation, species conservation, negative
   concentration warning, NaN/Inf detection

8. **Save results** to XDMF (for checkpointing) and VTK/PVD (for PyVista/ParaView):
   ```python
   with io.XDMFFile(comm, "results.xdmf", "w") as xdmf:
       xdmf.write_mesh(domain)
       xdmf.write_function(uh, 0.0)
   ```

### NEVER Do

- Generate code fragments (always complete scripts)
- Use `from dolfin import *` (legacy FEniCS, not FEniCSx)
- Use `dolfinx.io.gmshio` (v0.10 renamed to `dolfinx.io.gmsh`)
- Use raw Michaelis-Menten `c / (Km + c)` without regularization
- Use `cosh(Pe)/sinh(Pe)` for SUPG parameter (numerical overflow)
- Skip physical group definition before meshing
- Forget `gmsh.model.occ.synchronize()` after OCC operations
- Use equal-order elements (P1/P1) for flow without stabilization
- Omit the pressure reference point (at least one pressure BC needed)

### Pattern Library

Use code patterns from `references/fenicsx-patterns.md` as the primary source for all
generated code. These patterns have been validated against FEniCSx v0.10 API.

Use parameter values from `references/physics-models.md` lookup tables for bioreactor-specific
constants (fluid properties, diffusion coefficients, Michaelis-Menten kinetics by cell type).

### Additional Code Generation Rules

- **Taylor-Hood P2/P1** for Stokes/Navier-Stokes flow (standard inf-sup stable pair)
- **Newton continuation** for Re > 10: Solve Stokes for initial guess, ramp Re through
  intermediate values [1, 10, 50, target], use previous solution as initial guess
- **MUMPS** for < 50K DOFs (direct solver), **GMRES+ILU** for larger problems (iterative)
- **NonlinearProblem** from `dolfinx.fem.petsc`, **NewtonSolver** from `dolfinx.nls.petsc`
- **eps = 1e-10 * c_inlet** for Michaelis-Menten regularization parameter

---

## 14. Error Handling Protocol

### Simulation Error Handling

When errors occur during simulation, follow this diagnostic sequence:

| Error Stage | Symptom | First Action | Reference |
|---|---|---|---|
| Stage 0: Environment | ImportError, ModuleNotFoundError | Run pre-flight check | troubleshooting-guide.md Stage 0 |
| Stage 1: Geometry/Mesh | OCC exception, 0 volumes | Check STEP file, try defeaturing | troubleshooting-guide.md Stage 1 |
| Stage 2: Physics Setup | TypeError in BCs, shape mismatch | Check function space indices | troubleshooting-guide.md Stage 2-3 |
| Stage 3: Flow Solver | Newton divergence | Stokes initial guess -> ramp Re -> reduce relaxation | troubleshooting-guide.md Stage 4 |
| Stage 3: Flow Solver | MUMPS pivot error | Switch to iterative solver | troubleshooting-guide.md Stage 4 |
| Stage 3: Flow Solver | PETSc OOM | Reduce mesh, use iterative solver | troubleshooting-guide.md Stage 4 |
| Stage 4: Transport | Negative concentrations | Check SUPG active, refine mesh near sinks | troubleshooting-guide.md Stage 4 |
| Stage 5: Validation | Conservation violated > 1% | Refine mesh, check BCs, check form | troubleshooting-guide.md Stage 5 |
| Stage 6: Visualization | Empty plot, PyVista crash | Check data arrays, try headless mode | troubleshooting-guide.md Stage 5-6 |

### Solver Divergence Recovery Sequence

If the Newton solver diverges, try these fixes in order:

1. **Use Stokes solution as initial guess** (most common fix)
2. **Ramp Reynolds number** through intermediate values: [1, 10, 50, target]
3. **Reduce relaxation factor** to 0.5 (then 0.3, then 0.1)
4. **Switch to Picard iteration** instead of Newton
5. **Refine mesh** near boundaries and high-gradient regions
6. If all fail: suggest the user simplify the geometry or reduce Re

Always provide actionable diagnostic messages. Never just report an error code.

### Orchestration Error Handling

| Failure Mode | Detection | Response | Escalation |
|---|---|---|---|
| Agent timeout | Task tool timeout reached | Retry once with simplified prompt | Proceed with degraded mode |
| Malformed handoff | YAML parse failure or missing required fields | Retry agent with explicit format feedback | Proceed without that agent's input |
| Swarm confidence < 3/5 | `confidence_score` field in synthesis | Flag to user as informational | User decides whether to proceed |
| Context approaching limit | Phase number > 2 in FULL mode | Summarize previous handoffs to 200-token summaries | Skip remaining swarm invocations |
| Session directory missing | File check before each phase | Attempt recreation from state YAML | User informed; may need restart |
| Agent skill not found | SKILL.md existence check in Phase 0 | Degrade gracefully per Section 1 | User told to run sync-config.py |

---

## 15. Timeout Configuration

### Per-Agent Timeouts

| Agent | Timeout | Fallback on Timeout |
|---|---|---|
| brainstorming-pm | 15 min | Skip swarm; proceed with mathematician + reviewer only |
| cfd-mathematician | 8 min | Retry once with simplified prompt; if still fails, use reference file defaults |
| cfd-reviewer | 8 min | Retry once simplified; if still fails, APPROVED_WITH_WARNINGS + log "NO ADVERSARIAL REVIEW" |

### Per-Phase Execution Timeouts

| Problem Type | Expected Runtime | Execution Strategy |
|---|---|---|
| 2D validation (Tier 1) | < 2 minutes | Bash (timeout=300000) |
| 2D coupled (Tier 2) | 2-10 minutes | Bash (timeout=600000) |
| 3D coarse validation | 5-15 minutes | Bash (timeout=600000) if small; Write+manual otherwise |
| 3D production (Tier 3) | 30 min - 4 hours | Write script only. User runs manually. |
| 3D with convergence (Tier 4) | 4-24 hours | Write script only. User runs with MPI. |
| Mesh convergence study | 4x single-run time | Write script only for 3D. |

**Rule**: If estimated runtime exceeds 10 minutes, generate the script and instruct the
user to run it manually. Do not attempt to execute long-running simulations via Bash.

**For 3D simulations**: Always include MPI instructions:
```bash
mpirun -np 4 python simulation.py
```

---

## 16. Quality Gate Escalation Protocol

When the mathematician and reviewer cannot agree after 1 retry, present the disagreement
to the user in this structured format:

```
AGENT DISAGREEMENT: [Phase Name]

Mathematician recommends: [1-2 sentence summary]
Rationale: [1 sentence]

Reviewer objects: [1-2 sentence summary]
Concern: [1 sentence]

Options:
(A) Proceed with mathematician's recommendation (risk: [reviewer's concern])
(B) Proceed with reviewer's preferred approach (trade-off: [mathematician's concern])
(C) Use conservative approach: [orchestrator's auto-generated compromise]
(D) Provide your own parameters

Recommended: (C) [brief justification]
```

**Conservative compromise rule**: Always take the more cautious choice from both agents --
finer mesh, more stable solver, lower tolerance, more regularization.

**Escalation sequence**:
1. First rejection: Reviewer's `blocking_issues` become HARD CONSTRAINTS for mathematician retry
2. After 1 retry still REJECTED: Escalate to user with the structured format above
3. User selects option (A/B/C/D) or provides custom parameters
4. Orchestrator proceeds with user's choice

---

## 17. Orchestrator Context Management

The orchestrator loads reference file sections on a per-phase basis to control context
window consumption. Agent discussion details stay in handoff files on disk -- they are
NOT loaded back into the orchestrator context.

| Phase | Load Into Context | Summarize from Previous | Do NOT Load |
|---|---|---|---|
| 0 | environment-setup.md Section 4 | None | All other references |
| 1 | Agent handoffs (current phase) | Phase 0 summary | Full reference files |
| 2 | Phase 1 mesh-plan (full) + Phase 2 handoffs | Phase 0 summary | Phase 1 agent discussions |
| 3 | Phase 2 flow-result (full) + Phase 3 handoffs | Phases 0-1 summaries | Phase 2 agent discussions |
| 4 | Phase 3 transport-result | Phases 0-2 summaries | All agent discussions |
| 5 | Final handoffs only | All phases summarized | All discussions |

**Previous phase summaries**: Max 200 tokens per phase. Include only: key decisions,
approval status, validation metrics. Not full discussion content.

**Section 13 NEVER rules are always in context** (they are embedded inline in this
SKILL.md, not in a reference file).

---

## 18. Phase Dependencies

| Phase | Produces | Consumed By | Invalidation Rule |
|---|---|---|---|
| 0 | Session config, mode selection | All phases | Re-run Phase 0 if environment changes |
| 1 | mesh-plan handoff | Phase 2 (mesh code gen) | If Phase 1 re-runs, invalidate Phases 2-4 |
| 2 | flow-result handoff | Phase 3 (velocity field) | If Phase 2 re-runs, invalidate Phases 3-4 |
| 3 | transport-result handoff | Phase 4 (visualization) | If Phase 3 re-runs, invalidate Phase 4 |
| 4 | Visualizations, convergence study | Phase 5 (summary) | Can re-run independently |
| 5 | User summary | None | Can re-run independently |

**Rule**: If Phase N is re-executed, all phases > N are invalidated and must be re-run.

---

## 19. Complexity Tiers

### Tier 1: Validation (5-15 minutes) --> DIRECT Mode

2D Poiseuille flow. Tests the entire pipeline without any user-specific physics.
**Always start here** when the user has never run a FEniCSx simulation before.

See: `examples/01-2d-channel-flow.md`

### Tier 2: Single Physics (15-30 minutes) --> LITE Mode

2D user geometry + O2 transport with Michaelis-Menten and membrane BCs.
Demonstrates the full two-phase workflow on a simple geometry.

See: `examples/02-2d-oxygen-transport.md`

### Tier 3: Coupled Multiphysics (1-4 hours) --> FULL Mode

2D or 3D STEP geometry + full coupled physics. Production-quality simulation with
iterative solvers, Newton continuation, and validation suite.

See: `examples/03-3d-cartridge-template.md`

### Tier 4: Production 3D (4-24 hours) --> FULL Mode

Full 3D cartridge with mesh convergence study, publication-quality visualization,
and complete reproducibility metadata. Always run with MPI on a workstation.

**Progression rule**: Suggest the user start at Tier 1 and progress upward. Do not
jump to Tier 3 or 4 without first validating the environment at Tier 1.

**Mode mapping**:
| Tier | Default Mode | Rationale |
|---|---|---|
| 1 | DIRECT | Well-understood benchmark; agents add no value |
| 2 | LITE | Mathematician + reviewer sufficient for 2D |
| 3 | FULL | 3D complexity benefits from swarm discussion |
| 4 | FULL | Production quality requires thorough review |

---

## 20. Session State File Schema

State file location: `/tmp/cfd-bioreactor-state-{session-id}.yaml`

```yaml
session:
  session_id: "20260221-143000-12345"
  session_dir: "/tmp/cfd-bioreactor-session-20260221-143000-12345"
  mode: "LITE"                    # DIRECT | LITE | FULL
  tier: 2                         # 1-4
  current_phase: 2                # 0-5
  status: "running"               # running | paused | complete | failed

phases_completed: [0, 1]          # List of completed phase numbers
phases_valid: [0, 1]              # List of phases whose results are still valid

phase_results:
  phase_0:
    status: "complete"
    handoff_path: null
    timestamp: "2026-02-21T14:30:00Z"
  phase_1:
    status: "complete"
    handoff_path: "{session_dir}/handoffs/phase1-mesh-plan.yaml"
    reviewer_status: "APPROVED_WITH_WARNINGS"
    retries_used: 0
    timestamp: "2026-02-21T14:35:00Z"
  phase_2:
    status: "running"
    handoff_path: null
    reviewer_status: null
    retries_used: 0
    execution_retries_used: 0     # For self-correction loop

timestamps:
  started: "2026-02-21T14:30:00Z"
  last_updated: "2026-02-21T14:40:00Z"

errors: []                        # List of non-fatal errors/warnings encountered
```

---

## 21. Session Resume Protocol

On startup, check for existing state files:

1. Search for `/tmp/cfd-bioreactor-state-*.yaml` with `status != "complete"`
2. If found: display last completed phase, offer to resume
   ```
   Found incomplete session from [timestamp].
   Last completed phase: Phase [N] ([phase_name])
   Mode: [DIRECT/LITE/FULL], Tier: [N]
   Resume from Phase [N+1]? (yes/no/restart)
   ```
3. If user chooses resume:
   - **Phase 0**: Always re-run (fast, catches environment changes)
   - **Phases 1-4**: Resumable if handoff files exist at recorded paths
   - Verify handoff file existence before resuming
4. If user chooses restart: Create new session, archive old state file

---

## 22. Integration with Other Skills

| Need | Skill | How to Use |
|---|---|---|
| Quick feasibility estimate | `calculator` | Before full CFD: estimate Re, Pe, Da, O2 depletion depth. A 5-minute calculation can determine if CFD is even needed. |
| Literature parameter values | `bioinformatician` or `researcher` | Find Vmax, Km for specific cell types; diffusion coefficients in specific media. |
| Debug generated notebook | `notebook-debugger` | If user is running in Jupyter and encounters FEniCSx import/kernel issues. |
| Cross-check CFD results | `calculator` | After simulation: verify CFD results against analytical estimates. If they disagree by > 10x, investigate. |
| Mathematical analysis | `cfd-mathematician` | Via Task tool. Variational formulations, stability analysis, convergence estimates. |
| Adversarial review | `cfd-reviewer` | Via Task tool. Engineering review with severity-rated challenges and approval status. |
| Multi-perspective brainstorming | `brainstorming-pm` | Via Task tool. FULL mode only. Swarm discussion at 3 decision points. |

---

## 23. Reference Files

This skill includes the following reference files. Read them as needed during workflow
execution -- do not load all at once. See `references/agent-loading-guide.md` for
section-level loading maps per agent.

| File | When to Read | Purpose |
|---|---|---|
| `references/environment-setup.md` | Phase 0 (pre-flight) | Installation, pre-flight script, degradation modes |
| `references/physics-models.md` | Phases 2, 3 (code gen) | Equations, variational forms, parameter tables |
| `references/mesh-generation-guide.md` | Phase 1 (code gen) | STEP import, physical groups, refinement, quality |
| `references/fenicsx-patterns.md` | Phases 1-4 (all code gen) | FEniCSx v0.10 API patterns (the code library) |
| `references/validation-benchmarks.md` | Phases 2, 3 (validation) | Analytical solutions, convergence protocols |
| `references/troubleshooting-guide.md` | On error | Error catalog by stage, diagnostic commands |
| `references/orchestrator-handoff-schema.md` | Agent invocations | YAML handoff contracts for all agent communication |
| `references/agent-loading-guide.md` | Agent invocations | Maps agents to reference file sections to load |
| `references/swarm-framing-templates.md` | FULL mode swarms | Pre-written challenge templates for brainstorming-pm |
| `examples/environment.yml` | Phase 0 (if not installed) | Conda environment specification |

---

## 24. Notes

- **Version pinning**: All code patterns target FEniCSx v0.10. If a newer version is
  released, the version assertion will catch the mismatch. Update patterns before using
  with a new API version.

- **SI units everywhere**: All parameters are in SI (kg, m, s, mol, Pa) unless the user
  explicitly specifies otherwise. If a STEP file appears to be in mm (bounding box > 1 m),
  warn the user about potential unit mismatch.

- **Serial execution by default**: Generated scripts run in serial (single core). For 3D
  production runs, add MPI instructions. The scripts are MPI-aware (use `MPI.COMM_WORLD`
  throughout) and work correctly with `mpirun -np N`.

- **Reproducibility**: Every generated script includes a metadata header with all version
  numbers, parameters, and settings. Results can be reproduced by re-running the same
  script in the same environment.

- **This skill generates code but does not modify system packages**: It writes Python
  scripts to the user's working directory. It never installs packages, modifies conda
  environments, or changes system configuration.

- **Physical group convention**: inlet=1, outlet=2, walls=3, membrane=4, cell_region=5,
  fluid_volume=10. Consistent across all examples and generated scripts.

- **gmsh synchronize**: After ANY `gmsh.model.occ.*` operation (addBox, addCylinder,
  importShapes, etc.), ALWAYS call `gmsh.model.occ.synchronize()` before accessing
  entities, assigning physical groups, or meshing. Forgetting this is the single most
  common gmsh error.

- **v2.0 multi-agent orchestrator**: This is v2.0 of the cfd-bioreactor skill. v1.0 was
  a single-context generate-explain-validate skill; v2.0 is a multi-agent orchestrator
  with cfd-mathematician, cfd-reviewer, and brainstorming-pm coordination.
