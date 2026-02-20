---
name: cfd-bioreactor
version: "1.0"
description: >
  Use when user needs to simulate fluid flow through bioprocess cartridges,
  bioreactor geometries, or membrane devices using FEniCSx. Triggers on
  requests for CFD simulation, Navier-Stokes solving, species transport
  modeling, oxygen consumption modeling, or mesh generation from CAD files.
prerequisites:
  - FEniCSx v0.10+ installed (see references/environment-setup.md)
  - gmsh with OpenCASCADE kernel
  - PyVista for visualization
  - Conda environment recommended (see examples/environment.yml)
success_criteria:
  - Generated scripts run end-to-end without errors
  - Poiseuille flow validation achieves < 1% L2 error
  - Conservation checks pass (mass < 1e-6, species < 1%)
  - All visualizations render correctly
  - Scripts include complete reproducibility metadata
estimated_duration: 30min (2D validation) to 4hr+ (3D production)
allowed_tools: [Read, Write, Bash, Skill]
metadata:
  skill-author: David Angeles Albores
  category: scientific-simulation
  workflow: generate-explain-validate
  integrates-with: [calculator, bioinformatician]
---

# CFD Bioreactor Simulation Skill

Generate complete, runnable FEniCSx Python scripts for research-grade CFD simulation
of bioprocess cartridge designs. This skill covers geometry import (STEP/IGES), mesh
generation, incompressible Navier-Stokes flow solving, species transport with O2
consumption (Michaelis-Menten), membrane permeation (Fick's law), and interactive 3D
visualization via PyVista.

All generated code targets FEniCSx v0.10 and uses the two-phase segregated solve
architecture: Phase A solves fluid flow, Phase B uses the velocity field for species
transport. Each phase is independently validatable.

---

## 1. When to Use This Skill

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

## 2. When NOT to Use This Skill

Do NOT use this skill for:

- **Compressible flow or Mach effects**: Not relevant for bioreactors (Ma << 0.01)
- **Turbulence modeling**: Laminar flow assumed (Re < 100). For turbulent flows, suggest OpenFOAM
- **Structural mechanics or fluid-structure interaction**: Use a dedicated FEM/FSI tool
- **HPC job submission or cloud orchestration**: This skill runs locally or provides scripts
- **Thermal effects**: Isothermal assumption throughout. Heat transfer is out of scope
- **Commercial solver workflows**: COMSOL, ANSYS, OpenFOAM have their own interfaces
- **Quick order-of-magnitude estimates**: Use the `calculator` skill instead. A back-of-envelope
  calculation takes 5 minutes; a full CFD simulation takes 30 minutes to hours
- **Transient simulations with adaptive time-stepping**: Steady-state only in v1.0

**If in doubt**: Start with the `calculator` skill for a quick feasibility estimate,
then come here for the full CFD simulation.

---

## 3. Pre-Flight Validation Protocol (Step 0)

**MANDATORY**: Before executing ANY workflow step, run pre-flight validation.

### Action

1. Read `references/environment-setup.md` Section 4 (Pre-Flight Validation Script)
2. Generate `preflight_check.py` and execute it via Bash
3. Interpret results

### Decision Logic

| Pre-Flight Result | Action |
|---|---|
| All PASS | Proceed to workflow |
| FAIL on dolfinx | Provide conda/Docker install instructions from environment-setup.md. STOP. |
| FAIL on gmsh OCC | Can still do parametric geometry. Warn user that STEP import is disabled. |
| WARN on PyVista | Proceed; export to VTK instead of interactive plots. |
| WARN on MUMPS | Proceed; use iterative solver (GMRES+ILU) instead of direct solver. |
| FAIL on petsc4py | Cannot solve. Provide install instructions. STOP. |
| WARN on RAM < 8 GB | Proceed for 2D only. Warn that 3D will require coarse meshes. |

### Critical Rule

If pre-flight fails on dolfinx or petsc4py, do NOT attempt to generate simulation
scripts. The scripts will not run. Instead, help the user set up the environment first.

---

## 4. Workflow: Two-Phase Segregated Architecture

The simulation follows a two-phase segregated solve: Phase A solves fluid flow
independently, Phase B uses the frozen velocity field for species transport. This
architecture is standard for dilute species in bioreactors where concentration does
not affect fluid properties.

### Phase A: Fluid Flow (Steps 1-3)

#### Step 1: Geometry Import and Mesh Generation

**Tool**: Write (generate mesh script) + Bash (execute)

**Input**: STEP file path OR parametric dimensions from user

**Actions**:
1. Read `references/mesh-generation-guide.md` for meshing patterns
2. Generate a mesh script that:
   - Imports STEP geometry via `gmsh.model.occ.importShapes()` (with full error handling)
   - OR constructs parametric geometry via gmsh built-in/OCC kernel
   - Calls `gmsh.model.occ.synchronize()` after ALL OCC operations
   - Defines physical groups (inlet, outlet, walls, membrane, cell_region)
   - Applies mesh refinement near walls and membrane interfaces
   - Estimates memory requirements before meshing
   - Checks mesh quality after generation
   - Converts to DOLFINx via `dolfinx.io.gmsh.model_to_mesh()`
3. Execute the script (2D) or save for manual execution (3D)

**Quality Gate**:
- Mesh quality (min scaled Jacobian) > 0.1
- Element count within memory budget (check `references/mesh-generation-guide.md` Section 5)
- All physical groups assigned and non-empty
- Memory estimate < 70% available RAM

**Error handling**: See `references/troubleshooting-guide.md` Stage 1.

#### Step 2: Flow Solver Setup and Execution

**Tool**: Write (generate flow script) + Bash (execute for 2D; Write only for 3D)

**Input**: Mesh from Step 1, fluid properties from user

**Actions**:
1. Read `references/fenicsx-patterns.md` for solver patterns
2. Read `references/physics-models.md` for variational formulations
3. Generate a flow script that:
   - Defines Taylor-Hood P2/P1 function space (see patterns Section 3)
   - Applies boundary conditions from physical groups (see patterns Section 4)
   - For Re < 1: assembles Stokes variational form (linear solve)
   - For Re >= 1: Navier-Stokes with Newton continuation (see patterns Section 6):
     a. Solve Stokes for initial guess
     b. Ramp Re through intermediate values if Re > 10
     c. Use previous solution as initial guess for each step
   - Selects solver: MUMPS for < 50K DOFs, GMRES+ILU otherwise (see patterns Section 8)
   - Includes solver progress monitor for 3D (see patterns Section 9)
   - Saves velocity + pressure fields to XDMF checkpoint

**Quality Gate**:
- Solver converged (residual < tolerance)
- Mass conservation: |integral u.n ds_inlet + integral u.n ds_outlet| < 1e-6
- No NaN or Inf in solution arrays

**Error handling**: See `references/troubleshooting-guide.md` Stage 4.

#### Step 3: Flow Validation

**Tool**: Bash (run validation) + Read (check results)

**Actions**:
1. Check mass conservation (integral of u.n over all boundaries)
2. If geometry is a simple channel: compare with Poiseuille analytical solution
3. Visualize velocity field (PyVista contour plot)
4. Report validation results

**Quality Gate**:
- Conservation check passes
- No unphysical flow patterns (reversed flow at inlet, etc.)
- If Poiseuille comparison available: L2 error < 1% (note: P2 elements give ~1e-12 for
  quadratic solutions; see `references/validation-benchmarks.md` Benchmark 1)

**Decision Point**: If flow validation fails, diagnose and fix BEFORE proceeding to
Phase B. Common fixes: refine mesh, check BCs, adjust solver parameters.

### Phase B: Species Transport (Steps 4-7)

Phase B uses the velocity field from Phase A. The flow field is frozen -- there is no
back-coupling from concentration to flow (valid for dilute species).

#### Step 4: Species Transport Setup and Execution

**Tool**: Write (generate transport script) + Bash (execute for 2D)

**Input**: Velocity field from Phase A, species properties, sink parameters from user

**Actions**:
1. Read `references/fenicsx-patterns.md` Section 7 (SUPG) and Section 10 (quality checks)
2. Read `references/physics-models.md` for reaction models and parameter tables
3. Generate a transport script that:
   - Defines P1 function space for concentration
   - Estimates Peclet number: Pe = |u| * h / (2D)
   - Implements SUPG stabilization (MANDATORY for Pe > 1, safe for all Pe):
     ```python
     xi = conditional(gt(Pe, 1.0), 1.0 - 1.0/Pe, Pe/3.0)
     tau = h / (2.0 * u_mag) * xi
     v_supg = v + tau * dot(u, grad(v))
     ```
     IMPORTANT: Do NOT use the cosh/sinh formula (overflows for Pe > 710).
     Always use the numerically stable simplified form shown above.
   - Implements regularized Michaelis-Menten (ALWAYS regularized, NEVER raw form):
     ```python
     c_pos = (c + sqrt(c**2 + eps**2)) / 2.0
     R = -Vmax * c_pos / (Km + c_pos)
     ```
     with eps = 1e-10 * c_inlet
   - Implements membrane permeation as Robin BC (Fick's law)
   - Uses Newton solver (nonlinear due to Michaelis-Menten)
   - Import NonlinearProblem from `dolfinx.fem.petsc`, NewtonSolver from `dolfinx.nls.petsc`
   - Includes post-solve checks: negative concentration warning, species conservation

**Quality Gate**:
- Newton solver converged
- Species conservation: |inlet_flux + outlet_flux + membrane_flux + total_reaction| < 1%
  of total reaction
- min(c) check: warn if significantly negative (suggests insufficient mesh resolution)

#### Step 5: Transport Validation

**Tool**: Bash (run checks) + Read (results)

**Actions**:
1. Species conservation balance (all fluxes + reaction = 0)
2. Negative concentration report
3. If applicable: comparison with 1D analytical solution at centerline
   (see `references/validation-benchmarks.md` Benchmark 2)

**Quality Gate**:
- Conservation < 1% imbalance
- No large negative concentration regions

#### Step 6: Post-Processing and Visualization

**Tool**: Write (generate viz script) + Bash (execute)

**Actions**:
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

**Output**: Interactive visualizations + saved images + VTK files

#### Step 7: Mesh Convergence Study (Optional -- for Publication Quality)

**Tool**: Write (generate convergence script) + Bash (execute)

**Actions**:
1. Read `references/validation-benchmarks.md` Section 4 (convergence protocol)
2. Generate convergence study script:
   - 4 mesh refinement levels: h, h/2, h/4, h/8 (or 3 for 3D)
   - Solve on each level
   - Compute QoI: outlet average O2, max velocity, min O2 in cell region
   - Richardson extrapolation for grid-independent value
   - Convergence rate: p = log(e_{i-1}/e_i) / log(2)
   - Non-monotone convergence detection and warning

**Quality Gate**:
- Monotone convergence (error decreases with refinement)
- Convergence rate matches expected order (O(h^2) for P1 concentration, O(h^3) for P2 velocity)
- QoI changes < 1% between last two levels

---

## 5. Tool Selection Table

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

---

## 6. Code Generation Protocol

When generating FEniCSx Python scripts, ALWAYS follow these rules:

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

---

## 7. Error Handling Protocol

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

---

## 8. Timeout Configuration

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

## 9. Complexity Tiers

### Tier 1: Validation (5-15 minutes)

2D Poiseuille flow. Tests the entire pipeline without any user-specific physics.
**Always start here** when the user has never run a FEniCSx simulation before.

See: `examples/01-2d-channel-flow.md`

### Tier 2: Single Physics (15-30 minutes)

2D user geometry + O2 transport with Michaelis-Menten and membrane BCs.
Demonstrates the full two-phase workflow on a simple geometry.

See: `examples/02-2d-oxygen-transport.md`

### Tier 3: Coupled Multiphysics (1-4 hours)

2D or 3D STEP geometry + full coupled physics. Production-quality simulation with
iterative solvers, Newton continuation, and validation suite.

See: `examples/03-3d-cartridge-template.md`

### Tier 4: Production 3D (4-24 hours)

Full 3D cartridge with mesh convergence study, publication-quality visualization,
and complete reproducibility metadata. Always run with MPI on a workstation.

**Progression rule**: Suggest the user start at Tier 1 and progress upward. Do not
jump to Tier 3 or 4 without first validating the environment at Tier 1.

---

## 10. Integration with Other Skills

| Need | Skill | How to Use |
|---|---|---|
| Quick feasibility estimate | `calculator` | Before full CFD: estimate Re, Pe, Da, O2 depletion depth. A 5-minute calculation can determine if CFD is even needed. |
| Literature parameter values | `bioinformatician` or `researcher` | Find Vmax, Km for specific cell types; diffusion coefficients in specific media. |
| Debug generated notebook | `notebook-debugger` | If user is running in Jupyter and encounters FEniCSx import/kernel issues. |
| Cross-check CFD results | `calculator` | After simulation: verify CFD results against analytical estimates. If they disagree by > 10x, investigate. |

---

## 11. Reference Files

This skill includes the following reference files. Read them as needed during workflow
execution -- do not load all at once.

| File | When to Read | Purpose |
|---|---|---|
| `references/environment-setup.md` | Step 0 (pre-flight) | Installation, pre-flight script, degradation modes |
| `references/physics-models.md` | Steps 2, 4 (solver setup) | Equations, variational forms, parameter tables |
| `references/mesh-generation-guide.md` | Step 1 (geometry/mesh) | STEP import, physical groups, refinement, quality |
| `references/fenicsx-patterns.md` | Steps 1-6 (all code gen) | FEniCSx v0.10 API patterns (the code library) |
| `references/validation-benchmarks.md` | Steps 3, 5, 7 (validation) | Analytical solutions, convergence protocols |
| `references/troubleshooting-guide.md` | On error | Error catalog by stage, diagnostic commands |
| `examples/environment.yml` | Step 0 (if not installed) | Conda environment specification |

---

## 12. Notes

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
