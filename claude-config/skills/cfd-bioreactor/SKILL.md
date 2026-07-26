---
name: cfd-bioreactor
last_updated: 2026-06-09
description: >
  Use when user needs to simulate fluid flow through bioprocess cartridges,
  bioreactor geometries, or membrane devices using FEniCSx. Covers Navier-Stokes
  flow, O2 species transport with Michaelis-Menten kinetics, membrane permeation,
  mesh generation from CAD geometry, and interactive 3D visualization.

# v2.0 orchestrator handoff metadata
handoff:
  accepts_from:
    - "*"
  provides_to:
    - cfd-mathematician
    - cfd-reviewer
    - programming-pm
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

## Role and Delegation

You orchestrate the simulation, owning session setup, state, user communication, and code
generation. Code generation stays in your loop because it is a mechanical translation from a
validated mathematical specification into FEniCSx Python using the parameterized patterns in
`references/fenicsx-patterns.md` -- the mathematician produces the specification, you fill
the template. Delegate specialist work when the subtask is substantial and independent:

- `cfd-mathematician` -- variational formulations, function space selection, stability
  conditions (inf-sup), convergence rate estimation, dimensionless analysis.
- `cfd-reviewer` -- adversarial review of plans and generated code, error diagnosis,
  physical plausibility.
- Perspective agents (FULL mode) -- 3-5 parallel domain viewpoints plus a synthesis agent.
  See Multi-Perspective Analysis below.

Handle it directly when you could finish it in a handful of tool calls, when the work is
sequential, or when you need the result in your own loop. See
`../references/delegation-and-scope.md`.

**Graceful degradation**: if `cfd-mathematician` is unavailable, use reference file defaults
for mathematical decisions. If `cfd-reviewer` is unavailable, proceed with
APPROVED_WITH_WARNINGS and log "NO ADVERSARIAL REVIEW". If fewer than 3 perspective agents
return, skip swarm synthesis and continue as LITE mode for that phase.

**Platform constraint**: agents invoked via Task inherit Read, Write, Edit, Bash, Glob,
Grep, WebSearch, and WebFetch. They cannot use the Task tool (sub-agents cannot spawn
sub-agents) or AskUserQuestion (they cannot prompt the user). They load reference files
with Read and emit handoff YAML with Write.

**Prompt scaffolds**: `references/agent-prompt-templates.md` holds the Task-tool templates
-- Section 1 (base template for mathematician and reviewer), Section 2 (inline-reference
fallback, used when an agent reports `TOOL_UNAVAILABLE: Read`), Sections 3-4 (perspective
and synthesis agents). Fill the `{...}` placeholders and pass the result to Task. Which
reference sections each agent should load is in `references/agent-loading-guide.md`.

---

## Scope

Use this skill to simulate fluid flow through bioprocess cartridges, bioreactors, or
membrane devices; generate FEniCSx Python scripts for Navier-Stokes + species transport;
model O2 consumption with Michaelis-Menten kinetics or membrane permeation with Fick's law
Robin boundary conditions; import CAD geometry (STEP/IGES) for mesh generation via gmsh;
validate results against analytical solutions (Poiseuille flow, 1D diffusion-reaction);
create interactive 3D visualizations of flow and concentration fields; and run mesh
convergence studies for publication-quality results.

**Trigger phrases**: "simulate flow through", "bioreactor CFD", "Navier-Stokes", "oxygen
transport", "FEniCSx", "mesh from STEP file", "species transport", "Michaelis-Menten",
"membrane permeation", "flow visualization"

Out of scope: compressible flow and Mach effects (Ma << 0.01 here); turbulence modeling
(laminar flow assumed, Re < 100 -- suggest OpenFOAM for turbulent flows); structural
mechanics and fluid-structure interaction (use a dedicated FEM/FSI tool); thermal effects
(isothermal throughout); HPC job submission and cloud orchestration (this skill runs
locally or hands over scripts); commercial solver workflows (COMSOL, ANSYS, OpenFOAM have
their own interfaces); transient simulations with adaptive time-stepping (steady-state
only in v2.0). For a quick order-of-magnitude estimate use the `calculator` skill instead
-- a back-of-envelope calculation takes 5 minutes, a full CFD simulation 30 minutes to
hours.

---

## Modes and Complexity Tiers

Prefix responses with the current phase so the user can see where the workflow stands, for
example `[Phase 2/5 - Flow Planning] Mathematician analysis received. Invoking reviewer.`

| Mode | Agents used | Typical time | Triggered by |
|---|---|---|---|
| DIRECT | none | ~10 min | Tier 1 Poiseuille validation, or user override |
| LITE | mathematician + reviewer | ~20-25 min | Tier 2 2D problems without hidden complexity |
| FULL | perspectives + mathematician + reviewer | ~40-60 min | Tier 3-4, user override, or complexity flags in Tier 2 |

In DIRECT mode, suppress agent-related status messages; the user sees a simple linear
workflow (Phase 0 -> mesh + flow code gen -> validation -> Phase 5). Hidden complexity that
argues for upgrading a Tier 2 problem to FULL: Pe > 100, Re > 10 (Newton continuation
needed), more than one membrane surface, or channel width less than 10x the estimated
boundary layer thickness.

| Tier | Problem | Default mode | Example |
|---|---|---|---|
| 1 | 2D Poiseuille flow; tests the pipeline with no user-specific physics | DIRECT | `examples/01-2d-channel-flow.md` |
| 2 | 2D user geometry + O2 transport with Michaelis-Menten and membrane BCs | LITE | `examples/02-2d-oxygen-transport.md` |
| 3 | 2D or 3D STEP geometry + full coupled physics, iterative solvers, Newton continuation | FULL | `examples/03-3d-cartridge-template.md` |
| 4 | Full 3D cartridge with convergence study, publication visualization, reproducibility metadata; run with MPI | FULL | -- |

Suggest starting at Tier 1 when the user has never run a FEniCSx simulation, and progressing
upward rather than jumping straight to Tier 3 or 4.

**Run purpose**: `debug` (first run, error recovery, parameter exploration) or `production`
(publication quality). Tier 1-2 are always debug. Tier 3-4 start at debug and move to
production once a debug run has completed successfully. Debug runs do not need FULL mode --
downgrade to LITE rather than spend agent time planning a mesh that gets discarded.
Production runs benefit from reviewer input -- upgrade DIRECT to LITE. Show any
auto-adjustment with its rationale; the user can override. Mesh sizing and expected runtimes
per tier and run purpose are in `references/runtime-and-mesh-sizing.md`.

---

## Phase 0 -- Session Setup

**Pre-flight environment check.** Read `references/environment-setup.md` Section 4,
generate `preflight_check.py`, and run it via Bash. Interpret the result:

| Pre-Flight Result | Action |
|---|---|
| All PASS | Proceed to workflow |
| FAIL on dolfinx | Provide micromamba/Docker install instructions from environment-setup.md. Stop. |
| FAIL on petsc4py | Cannot solve. Provide install instructions. Stop. |
| FAIL on gmsh OCC | Parametric geometry still works; warn that STEP import is disabled. |
| WARN on PyVista | Proceed; export to VTK instead of interactive plots. |
| WARN on MUMPS | Proceed; use iterative solver (GMRES+ILU) instead of direct solver. |
| WARN on RAM < 8 GB | Proceed for 2D only; warn that 3D will require coarse meshes. |

When dolfinx or petsc4py is missing, help the user set up the environment instead of
generating simulation scripts.

**Skill base path.** The directory containing this SKILL.md is the skill base path:
`SKILL_BASE = "/Users/davidangelesalbores/repos/claude/claude-config/skills/cfd-bioreactor"`.
Store it as `skill_base_path` in the session state file; orchestrator Reads and agent
invocations use absolute paths of the form `{skill_base_path}/references/{filename}`.

**Session directory.** Create `/tmp/cfd-bioreactor-session-{YYYYMMDD-HHMMSS}-{PID}/` with
subdirectories `handoffs/`, `scripts/`, `logs/`, `perspectives/`, and initialize the state
file (schema and resume behavior: `references/session-state-and-handoff.md`).

**Problem parsing and mode selection.** Extract geometry type, fluid properties, species
parameters, boundary conditions, and target tier from the request. Compute Re, Pe, and Da
where parameters allow. Select mode and run purpose, then present them with the tier and
time estimates so the user can override.

---

## Phases 1-3 -- Planning, Code Generation, Execution

Three deliverables, each with its own handoff YAML in `{session_dir}/handoffs/` (contracts
and field names: `references/orchestrator-handoff-schema.md`). In FULL mode each planning
step is preceded by a perspective swarm; in LITE mode it is mathematician then reviewer; in
DIRECT mode there is no agent involvement. Advance a phase when its conditions hold, noting
in the handoff anything that did not. Re-running phase N invalidates every later phase, so
their results have to be regenerated.

### Phase 1 -- Mesh Plan

Ask `cfd-mathematician` for element type and order, required mesh density, and boundary
layer resolution, given the geometry and target error tolerance. Ask `cfd-reviewer` to
challenge mesh quality thresholds, refinement adequacy, memory feasibility, and STEP unit
assumptions. Before generating mesh code, confirm the plan has a handoff YAML with its
required fields, a memory estimate below 70% of available RAM, element count > 0, a
specified mesh quality threshold (min scaled Jacobian > 0.1), and reviewer
`approval_status` of APPROVED or APPROVED_WITH_WARNINGS with no CRITICAL-severity
challenges.

Handoff: `{session_dir}/handoffs/phase1-mesh-plan.yaml`

### Phase 2 -- Flow

Ask `cfd-mathematician` for the weak form, verification of Taylor-Hood inf-sup stability,
expected convergence order, and a continuation strategy if Re > 10. Ask `cfd-reviewer` to
challenge BC physical consistency, solver convergence criteria, the presence of a pressure
reference point, the Newton vs. Picard choice, and memory estimates.

Generate the mesh and flow scripts per the Code Generation Protocol, then execute
(Bash for 2D, script-for-manual-run for 3D). Flow results are valid when the solver
converged (residual < tolerance), mass conservation holds
(`|integral u.n ds_inlet + integral u.n ds_outlet| < 1e-6`), no NaN or Inf appears in the
solution arrays, there are no unphysical flow patterns (reversed flow at the inlet, etc.),
and -- where a Poiseuille comparison applies -- L2 error < 1% (P2 elements give ~1e-12 for
quadratic solutions; see `references/validation-benchmarks.md` Benchmark 1).

Handoff: `{session_dir}/handoffs/phase2-flow-result.yaml`

### Phase 3 -- Transport

Ask `cfd-mathematician` for the transport variational form with SUPG, analysis of the
stabilization parameter (Pe-dependent xi formula), the convergence properties of the
regularized Michaelis-Menten term, and the expected transport convergence order. Ask
`cfd-reviewer` to challenge the SUPG parameter for overflow risk (Pe > 710 with
cosh/sinh), the adequacy of the MM regularization epsilon, the Robin BC formulation for
the membrane, and Newton convergence expectations.

Generate and execute the transport script. Transport results are valid when the Newton
solver converged, species conservation holds
(`|inlet_flux + outlet_flux + membrane_flux + total_reaction| < 1%` of total reaction), and
min(c) is not significantly negative -- negative concentrations suggest insufficient mesh
resolution near sinks.

Handoff: `{session_dir}/handoffs/phase3-transport-result.yaml`

### When the reviewer rejects a plan

The reviewer's `blocking_issues` become HARD CONSTRAINTS for one mathematician retry. If
the plan is still REJECTED, present the disagreement to the user -- each side's
recommendation and concern in a sentence -- with the options of following the
mathematician, following the reviewer, taking a conservative compromise, or supplying their
own parameters. The conservative compromise is the more cautious choice from both sides:
finer mesh, more stable solver, lower tolerance, more regularization.

---

## Multi-Perspective Analysis (FULL Mode)

Spawn perspective agents directly rather than delegating to `brainstorming-pm`: sub-agents
cannot use the Task tool, so brainstorming-pm invoked as a sub-agent would collapse to
single-context sequential operation and lose the diversity that makes swarm discussion
worthwhile. Read the relevant challenge template from
`references/swarm-framing-templates.md` (Swarm 1 for mesh, 2 for flow, 3 for transport),
fill in the problem parameters, and spawn 3-5
perspective agents as parallel Task calls using the template in
`references/agent-prompt-templates.md` Section 3. Then read the outputs from
`{session_dir}/perspectives/phase{N}-*.md` and spawn one synthesis agent (Section 4 of the
same file) to produce the swarm-synthesis handoff. The five perspectives (Numerical
Analyst, Mesh Engineer, Physical Modeler, Computational Pragmatist, Validation
Strategist), their prompts, and the phase-to-perspective mapping are in
`references/swarm-framing-templates.md` Sections 6a-6b. Keep spawn counts to what the
phase needs; see `../references/delegation-and-scope.md`.

If fewer than 3 perspectives return, skip synthesis and proceed as LITE mode for that phase.
If the synthesis agent fails, read the perspective files and summarize them yourself. Swarm
output is worth using when it contains at least two specific numerical recommendations in
`convergent_insights` and one concrete alternative in `divergent_alternatives`; otherwise log
that it fell below the threshold and proceed without it. A synthesis `confidence_score` below
3 of 5 is worth mentioning to the user as informational.

---

## Phases 4-5 -- Post-Processing and Summary

Both are orchestrator-owned. Read `references/fenicsx-patterns.md` Section 12 (PyVista)
and generate a visualization script covering velocity magnitude contours, the O2
concentration field with depletion zones highlighted, streamlines from the inlet (3D),
cross-section slices at several positions (3D), a centerline line probe of O2 vs.
position, PNG screenshots, and VTK output for ParaView. Include the headless fallback:
`if not pyvista.system_supports_plotting(): pyvista.OFF_SCREEN = True`. Execute 2D
visualizations; save 3D scripts for manual execution.

The mesh convergence study is optional -- `references/runtime-and-mesh-sizing.md` covers
when it is worth running and how to judge the result. When it runs, follow the Mesh
Convergence Study Protocol in `references/validation-benchmarks.md`: solve on 4 levels (h, h/2, h/4,
h/8; 3 levels for 3D), compute the QoIs (outlet average O2, max velocity, min O2 in the
cell region), Richardson-extrapolate the grid-independent value, and report the convergence
rate `p = log(e_{i-1}/e_i) / log(2)`, flagging non-monotone convergence.

For the Phase 5 summary, collect the handoff documents from `{session_dir}/handoffs/` and
report what was simulated (geometry, physics, parameters), the key agent decisions
(mathematician recommendations, reviewer concerns and how they were resolved), validation
results, and the files produced (scripts, XDMF checkpoints, VTK files, PNG images). Offer to
preserve or clean up the session directory, and set `status: "complete"` in the state file.

---

## Code Generation Protocol

These rules guard against code-level errors that produce silently wrong results, so they stay
inline rather than in a reference file.

### Include in every script

1. **Version assertion header** at the top:
   ```python
   import importlib.metadata as _meta
   _ver = tuple(int(x) for x in _meta.version("fenics-dolfinx").split(".")[:2])
   assert _ver >= (0, 10), f"Requires FEniCSx >= 0.10. See references/environment-setup.md."
   ```

2. **Reproducibility metadata header**: date, dolfinx version, basix version, gmsh
   version, numpy version, Python version, OS, all physical parameters.
3. **Complete, self-contained scripts** that run standalone from the command line, not
   fragments requiring manual assembly.
4. **Inline comments** explaining each step: what it does, why it is needed, what the
   expected output is.
5. **Regularized Michaelis-Menten**: `c_pos = (c + sqrt(c^2 + eps^2)) / 2` with
   `eps = 1e-10 * c_inlet`, instead of the raw `c / (Km + c)` form.
6. **SUPG stabilization** for all species transport (check Pe first, but enable by
   default), using the numerically stable form
   `xi = conditional(gt(Pe, 1.0), 1.0 - 1.0/Pe, Pe/3.0)` rather than cosh/sinh, which
   overflows for Pe > 710.
7. **Post-solve quality checks**: mass conservation, species conservation, negative
   concentration warning, NaN/Inf detection.
8. **Save results** to XDMF (checkpointing) and VTK/PVD (PyVista/ParaView):
   ```python
   with io.XDMFFile(comm, "results.xdmf", "w") as xdmf:
       xdmf.write_mesh(domain)
       xdmf.write_function(uh, 0.0)
   ```

### Avoid

- `from dolfin import *` (legacy FEniCS, not FEniCSx)
- `dolfinx.io.gmshio` (v0.10 renamed it to `dolfinx.io.gmsh`)
- Raw Michaelis-Menten `c / (Km + c)` without regularization
- `cosh(Pe)/sinh(Pe)` for the SUPG parameter (numerical overflow)
- Skipping physical group definition before meshing
- Omitting `gmsh.model.occ.synchronize()` after OCC operations
- Equal-order elements (P1/P1) for flow without stabilization
- Omitting the pressure reference point (at least one pressure BC is needed)

### Mesh script sequence

Load the reference sections listed in `references/agent-loading-guide.md` for mesh code
generation, then generate in this order.

1. Size the mesh per tier and `run_purpose` (`references/runtime-and-mesh-sizing.md`): coarse
   range for debug, fine range for production. Keep >= 3 elements across the thinnest
   geometric feature (membranes, thin channels); if the debug range cannot, raise the count
   to the minimum that can and warn the user.
2. Import STEP geometry via `gmsh.model.occ.importShapes()` (with error handling) or
   construct parametric geometry via the gmsh built-in/OCC kernel.
3. Call `gmsh.model.occ.synchronize()` after all OCC operations.
4. Define physical groups (inlet, outlet, walls, membrane, cell_region).
5. Apply mesh refinement near walls and membrane interfaces.
6. Estimate memory requirements before meshing.
7. Check mesh quality after generation.
8. Convert to DOLFINx via `dolfinx.io.gmsh.model_to_mesh()`.

### Flow script sequence

1. Define the Taylor-Hood P2/P1 function space (standard inf-sup stable pair).
2. Apply boundary conditions from the physical groups.
3. For Re < 1, assemble the Stokes variational form (linear solve). For Re >= 1, use
   Navier-Stokes with Newton continuation: solve Stokes for the initial guess, ramp Re
   through intermediate values [1, 10, 50, target] if Re > 10, using the previous
   solution as the initial guess at each step.
4. Select the solver: MUMPS (direct) below 50K DOFs, GMRES+ILU (iterative) above.
5. Include a solver progress monitor for 3D.
6. Save velocity and pressure fields to an XDMF checkpoint.

### Transport script sequence

1. Define the P1 function space for concentration, with any transport-specific
   refinement (near-membrane refinement for SUPG accuracy) sized per `run_purpose`.
2. Estimate the Peclet number: `Pe = |u| * h / (2D)`.
3. Implement SUPG stabilization and regularized Michaelis-Menten per the rules above.
4. Implement membrane permeation as a Robin BC (Fick's law).
5. Solve with Newton (the MM term makes it nonlinear), importing `NonlinearProblem` from
   `dolfinx.fem.petsc` and `NewtonSolver` from `dolfinx.nls.petsc`.
6. Include post-solve checks: negative concentration warning, species conservation.

### Pattern library and syntax check

`references/fenicsx-patterns.md` is the primary source for all generated code; its patterns
are validated against the FEniCSx v0.10 API. Bioreactor-specific constants (fluid
properties, diffusion coefficients, Michaelis-Menten kinetics by cell type) come from the
lookup tables in `references/physics-models.md`. Execution strategy by problem size, including
Bash timeouts and MPI instructions, is in `references/runtime-and-mesh-sizing.md`. After
writing a script, syntax-check it:
`python3 -c "import ast; ast.parse(open('script.py').read())"`

---

## Error Handling

`references/troubleshooting-guide.md` catalogs errors by stage: Stage 0 environment
(ImportError, ModuleNotFoundError -- run the pre-flight check), Stage 1 geometry and mesh
(OCC exceptions, 0 volumes after STEP import, mesh quality, meshing OOM), Stage 2-3
physics setup (TypeError in BCs, shape mismatch, missing pressure reference), Stage 4
solver execution (Newton divergence, MUMPS pivot errors, PETSc OOM, negative
concentrations), Stage 5 validation (conservation violated), Stage 5-6 visualization
(empty plot, PyVista crash).

For Newton divergence, the recovery order is: Stokes solution as initial guess (most
common fix), ramp Re through [1, 10, 50, target], reduce the relaxation factor to 0.5 then
0.3 then 0.1, switch to Picard iteration, refine the mesh near boundaries and
high-gradient regions, and if all of that fails suggest simplifying the geometry or
reducing Re. Give actionable diagnostic messages rather than bare error codes.

### Self-correction loop

When execution or validation fails, collect the error output (traceback, log messages,
solver convergence history) and classify it: `import_error`, `mesh_error`,
`solver_divergence`, `assertion_failure`, `numerical_instability`, `wall_clock_timeout`.

- Simple errors (import path, syntax, assertion): fix directly and retry, up to 2 direct
  retries.
- Complex errors (solver divergence, numerical instability): invoke `cfd-reviewer` in
  error diagnosis mode with the error output, the `error_history` from previous attempts,
  and the mathematical specification. The reviewer classifies the root cause and
  recommends a fix; `error_history` keeps it from re-recommending fixes that already
  failed.
- Wall-clock timeout past the kill threshold in `references/runtime-and-mesh-sizing.md`: on a
  debug run the mesh is likely too fine, so reduce the element count by 50% and retry without
  involving the reviewer. On a production run, invoke `cfd-reviewer` -- common causes are
  insufficient Newton continuation steps or solver parameters needing tuning.

Regenerate the code with the fix applied and append to `error_history` (error_type,
error_message, fix_attempted, fix_outcome). After 3 total retries, stop and report the error
type, each attempt and its outcome, and a recommended next step.

On the orchestration side, a malformed or unparseable handoff warrants one retry of that
agent with explicit format feedback, then proceeding without its input. A missing session
directory can be recreated from the state YAML; tell the user if it cannot be.

---

## Context Management

Load reference file sections per phase rather than all at once. Agent discussion detail stays
in the handoff files on disk; carry forward only key decisions, approval status, and
validation metrics (roughly 200 tokens per completed phase), so each phase works from the
previous phase's handoff plus its own.

---

## Integration with Other Skills

| Need | Skill | How to Use |
|---|---|---|
| Quick feasibility estimate | `calculator` | Before full CFD: estimate Re, Pe, Da, O2 depletion depth. A 5-minute calculation can determine if CFD is even needed. |
| Cross-check CFD results | `calculator` | After simulation: verify results against analytical estimates. Disagreement by > 10x warrants investigation. |
| Literature parameter values | `bioinformatician` or `researcher` | Find Vmax, Km for specific cell types; diffusion coefficients in specific media. |
| Debug generated notebook | `notebook-debugger` | If the user is running in Jupyter and hits FEniCSx import/kernel issues. |
| Transition to software project | `programming-pm` | When simulation scripts evolve into a maintained software project (library packaging, tests/CI, multi-file project, parameter sweep framework). Criteria: `references/session-state-and-handoff.md`. |

---

## Reference Files

Read these as needed rather than all at once; section-level loading maps per agent are in
`references/agent-loading-guide.md`.

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
| `references/swarm-framing-templates.md` | FULL mode swarms | Pre-written challenge templates for perspective agents |
| `references/agent-prompt-templates.md` | Agent invocations | Task-tool prompt scaffolds (base, fallback, perspective, synthesis) |
| `references/runtime-and-mesh-sizing.md` | Phases 0-4 | Mesh sizing, wall-clock limits, execution strategy, convergence diagnostics |
| `references/session-state-and-handoff.md` | Phase 0, Phase 5 | State file schema, resume protocol, programming-pm handoff criteria |
| `examples/environment.yml` | Phase 0 (if not installed) | Conda environment specification |

---

## Notes

- **Version pinning**: all code patterns target FEniCSx v0.10; the version assertion catches
  a mismatch. Update the patterns before using them with a new API version.
- **SI units everywhere**: parameters are in SI (kg, m, s, mol, Pa) unless the user says
  otherwise. If a STEP file appears to be in mm (bounding box > 1 m), warn the user about
  the potential unit mismatch.
- **Physical group convention**: inlet=1, outlet=2, walls=3, membrane=4, cell_region=5,
  fluid_volume=10. Consistent across all examples and generated scripts.
- **gmsh synchronize**: after any `gmsh.model.occ.*` operation (addBox, addCylinder,
  importShapes, etc.), call `gmsh.model.occ.synchronize()` before accessing entities,
  assigning physical groups, or meshing. Forgetting this is the single most common gmsh
  error.
- **Serial execution by default**: generated scripts run on a single core but are MPI-aware
  (they use `MPI.COMM_WORLD` throughout) and work with `mpirun -np N`. Add MPI instructions
  for 3D production runs.
- **Reproducibility**: every generated script carries a metadata header with version numbers,
  parameters, and settings, so re-running it in the same environment reproduces the results.
- **Generates code, does not modify the system**: this skill writes Python scripts to the
  user's working directory. It never installs packages, modifies micromamba environments, or
  changes system configuration.
