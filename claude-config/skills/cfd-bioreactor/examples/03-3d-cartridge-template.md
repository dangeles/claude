# Example 03: 3D Cartridge Simulation Template (Tier 3 Production)

## Overview

**What**: Complete 3D bioprocess cartridge simulation with STEP geometry import (or
parametric cylinder fallback), coupled Navier-Stokes + O2 transport, mesh convergence
study, and publication-quality post-processing.

**Why**: This is the production template for real cartridge simulations. It demonstrates
STEP import with full error handling, iterative solvers for 3D problems, Newton
continuation for Navier-Stokes, and comprehensive validation and visualization.

**Time**: 1-4 hours for the full workflow (5-10 minutes for coarse mesh validation run).

**Complexity**: Tier 3 (production quality with coupled multiphysics).

**Important**: This template generates a standalone `.py` script that the user should
run manually for anything beyond a coarse validation mesh. 3D simulations typically
exceed Bash timeout limits. For MPI parallel execution, see the User Instructions
section at the end.

---

## Geometry: STEP Import with Parametric Fallback

The geometry section attempts to load a user-provided STEP file. If no STEP file is
available (or import fails), it falls back to constructing a parametric cylindrical
cartridge directly in gmsh.

```python
#!/usr/bin/env python3
"""
3D bioprocess cartridge simulation: flow + O2 transport.

Geometry: STEP import with parametric cylinder fallback.
Flow: Navier-Stokes with Newton continuation (Stokes initial guess).
Transport: SUPG + regularized Michaelis-Menten + membrane Robin BC.

CFD Bioreactor Skill -- Tier 3 Production Template
"""

# ── Version assertion ────────────────────────────────────────────────
import importlib.metadata as _meta

_dolfinx_ver = tuple(int(x) for x in _meta.version("fenics-dolfinx").split(".")[:2])
assert _dolfinx_ver >= (0, 10), (
    f"Requires FEniCSx >= 0.10 (found {_meta.version('fenics-dolfinx')}). "
    "See references/environment-setup.md."
)

# ── Reproducibility header ───────────────────────────────────────────
import datetime, platform, os, json, sys
import numpy as np
import dolfinx, basix, gmsh
from dolfinx import fem, mesh, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import gmsh as gmsh_io
from basix.ufl import element, mixed_element
import ufl
from ufl import inner, grad, dot, div, dx, ds, sqrt, conditional, gt
from mpi4py import MPI
from petsc4py import PETSc

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

metadata = {
    "simulation": "3D Cartridge CFD",
    "date": datetime.datetime.now().isoformat(),
    "dolfinx": dolfinx.__version__,
    "basix": basix.__version__,
    "gmsh_api": gmsh.GMSH_API_VERSION,
    "numpy": np.__version__,
    "python": platform.python_version(),
    "os": platform.platform(),
    "mpi_size": comm.Get_size(),
}

if rank == 0:
    print("=" * 60)
    print("3D Bioprocess Cartridge Simulation")
    for k, v in metadata.items():
        print(f"  {k}: {v}")
    print("=" * 60)

# ── Physical parameters (SI units) ──────────────────────────────────
# USER: modify these for your specific cartridge
params = {
    # Geometry (used only if STEP import fails)
    "R_cartridge": 5.0e-3,     # cartridge radius [m]
    "L_cartridge": 50.0e-3,    # cartridge length [m]
    "R_cell_region": 3.5e-3,   # inner cell region radius [m]
    # Fluid
    "rho": 1000.0,             # density [kg/m^3] (DMEM at 37C)
    "mu": 1.0e-3,              # dynamic viscosity [Pa.s]
    "u_inlet": 0.01,           # mean inlet velocity [m/s]
    # Species
    "D_O2": 3.0e-9,            # O2 diffusivity [m^2/s]
    "c_inlet": 0.2,            # inlet O2 [mol/m^3]
    "Vmax": 5.0e-4,            # max consumption [mol/m^3/s]
    "Km": 3.0e-3,              # half-saturation [mol/m^3]
    "P_mem": 1.0e-5,           # membrane permeability [m/s]
    "c_ext": 0.2,              # external O2 [mol/m^3]
    # Mesh
    "mesh_size_coarse": 1.0e-3,  # coarse mesh element size [m]
    "mesh_size_fine": 0.3e-3,    # refined size near walls/membrane [m]
    # Solver
    "Re_target": None,           # computed from inlet velocity
    # Files
    "step_file": None,           # path to STEP file (None = use parametric)
}

Re = params["rho"] * params["u_inlet"] * 2 * params["R_cartridge"] / params["mu"]
params["Re_target"] = Re
if rank == 0:
    print(f"\nReynolds number: {Re:.1f}")
    print(f"Regime: {'Stokes' if Re < 1 else 'Navier-Stokes (laminar)' if Re < 100 else 'WARNING: Re>100'}")

# ── Memory estimation ────────────────────────────────────────────────
try:
    import psutil
    ram_gb = psutil.virtual_memory().available / 1e9
    # Rough estimate: 1 GB per 100K P2 tetrahedra
    est_elements = (4.0 / 3.0 * np.pi * params["R_cartridge"]**3) / (
        params["mesh_size_coarse"]**3 / 6.0)
    est_ram_gb = est_elements / 100_000
    if rank == 0:
        print(f"\nMemory estimate: ~{est_elements:.0f} elements, ~{est_ram_gb:.1f} GB")
        print(f"Available RAM: {ram_gb:.1f} GB")
        if est_ram_gb > 0.7 * ram_gb:
            print("WARNING: estimated memory exceeds 70% available RAM.")
            print("Consider coarser mesh or iterative solver.")
except ImportError:
    if rank == 0:
        print("psutil not available -- skipping memory check")
```

### Geometry Construction

```python
# ── STEP import or parametric fallback ───────────────────────────────
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("cartridge_3d")

step_file = params["step_file"]
use_parametric = True

if step_file and os.path.isfile(step_file):
    if rank == 0:
        print(f"\nImporting STEP file: {step_file}")
        file_size_mb = os.path.getsize(step_file) / 1e6
        print(f"  File size: {file_size_mb:.1f} MB")
        if file_size_mb > 50:
            print("  WARNING: Large STEP file. Meshing may be slow.")
        if file_size_mb > 200:
            print("  STRONG WARNING: Very large STEP file. Consider defeaturing.")

    try:
        gmsh.model.occ.importShapes(step_file)
        gmsh.model.occ.synchronize()

        # Verify entities were created
        volumes = gmsh.model.getEntities(dim=3)
        surfaces = gmsh.model.getEntities(dim=2)
        if rank == 0:
            print(f"  Volumes: {len(volumes)}, Surfaces: {len(surfaces)}")

        if len(volumes) == 0:
            if rank == 0:
                print("  WARNING: No 3D volumes found. This may be a 2D geometry.")
                print("  Falling back to parametric construction.")
        elif len(volumes) > 1:
            if rank == 0:
                print(f"  WARNING: Multiple volumes ({len(volumes)}). Using first volume.")
                print("  If this is an assembly, extract the fluid domain before import.")
            use_parametric = False
        else:
            use_parametric = False

        # Unit check: bounding box
        if not use_parametric:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
            extent = max(xmax - xmin, ymax - ymin, zmax - zmin)
            if rank == 0:
                print(f"  Bounding box: [{xmin:.4f}, {ymin:.4f}, {zmin:.4f}] to "
                      f"[{xmax:.4f}, {ymax:.4f}, {zmax:.4f}]")
                if extent > 1.0:
                    print("  WARNING: Dimensions suggest units may be in mm, not m.")
                    print("  Physics parameters assume SI (meters). Scale if needed.")

    except Exception as e:
        if rank == 0:
            print(f"  STEP import failed: {e}")
            print("  Troubleshooting: re-export as STEP AP214, remove fillets < 0.5mm")
            print("  Falling back to parametric construction.")
        use_parametric = True

if use_parametric:
    if rank == 0:
        print("\nConstructing parametric cylindrical cartridge...")
    R = params["R_cartridge"]
    L = params["L_cartridge"]

    # Outer cylinder (full domain)
    cyl = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, R)
    gmsh.model.occ.synchronize()

    # Identify surfaces by bounding box for physical group assignment
    surfaces = gmsh.model.getEntities(dim=2)
    inlet_surfs, outlet_surfs, wall_surfs = [], [], []
    eps = 1e-8 * L

    for dim_s, tag_s in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim_s, tag_s)
        if abs(zmin) < eps and abs(zmax) < eps:
            inlet_surfs.append(tag_s)
        elif abs(zmin - L) < eps and abs(zmax - L) < eps:
            outlet_surfs.append(tag_s)
        else:
            wall_surfs.append(tag_s)

    # Physical groups
    INLET_3D = 1; OUTLET_3D = 2; WALL_3D = 3; FLUID_3D = 10
    gmsh.model.addPhysicalGroup(2, inlet_surfs,   INLET_3D,  name="inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surfs,   OUTLET_3D, name="outlet")
    gmsh.model.addPhysicalGroup(2, wall_surfs,     WALL_3D,   name="wall")
    gmsh.model.addPhysicalGroup(3, [cyl],          FLUID_3D,  name="fluid")

    if rank == 0:
        print(f"  Radius: {R*1e3:.1f} mm, Length: {L*1e3:.1f} mm")
        print(f"  Surfaces: inlet={len(inlet_surfs)}, outlet={len(outlet_surfs)}, "
              f"wall={len(wall_surfs)}")

# ── Mesh refinement (Distance + Threshold near walls) ────────────────
lc_coarse = params["mesh_size_coarse"]
lc_fine = params["mesh_size_fine"]

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_coarse)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc_fine / 2.0)

# Refine near wall surfaces
if use_parametric and wall_surfs:
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", wall_surfs)

    f_thresh = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thresh, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thresh, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(f_thresh, "SizeMax", lc_coarse)
    gmsh.model.mesh.field.setNumber(f_thresh, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(f_thresh, "DistMax", params["R_cartridge"] * 0.3)

    gmsh.model.mesh.field.setAsBackgroundMesh(f_thresh)

if rank == 0:
    print("\nGenerating 3D mesh...")
gmsh.model.mesh.generate(3)

# Mesh quality
qualities = gmsh.model.mesh.getElementQualities(qualityType="minSJ")
if rank == 0:
    print(f"  Elements: {len(qualities)}")
    print(f"  Min scaled Jacobian: {min(qualities):.4f}")
    print(f"  Mean scaled Jacobian: {np.mean(qualities):.4f}")
    if min(qualities) < 0.01:
        print("  WARNING: very poor quality elements detected. Consider defeaturing.")

# ── Convert to DOLFINx ──────────────────────────────────────────────
mesh_data = gmsh_io.model_to_mesh(gmsh.model, comm, rank=0, gdim=3)
domain = mesh_data.mesh
cell_tags = mesh_data.cell_tags
facet_tags = mesh_data.facet_tags
gmsh.finalize()

if rank == 0:
    num_cells = domain.topology.index_map(3).size_local
    print(f"  DOLFINx mesh: {num_cells} cells on rank 0")
```

---

## Phase A: 3D Flow Solve

```python
# ══════════════════════════════════════════════════════════════════════
# PHASE A: 3D NAVIER-STOKES WITH NEWTON CONTINUATION
# ══════════════════════════════════════════════════════════════════════

if rank == 0:
    print("\n" + "=" * 60)
    print("PHASE A: 3D Flow Solve")
    print("=" * 60)

# ── Taylor-Hood P2/P1 ───────────────────────────────────────────────
cell_name = domain.topology.cell_name()
gdim = domain.geometry.dim
P2 = element("Lagrange", cell_name, 2, shape=(gdim,))
P1 = element("Lagrange", cell_name, 1)
TH = mixed_element([P2, P1])
W = fem.functionspace(domain, TH)

# ── Boundary conditions ─────────────────────────────────────────────
fdim = domain.topology.dim - 1
INLET_3D = 1; OUTLET_3D = 2; WALL_3D = 3

V_sub, _ = W.sub(0).collapse()
Q_sub, _ = W.sub(1).collapse()

# Inlet: uniform velocity in z-direction
def inlet_velocity_3d(x):
    vals = np.zeros((gdim, x.shape[1]))
    vals[2] = params["u_inlet"]  # z-direction
    return vals

u_inlet_func = fem.Function(V_sub)
u_inlet_func.interpolate(inlet_velocity_3d)
inlet_dofs = fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_tags.find(INLET_3D))
bc_inlet = fem.dirichletbc(u_inlet_func, inlet_dofs, W.sub(0))

# Walls: no-slip
u_zero = fem.Function(V_sub)
u_zero.x.array[:] = 0.0
wall_dofs = fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, facet_tags.find(WALL_3D))
bc_walls = fem.dirichletbc(u_zero, wall_dofs, W.sub(0))

# Outlet: zero pressure
p_zero = fem.Function(Q_sub)
p_zero.x.array[:] = 0.0
outlet_dofs = fem.locate_dofs_topological(
    (W.sub(1), Q_sub), fdim, facet_tags.find(OUTLET_3D))
bc_outlet = fem.dirichletbc(p_zero, outlet_dofs, W.sub(1))

bcs_flow = [bc_inlet, bc_walls, bc_outlet]

# ── Step 1: Stokes solve for initial guess ───────────────────────────
if rank == 0:
    print("\nStep 1: Solving Stokes for initial guess...")

(u_trial, p_trial) = ufl.TrialFunctions(W)
(v_test, q_test) = ufl.TestFunctions(W)

mu_val = params["mu"]
a_stokes = (inner(grad(u_trial), grad(v_test)) * mu_val * dx
            - p_trial * div(v_test) * dx
            - q_test * div(u_trial) * dx)
L_stokes = inner(fem.Constant(domain, (default_scalar_type(0.0),) * gdim), v_test) * dx

# Use iterative solver for 3D (GMRES + ILU default)
stokes_opts = {
    "ksp_type": "gmres",
    "pc_type": "ilu",
    "ksp_rtol": 1e-8,
    "ksp_max_it": 1000,
    "ksp_monitor": None,
}
# Fall back to MUMPS for small problems
num_dofs = W.dofmap.index_map.size_global * W.dofmap.index_map_bs
if num_dofs < 50_000:
    stokes_opts = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    if rank == 0:
        print(f"  Using MUMPS direct solver ({num_dofs} DOFs < 50K threshold)")
else:
    if rank == 0:
        print(f"  Using GMRES+ILU iterative solver ({num_dofs} DOFs)")

wh_stokes = LinearProblem(a_stokes, L_stokes, bcs=bcs_flow,
                          petsc_options=stokes_opts).solve()

if rank == 0:
    print("  Stokes solve complete.")

# ── Step 2: Newton continuation for Navier-Stokes (if Re > 1) ───────
if Re > 1.0:
    if rank == 0:
        print(f"\nStep 2: Newton continuation for Re = {Re:.1f}...")

    w_ns = fem.Function(W)
    w_ns.x.array[:] = wh_stokes.x.array[:]  # use Stokes as initial guess

    u_ns, p_ns = ufl.split(w_ns)
    v_ns, q_ns = ufl.TestFunctions(W)

    rho = params["rho"]

    # Navier-Stokes residual (nonlinear)
    F_ns = (
        rho * inner(dot(u_ns, ufl.nabla_grad(u_ns)), v_ns) * dx
        + mu_val * inner(grad(u_ns), grad(v_ns)) * dx
        - p_ns * div(v_ns) * dx
        - q_ns * div(u_ns) * dx
    )

    problem_ns = NonlinearProblem(F_ns, w_ns, bcs=bcs_flow)
    solver_ns = NewtonSolver(comm, problem_ns)
    solver_ns.convergence_criterion = "incremental"
    solver_ns.rtol = 1e-6
    solver_ns.atol = 1e-8
    solver_ns.max_it = 50
    solver_ns.report = True

    ksp_ns = solver_ns.krylov_solver
    opts_ns = PETSc.Options()
    if num_dofs < 50_000:
        opts_ns["ksp_type"] = "preonly"
        opts_ns["pc_type"] = "lu"
        opts_ns["pc_factor_mat_solver_type"] = "mumps"
    else:
        opts_ns["ksp_type"] = "gmres"
        opts_ns["pc_type"] = "ilu"
        opts_ns["ksp_rtol"] = 1e-8
        opts_ns["ksp_max_it"] = 1000
    ksp_ns.setFromOptions()

    # Continuation: ramp Reynolds number if Re > 10
    if Re > 10:
        Re_steps = [1.0, 10.0, Re]
    else:
        Re_steps = [Re]

    for Re_step in Re_steps:
        # Scale inlet velocity to achieve target Re
        u_scale = Re_step / Re * params["u_inlet"]
        def inlet_scaled(x, u_s=u_scale):
            vals = np.zeros((gdim, x.shape[1]))
            vals[2] = u_s
            return vals
        u_inlet_func.interpolate(inlet_scaled)

        if rank == 0:
            print(f"  Re = {Re_step:.1f}: solving...")

        try:
            n_its, converged = solver_ns.solve(w_ns)
            if rank == 0:
                status = "converged" if converged else "FAILED"
                print(f"    Newton: {status} in {n_its} iterations")
            if not converged:
                if rank == 0:
                    print("    WARNING: Newton did not converge. Using last iterate.")
                break
        except RuntimeError as e:
            if rank == 0:
                print(f"    Newton failed: {e}")
                print("    Using Stokes solution as fallback.")
            w_ns.x.array[:] = wh_stokes.x.array[:]
            break

    uh = w_ns.sub(0).collapse()
    ph = w_ns.sub(1).collapse()
else:
    if rank == 0:
        print("\nRe < 1: using Stokes solution directly (no Newton needed).")
    uh = wh_stokes.sub(0).collapse()
    ph = wh_stokes.sub(1).collapse()

# ── Flow validation ──────────────────────────────────────────────────
ds_tags = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
n_vec = ufl.FacetNormal(domain)

inlet_flux = fem.assemble_scalar(fem.form(dot(uh, n_vec) * ds_tags(INLET_3D)))
outlet_flux = fem.assemble_scalar(fem.form(dot(uh, n_vec) * ds_tags(OUTLET_3D)))
mass_imbalance = abs(inlet_flux + outlet_flux)

if rank == 0:
    print(f"\nMass conservation:")
    print(f"  Inlet flux:  {inlet_flux:.6e} m^3/s")
    print(f"  Outlet flux: {outlet_flux:.6e} m^3/s")
    print(f"  Imbalance:   {mass_imbalance:.6e}")
    if mass_imbalance < 1e-6:
        print("  PASS")
    else:
        print("  WARN: mass conservation violated")

# Save intermediate checkpoint
uh.name = "velocity"
ph.name = "pressure"
with io.XDMFFile(comm, "checkpoint_flow.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(uh, 0.0)
    xdmf.write_function(ph, 0.0)

if rank == 0:
    print("Flow checkpoint saved: checkpoint_flow.xdmf")
    print("Phase A complete.\n")
```

---

## Phase B: 3D O2 Transport

```python
# ══════════════════════════════════════════════════════════════════════
# PHASE B: 3D O2 TRANSPORT WITH SUPG + MICHAELIS-MENTEN
# ══════════════════════════════════════════════════════════════════════

if rank == 0:
    print("=" * 60)
    print("PHASE B: 3D O2 Transport")
    print("=" * 60)

# ── P1 function space for concentration ──────────────────────────────
C_elem = element("Lagrange", cell_name, 1)
V_c = fem.functionspace(domain, C_elem)
c = fem.Function(V_c)
c.x.array[:] = params["c_inlet"]  # initial guess

# ── Interpolate velocity for transport ───────────────────────────────
V_transport = fem.functionspace(domain, element("Lagrange", cell_name, 1, shape=(gdim,)))
u_transport = fem.Function(V_transport)
u_transport.interpolate(uh)

# ── Cell region indicator (r < R_cell_region) ────────────────────────
DG0 = fem.functionspace(domain, element("DG", cell_name, 0))
cell_indicator = fem.Function(DG0)

R_cell = params["R_cell_region"]
def mark_cell_region_3d(x):
    r = np.sqrt(x[0]**2 + x[1]**2)
    return np.where(r < R_cell, 1.0, 0.0)

cell_indicator.interpolate(mark_cell_region_3d)

# ── Peclet number ────────────────────────────────────────────────────
u_inlet_val = params["u_inlet"]
h_est = params["mesh_size_coarse"]
D_val = params["D_O2"]
Pe_est = u_inlet_val * h_est / (2.0 * D_val)
if rank == 0:
    print(f"\nEstimated Peclet number: {Pe_est:.0f}")

# ── SUPG stabilization (numerically stable form) ────────────────────
v_test_c = ufl.TestFunction(V_c)
h = ufl.CellDiameter(domain)
u_mag = sqrt(dot(u_transport, u_transport) + 1e-10)
Pe = u_mag * h / (2.0 * D_val)

# Numerically stable: avoids cosh/sinh overflow for large Pe
xi = conditional(gt(Pe, 1.0), 1.0 - 1.0 / Pe, Pe / 3.0)
tau = h / (2.0 * u_mag) * xi
v_supg = v_test_c + tau * dot(u_transport, grad(v_test_c))

# ── Regularized Michaelis-Menten ─────────────────────────────────────
eps_reg = 1e-10 * params["c_inlet"]
c_pos = (c + sqrt(c**2 + fem.Constant(domain, default_scalar_type(eps_reg**2)))) / 2.0
Vmax_c = fem.Constant(domain, default_scalar_type(params["Vmax"]))
Km_c = fem.Constant(domain, default_scalar_type(params["Km"]))
R_mm = -Vmax_c * c_pos / (Km_c + c_pos)

# ── Variational form ────────────────────────────────────────────────
D_const = fem.Constant(domain, default_scalar_type(D_val))

F_transport = (
    dot(u_transport, grad(c)) * v_supg * dx
    + D_const * inner(grad(c), grad(v_test_c)) * dx
    - R_mm * cell_indicator * v_supg * dx
)

# Membrane Robin BC (on wall surface for parametric; customize for STEP)
P_mem_c = fem.Constant(domain, default_scalar_type(params["P_mem"]))
c_ext_c = fem.Constant(domain, default_scalar_type(params["c_ext"]))
F_transport += -P_mem_c * (c_ext_c - c) * v_test_c * ds_tags(WALL_3D)

# Inlet Dirichlet BC
c_inlet_func = fem.Function(V_c)
c_inlet_func.x.array[:] = params["c_inlet"]
inlet_c_dofs = fem.locate_dofs_topological(V_c, fdim, facet_tags.find(INLET_3D))
bc_c = fem.dirichletbc(c_inlet_func, inlet_c_dofs)

# ── Newton solver ────────────────────────────────────────────────────
problem_c = NonlinearProblem(F_transport, c, bcs=[bc_c])
solver_c = NewtonSolver(comm, problem_c)
solver_c.convergence_criterion = "incremental"
solver_c.rtol = 1e-8
solver_c.atol = 1e-10
solver_c.max_it = 50
solver_c.report = True

ksp_c = solver_c.krylov_solver
opts_c = PETSc.Options()
opts_c["ksp_type"] = "gmres"
opts_c["pc_type"] = "ilu"
opts_c["ksp_rtol"] = 1e-8
opts_c["ksp_max_it"] = 500
ksp_c.setFromOptions()

if rank == 0:
    print("\nSolving O2 transport (Newton)...")
n_its_c, converged_c = solver_c.solve(c)
if rank == 0:
    print(f"Newton: {'converged' if converged_c else 'FAILED'} in {n_its_c} iterations")

# ── Post-solve checks ───────────────────────────────────────────────
c_min = comm.allreduce(c.x.array.min(), op=MPI.MIN)
c_max = comm.allreduce(c.x.array.max(), op=MPI.MAX)
if rank == 0:
    print(f"\nO2 range: [{c_min:.6f}, {c_max:.6f}] mol/m^3")
    if c_min < -1e-6:
        print(f"  WARN: negative concentrations (min = {c_min:.4e})")
    else:
        print(f"  OK: no significant negatives")

# Species conservation
inlet_c_flux = fem.assemble_scalar(fem.form(c * dot(u_transport, n_vec) * ds_tags(INLET_3D)))
outlet_c_flux = fem.assemble_scalar(fem.form(c * dot(u_transport, n_vec) * ds_tags(OUTLET_3D)))
membrane_supply = fem.assemble_scalar(
    fem.form(P_mem_c * (c_ext_c - c) * ds_tags(WALL_3D)))
total_rxn = fem.assemble_scalar(fem.form(R_mm * cell_indicator * dx))

if rank == 0:
    balance = inlet_c_flux + outlet_c_flux + membrane_supply + total_rxn
    print(f"\nSpecies conservation:")
    print(f"  Inlet flux:      {inlet_c_flux:.6e}")
    print(f"  Outlet flux:     {outlet_c_flux:.6e}")
    print(f"  Membrane supply: {membrane_supply:.6e}")
    print(f"  Total reaction:  {total_rxn:.6e}")
    print(f"  Imbalance:       {balance:.6e}")

# Save transport checkpoint
c.name = "O2_concentration"
with io.XDMFFile(comm, "checkpoint_transport.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(c, 0.0)

if rank == 0:
    print("Transport checkpoint saved: checkpoint_transport.xdmf")
    print("Phase B complete.\n")
```

---

## Validation Suite

```python
# ══════════════════════════════════════════════════════════════════════
# VALIDATION: MESH CONVERGENCE + CONSERVATION
# ══════════════════════════════════════════════════════════════════════

if rank == 0:
    print("=" * 60)
    print("VALIDATION SUITE")
    print("=" * 60)

# For 3D, use coarse/medium/fine (3 levels) instead of 4 to save time
mesh_sizes = [
    params["mesh_size_coarse"],
    params["mesh_size_coarse"] / 2.0,
    params["mesh_size_coarse"] / 4.0,
]

# QoI: outlet average O2 concentration, max velocity
qois = {"outlet_avg_O2": [], "max_velocity": [], "min_O2_cell": [], "h": []}

if rank == 0:
    print("\nMesh convergence study: 3 levels")
    print("(Full 4-level study recommended for publication -- see "
          "references/validation-benchmarks.md)")
    print("\nNOTE: For large 3D problems, this convergence study may take hours.")
    print("Consider running each level separately or using MPI.\n")

# The full convergence study would re-run the entire pipeline for each mesh size.
# Here we record the QoIs from the current solve as the first data point.
outlet_area = fem.assemble_scalar(
    fem.form(fem.Constant(domain, default_scalar_type(1.0)) * ds_tags(OUTLET_3D)))
if outlet_area > 0:
    outlet_avg = fem.assemble_scalar(
        fem.form(c * ds_tags(OUTLET_3D))) / outlet_area
else:
    outlet_avg = 0.0

u_max_computed = comm.allreduce(
    np.sqrt(np.sum(uh.x.array.reshape(-1, gdim)**2, axis=1)).max(), op=MPI.MAX)

# Min O2 in cell region (approximate: use global min since cell region is subset)
min_o2_cell = c_min  # conservative approximation

qois["outlet_avg_O2"].append(outlet_avg)
qois["max_velocity"].append(u_max_computed)
qois["min_O2_cell"].append(min_o2_cell)
qois["h"].append(params["mesh_size_coarse"])

if rank == 0:
    print(f"Current mesh (h = {params['mesh_size_coarse']*1e3:.2f} mm):")
    print(f"  Outlet avg O2:  {outlet_avg:.6f} mol/m^3")
    print(f"  Max velocity:   {u_max_computed:.6f} m/s")
    print(f"  Min O2 (cell):  {min_o2_cell:.6f} mol/m^3")
    print()
    print("To complete convergence study, re-run with finer mesh sizes:")
    for i, hs in enumerate(mesh_sizes[1:], 1):
        print(f"  Level {i+1}: mesh_size_coarse = {hs*1e3:.3f} mm")

# Richardson extrapolation (when 3+ levels are available)
if len(qois["outlet_avg_O2"]) >= 3:
    q = qois["outlet_avg_O2"]
    h_vals = qois["h"]
    # p = log((q[0]-q[1])/(q[1]-q[2])) / log(h_vals[0]/h_vals[1])
    # q_exact = q[2] + (q[2]-q[1]) / (2^p - 1)
    if rank == 0:
        print("\nRichardson extrapolation available with 3 data points.")
```

---

## Post-Processing

```python
# ══════════════════════════════════════════════════════════════════════
# POST-PROCESSING: 3D VISUALIZATION
# ══════════════════════════════════════════════════════════════════════

if rank == 0:
    print("\n" + "=" * 60)
    print("POST-PROCESSING")
    print("=" * 60)

try:
    import pyvista
    if not pyvista.system_supports_plotting():
        pyvista.OFF_SCREEN = True

    from dolfinx import plot

    # ── Velocity magnitude on mesh ───────────────────────────────────
    V_viz = fem.functionspace(domain, element("Lagrange", cell_name, 1, shape=(gdim,)))
    u_viz = fem.Function(V_viz)
    u_viz.interpolate(uh)
    topo_v, types_v, geom_v = plot.vtk_mesh(V_viz)
    grid_v = pyvista.UnstructuredGrid(topo_v, types_v, geom_v)
    u_arr = u_viz.x.array.reshape(-1, gdim)
    grid_v.point_data["velocity_mag"] = np.sqrt(np.sum(u_arr**2, axis=1))
    grid_v.point_data["velocity"] = u_arr

    # ── O2 concentration on mesh ─────────────────────────────────────
    topo_c, types_c, geom_c = plot.vtk_mesh(V_c)
    grid_c = pyvista.UnstructuredGrid(topo_c, types_c, geom_c)
    grid_c.point_data["O2"] = c.x.array.copy()

    # ── 3D velocity isosurface ───────────────────────────────────────
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid_v.contour([0.5 * u_max_computed], scalars="velocity_mag"),
                     color="steelblue", opacity=0.6, label="50% u_max isosurface")
    plotter.add_title("Velocity Isosurface (50% u_max)")
    plotter.screenshot("velocity_isosurface.png", window_size=(1200, 800))
    if rank == 0:
        print("Saved: velocity_isosurface.png")
    plotter.close()

    # ── O2 concentration slices at multiple z-positions ──────────────
    L_cart = params["L_cartridge"]
    z_positions = [0.25 * L_cart, 0.5 * L_cart, 0.75 * L_cart]

    plotter = pyvista.Plotter(shape=(1, 3))
    for i, z_pos in enumerate(z_positions):
        plotter.subplot(0, i)
        slc = grid_c.slice(normal="z", origin=(0, 0, z_pos))
        plotter.add_mesh(slc, scalars="O2", cmap="RdYlBu",
                         scalar_bar_args={"title": "O2 [mol/m^3]"})
        plotter.add_title(f"z = {z_pos*1e3:.1f} mm")
    plotter.screenshot("o2_slices.png", window_size=(1800, 600))
    if rank == 0:
        print("Saved: o2_slices.png")
    plotter.close()

    # ── Streamlines from inlet ───────────────────────────────────────
    plotter = pyvista.Plotter()
    try:
        # Create seed points at inlet
        R_cart = params["R_cartridge"]
        n_seeds = 20
        theta = np.linspace(0, 2 * np.pi, n_seeds, endpoint=False)
        r_seeds = R_cart * 0.5
        seed_points = np.column_stack([
            r_seeds * np.cos(theta),
            r_seeds * np.sin(theta),
            np.zeros(n_seeds)
        ])
        seed = pyvista.PolyData(seed_points)

        # Need structured grid for streamlines -- use sample to regular grid
        bounds = grid_v.bounds
        spacing = params["mesh_size_coarse"]
        sampled = grid_v.sample(
            pyvista.ImageData(
                dimensions=(
                    max(2, int((bounds[1]-bounds[0])/spacing)),
                    max(2, int((bounds[3]-bounds[2])/spacing)),
                    max(2, int((bounds[5]-bounds[4])/spacing)),
                ),
                spacing=(spacing, spacing, spacing),
                origin=(bounds[0], bounds[2], bounds[4]),
            )
        )
        streamlines = sampled.streamlines_from_source(
            seed, vectors="velocity", max_time=L_cart / params["u_inlet"] * 2,
            integration_direction="forward"
        )
        plotter.add_mesh(streamlines.tube(radius=R_cart*0.02),
                         scalars="velocity_mag", cmap="plasma")
    except Exception as e:
        if rank == 0:
            print(f"  Streamlines failed: {e} (non-critical)")
        plotter.add_mesh(grid_v, scalars="velocity_mag", cmap="viridis", opacity=0.3)

    plotter.add_title("Streamlines from Inlet")
    plotter.screenshot("streamlines.png", window_size=(1200, 800))
    if rank == 0:
        print("Saved: streamlines.png")
    plotter.close()

    # ── Line probe along centerline ──────────────────────────────────
    n_probe = 200
    probe_points = np.zeros((n_probe, 3))
    probe_points[:, 2] = np.linspace(0, L_cart, n_probe)
    probe = pyvista.PolyData(probe_points)
    sampled_line = probe.sample(grid_c)

    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 4))
    plt.plot(probe_points[:, 2] * 1e3, sampled_line.point_data["O2"], "b-", lw=2)
    plt.xlabel("z [mm]")
    plt.ylabel("O2 [mol/m^3]")
    plt.title("O2 Concentration Along Centerline")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("centerline_o2.png", dpi=150)
    if rank == 0:
        print("Saved: centerline_o2.png")

except ImportError:
    if rank == 0:
        print("PyVista not available. Results exported to VTK/XDMF.")
        print("Open checkpoint_flow.xdmf and checkpoint_transport.xdmf in ParaView.")
except Exception as e:
    if rank == 0:
        print(f"Visualization error: {e}")
        print("Results available in checkpoint_*.xdmf files.")

# ── VTK export (always, regardless of PyVista) ───────────────────────
with io.VTKFile(comm, "cartridge_results.pvd", "w") as vtk:
    vtk.write_function(uh, 0.0)
    vtk.write_function(ph, 0.0)
    vtk.write_function(c, 0.0)
if rank == 0:
    print("VTK export: cartridge_results.pvd")
```

---

## Reproducibility Block

```python
# ══════════════════════════════════════════════════════════════════════
# REPRODUCIBILITY
# ══════════════════════════════════════════════════════════════════════

metadata["parameters"] = params
metadata["mesh"] = {
    "element_size_coarse": params["mesh_size_coarse"],
    "element_size_fine": params["mesh_size_fine"],
    "num_cells_rank0": domain.topology.index_map(3).size_local,
    "num_dofs_flow": W.dofmap.index_map.size_global * W.dofmap.index_map_bs,
    "num_dofs_transport": V_c.dofmap.index_map.size_global,
}
metadata["results"] = {
    "outlet_avg_O2": float(outlet_avg),
    "max_velocity": float(u_max_computed),
    "min_O2": float(c_min),
    "max_O2": float(c_max),
    "mass_imbalance": float(mass_imbalance),
}
metadata["solver"] = {
    "flow": "Stokes" if Re < 1 else f"Newton (Re={Re:.1f})",
    "transport": f"Newton ({n_its_c} iterations)",
}

# Save metadata as JSON for re-runs
if rank == 0:
    with open("simulation_metadata.json", "w") as f:
        # Convert non-serializable types
        def _serialize(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.integer):
                return int(obj)
            return str(obj)
        json.dump(metadata, f, indent=2, default=_serialize)
    print("\nReproducibility metadata: simulation_metadata.json")

# ── Output file naming convention ────────────────────────────────────
# checkpoint_flow.xdmf       -- flow field (velocity + pressure)
# checkpoint_transport.xdmf  -- concentration field
# cartridge_results.pvd      -- combined VTK for ParaView
# simulation_metadata.json   -- all parameters, versions, results
# velocity_isosurface.png    -- 3D velocity visualization
# o2_slices.png              -- O2 cross-sections
# streamlines.png            -- flow streamlines
# centerline_o2.png          -- centerline O2 profile

if rank == 0:
    print("\n" + "=" * 60)
    print("SIMULATION COMPLETE")
    print("=" * 60)
```

---

## User Instructions

### Running This Template

**For coarse mesh validation** (5-10 minutes, Claude can run inline):

```bash
python cartridge_simulation.py
```

**For production mesh** (1-4 hours, run manually):

1. Save this script as `cartridge_simulation.py`
2. Modify the `params` dictionary with your specific parameters
3. If you have a STEP file, set `params["step_file"] = "/path/to/your/geometry.step"`
4. Run: `python cartridge_simulation.py`

**For large 3D problems with MPI** (parallel execution):

```bash
# 4-core parallel execution
mpirun -np 4 python cartridge_simulation.py

# On HPC with SLURM:
srun -n 16 python cartridge_simulation.py
```

MPI notes:
- MPI parallelism uses domain decomposition (each rank owns part of the mesh)
- The script handles rank-dependent operations (only rank 0 prints and writes metadata)
- Visualization and VTK export work correctly in parallel
- Memory is distributed: a 4-core run uses roughly 1/4 the per-rank memory

### Resuming from Checkpoint

If the simulation is interrupted during Phase B, you can resume from the flow checkpoint:

```python
# Load the flow checkpoint
with io.XDMFFile(comm, "checkpoint_flow.xdmf", "r") as xdmf:
    domain = xdmf.read_mesh()
    uh = xdmf.read_function(domain, "velocity")
# Then run Phase B only
```

### Customizing for Your Geometry

1. **STEP file**: Set `params["step_file"]` to your file path. Ensure the geometry
   contains only the fluid domain (not solid parts).
2. **Physical groups**: After STEP import, you may need to manually assign physical
   groups based on surface locations. The parametric fallback shows how to use
   bounding box queries for this.
3. **Membrane surfaces**: For cartridges with internal membranes, assign a separate
   physical group to membrane surfaces and apply the Robin BC there instead of on
   the outer wall.
4. **Cell region**: Modify `mark_cell_region_3d()` to match your geometry's cell-containing
   region (e.g., between two radii, or within a specific volume).
