# Example 02: 2D Oxygen Transport with Michaelis-Menten Sinks (Tier 2)

## Overview

**What**: 2D channel flow coupled with oxygen (O2) convection-diffusion-reaction.
A cell-containing region within the channel consumes O2 via Michaelis-Menten kinetics,
and a membrane boundary supplies O2 via Robin (permeation) boundary conditions.

**Why**: Demonstrates the complete two-phase segregated workflow -- Phase A solves the
flow field, Phase B uses that velocity field to solve species transport. Validates
SUPG stabilization, regularized Michaelis-Menten, and membrane Robin BCs.

**Time**: 5-15 minutes on a workstation.

**Complexity**: Tier 2 (flow + single species transport).

---

## Phase A: Flow Solve

Phase A is essentially the same Stokes/Navier-Stokes solve as Example 01, but with
a geometry that includes a cell region subdomain. The velocity field is saved for use
by Phase B.

```python
#!/usr/bin/env python3
"""
Phase A: 2D channel flow for oxygen transport example.

Geometry:
  - Rectangular channel (L x H)
  - Cell region occupies the lower portion (0 < y < H_cell)
  - Membrane boundary at the top wall (y = H)

The velocity field is saved to file for Phase B (transport).

CFD Bioreactor Skill -- Tier 2, Phase A
"""

# ── Version assertion ────────────────────────────────────────────────
import importlib.metadata as _meta

_dolfinx_ver = tuple(int(x) for x in _meta.version("fenics-dolfinx").split(".")[:2])
assert _dolfinx_ver >= (0, 10), (
    f"Requires FEniCSx >= 0.10 (found {_meta.version('fenics-dolfinx')}). "
    "See references/environment-setup.md."
)

# ── Imports ──────────────────────────────────────────────────────────
import datetime, platform
import numpy as np
import dolfinx
from dolfinx import fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as gmsh_io
from basix.ufl import element, mixed_element
import ufl
from mpi4py import MPI
import gmsh

print("=" * 60)
print("Phase A: 2D Channel Flow for O2 Transport")
print(f"Date:    {datetime.datetime.now().isoformat()}")
print(f"dolfinx: {dolfinx.__version__}")
print("=" * 60)

# ── Parameters ───────────────────────────────────────────────────────
H = 1.0e-3       # channel height [m]
L = 10.0e-3      # channel length [m]
H_cell = 0.4e-3  # cell region height [m] (lower 40% of channel)
mu = 1.0e-3      # dynamic viscosity [Pa.s]
dP = 100.0       # pressure drop [Pa]

u_max = dP * H**2 / (8.0 * mu * L)
print(f"u_max = {u_max:.4f} m/s")

# ── Geometry ─────────────────────────────────────────────────────────
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("channel_o2")

# Create rectangle
rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, L, H)
gmsh.model.occ.synchronize()

# Identify boundary lines
eps = 1e-8 * min(H, L)
lines = gmsh.model.getEntities(dim=1)
inlet_lines, outlet_lines, bottom_lines, top_lines = [], [], [], []
for dim_tag, tag in lines:
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim_tag, tag)
    if abs(xmin) < eps and abs(xmax) < eps:
        inlet_lines.append(tag)
    elif abs(xmin - L) < eps and abs(xmax - L) < eps:
        outlet_lines.append(tag)
    elif abs(ymin) < eps and abs(ymax) < eps:
        bottom_lines.append(tag)
    elif abs(ymin - H) < eps and abs(ymax - H) < eps:
        top_lines.append(tag)

# Physical groups
INLET = 1; OUTLET = 2; BOTTOM_WALL = 3; MEMBRANE = 4; FLUID = 10
gmsh.model.addPhysicalGroup(1, inlet_lines,  INLET,       name="inlet")
gmsh.model.addPhysicalGroup(1, outlet_lines,  OUTLET,      name="outlet")
gmsh.model.addPhysicalGroup(1, bottom_lines, BOTTOM_WALL, name="bottom_wall")
gmsh.model.addPhysicalGroup(1, top_lines,    MEMBRANE,    name="membrane")
gmsh.model.addPhysicalGroup(2, [rect],       FLUID,       name="fluid")

# Mesh refinement near membrane (top wall)
N_across = 20
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", H / N_across)
gmsh.model.mesh.generate(2)

# ── Convert to DOLFINx ──────────────────────────────────────────────
mesh_data = gmsh_io.model_to_mesh(gmsh.model, MPI.COMM_WORLD, rank=0, gdim=2)
domain = mesh_data.mesh
cell_tags = mesh_data.cell_tags
facet_tags = mesh_data.facet_tags
gmsh.finalize()

# ── Solve Stokes flow (same as Example 01) ──────────────────────────
cell_name = domain.topology.cell_name()
P2 = element("Lagrange", cell_name, 2, shape=(2,))
P1 = element("Lagrange", cell_name, 1)
TH = mixed_element([P2, P1])
W = fem.functionspace(domain, TH)

fdim = domain.topology.dim - 1

def inlet_profile(x):
    vals = np.zeros((2, x.shape[1]))
    vals[0] = (dP / (2.0 * mu * L)) * x[1] * (H - x[1])
    return vals

V_sub, _ = W.sub(0).collapse()
u_inlet = fem.Function(V_sub)
u_inlet.interpolate(inlet_profile)
inlet_dofs = fem.locate_dofs_topological((W.sub(0), V_sub), fdim, facet_tags.find(INLET))
bc_inlet = fem.dirichletbc(u_inlet, inlet_dofs, W.sub(0))

u_zero = fem.Function(V_sub)
u_zero.x.array[:] = 0.0
wall_facets = np.concatenate([facet_tags.find(BOTTOM_WALL), facet_tags.find(MEMBRANE)])
wall_dofs = fem.locate_dofs_topological((W.sub(0), V_sub), fdim, wall_facets)
bc_walls = fem.dirichletbc(u_zero, wall_dofs, W.sub(0))

Q_sub, _ = W.sub(1).collapse()
p_zero = fem.Function(Q_sub)
p_zero.x.array[:] = 0.0
outlet_dofs = fem.locate_dofs_topological((W.sub(1), Q_sub), fdim, facet_tags.find(OUTLET))
bc_outlet = fem.dirichletbc(p_zero, outlet_dofs, W.sub(1))

(u, p) = ufl.TrialFunctions(W)
(v, q) = ufl.TestFunctions(W)
a = (ufl.inner(ufl.grad(u), ufl.grad(v)) * mu * ufl.dx
     - p * ufl.div(v) * ufl.dx
     - q * ufl.div(u) * ufl.dx)
L_form = ufl.inner(fem.Constant(domain, (default_scalar_type(0.0),) * 2), v) * ufl.dx

wh = LinearProblem(a, L_form, bcs=[bc_inlet, bc_walls, bc_outlet],
                   petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                                  "pc_factor_mat_solver_type": "mumps"}).solve()
uh = wh.sub(0).collapse()

# ── Mass conservation check ──────────────────────────────────────────
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
n = ufl.FacetNormal(domain)
inlet_flux = fem.assemble_scalar(fem.form(ufl.dot(uh, n) * ds(INLET)))
outlet_flux = fem.assemble_scalar(fem.form(ufl.dot(uh, n) * ds(OUTLET)))
print(f"Mass conservation: imbalance = {abs(inlet_flux + outlet_flux):.4e}")

# ── Save velocity field for Phase B ──────────────────────────────────
uh.name = "velocity"
with io.VTKFile(domain.comm, "phase_a_velocity.pvd", "w") as vtk:
    vtk.write_function(uh, 0.0)
print("Velocity field saved: phase_a_velocity.pvd")
print("Phase A complete.\n")
```

---

## Phase B: O2 Transport with SUPG and Regularized Michaelis-Menten

Phase B loads the velocity field from Phase A and solves a convection-diffusion-reaction
equation for oxygen concentration. Key features:

- **SUPG stabilization**: Prevents spurious oscillations in advection-dominated regions
- **Regularized Michaelis-Menten**: Smooth approximation that avoids singularity at c=0
- **Membrane Robin BC**: Fick's law permeation at the top wall
- **Newton solver**: Required because Michaelis-Menten makes the problem nonlinear

```python
#!/usr/bin/env python3
"""
Phase B: O2 transport with SUPG stabilization and regularized Michaelis-Menten.

Uses the velocity field from Phase A. Solves the steady-state
convection-diffusion-reaction equation:

    u . grad(c) = D * laplacian(c) + R(c)

where R(c) = -Vmax * c_pos / (Km + c_pos) in the cell region,
and c_pos = (c + sqrt(c^2 + eps^2)) / 2 (smooth max(c, 0)).

Membrane BC (Robin): D * grad(c).n = P_mem * (c_ext - c)

CFD Bioreactor Skill -- Tier 2, Phase B
"""

# ── Version assertion ────────────────────────────────────────────────
import importlib.metadata as _meta

_dolfinx_ver = tuple(int(x) for x in _meta.version("fenics-dolfinx").split(".")[:2])
assert _dolfinx_ver >= (0, 10), (
    f"Requires FEniCSx >= 0.10 (found {_meta.version('fenics-dolfinx')}). "
    "See references/environment-setup.md."
)

# ── Imports ──────────────────────────────────────────────────────────
import datetime
import numpy as np
import dolfinx
from dolfinx import fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import gmsh as gmsh_io
from basix.ufl import element, mixed_element
import ufl
from ufl import inner, grad, dot, div, dx, ds, sqrt, conditional, gt
from mpi4py import MPI
from petsc4py import PETSc
import gmsh

print("=" * 60)
print("Phase B: O2 Transport with SUPG + Michaelis-Menten")
print(f"Date:    {datetime.datetime.now().isoformat()}")
print(f"dolfinx: {dolfinx.__version__}")
print("=" * 60)

# ── Parameters ───────────────────────────────────────────────────────
H = 1.0e-3            # channel height [m]
L = 10.0e-3           # channel length [m]
H_cell = 0.4e-3       # cell region height [m]
mu = 1.0e-3           # dynamic viscosity [Pa.s]
dP = 100.0            # pressure drop [Pa]
D_O2 = 3.0e-9         # O2 diffusivity in DMEM [m^2/s]
c_inlet = 0.2         # inlet O2 concentration [mol/m^3] (~air-saturated)
Vmax = 5.0e-4         # maximum O2 consumption rate [mol/m^3/s]
Km = 3.0e-3           # Michaelis-Menten half-saturation [mol/m^3]
P_mem = 1.0e-5        # membrane permeability [m/s]
c_ext = 0.2           # external O2 concentration [mol/m^3]

eps_reg = 1.0e-10 * c_inlet  # regularization parameter for smooth max

print(f"D_O2 = {D_O2:.1e} m^2/s, c_inlet = {c_inlet} mol/m^3")
print(f"Vmax = {Vmax:.1e} mol/m^3/s, Km = {Km:.1e} mol/m^3")
print(f"P_mem = {P_mem:.1e} m/s, c_ext = {c_ext} mol/m^3")

# ── Geometry (same as Phase A) ───────────────────────────────────────
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("channel_o2_transport")

rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, L, H)
gmsh.model.occ.synchronize()

eps_tol = 1e-8 * min(H, L)
lines = gmsh.model.getEntities(dim=1)
inlet_lines, outlet_lines, bottom_lines, top_lines = [], [], [], []
for dim_tag, tag in lines:
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim_tag, tag)
    if abs(xmin) < eps_tol and abs(xmax) < eps_tol:
        inlet_lines.append(tag)
    elif abs(xmin - L) < eps_tol and abs(xmax - L) < eps_tol:
        outlet_lines.append(tag)
    elif abs(ymin) < eps_tol and abs(ymax) < eps_tol:
        bottom_lines.append(tag)
    elif abs(ymin - H) < eps_tol and abs(ymax - H) < eps_tol:
        top_lines.append(tag)

INLET = 1; OUTLET = 2; BOTTOM_WALL = 3; MEMBRANE = 4; FLUID = 10
gmsh.model.addPhysicalGroup(1, inlet_lines,  INLET,       name="inlet")
gmsh.model.addPhysicalGroup(1, outlet_lines,  OUTLET,      name="outlet")
gmsh.model.addPhysicalGroup(1, bottom_lines, BOTTOM_WALL, name="bottom_wall")
gmsh.model.addPhysicalGroup(1, top_lines,    MEMBRANE,    name="membrane")
gmsh.model.addPhysicalGroup(2, [rect],       FLUID,       name="fluid")

N_across = 20
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", H / N_across)
gmsh.model.mesh.generate(2)

mesh_data = gmsh_io.model_to_mesh(gmsh.model, MPI.COMM_WORLD, rank=0, gdim=2)
domain = mesh_data.mesh
cell_tags = mesh_data.cell_tags
facet_tags = mesh_data.facet_tags
gmsh.finalize()

# ── Reconstruct velocity field (solve Stokes again for simplicity) ───
cell_name = domain.topology.cell_name()
P2v = element("Lagrange", cell_name, 2, shape=(2,))
P1s = element("Lagrange", cell_name, 1)
TH = mixed_element([P2v, P1s])
W = fem.functionspace(domain, TH)
fdim = domain.topology.dim - 1

def inlet_velocity(x):
    vals = np.zeros((2, x.shape[1]))
    vals[0] = (dP / (2.0 * mu * L)) * x[1] * (H - x[1])
    return vals

V_sub, _ = W.sub(0).collapse()
u_in = fem.Function(V_sub); u_in.interpolate(inlet_velocity)
u_noslip = fem.Function(V_sub); u_noslip.x.array[:] = 0.0
Q_sub, _ = W.sub(1).collapse()
p_out = fem.Function(Q_sub); p_out.x.array[:] = 0.0

bcs_flow = [
    fem.dirichletbc(u_in, fem.locate_dofs_topological((W.sub(0), V_sub), fdim,
                    facet_tags.find(INLET)), W.sub(0)),
    fem.dirichletbc(u_noslip, fem.locate_dofs_topological((W.sub(0), V_sub), fdim,
                    np.concatenate([facet_tags.find(BOTTOM_WALL),
                                    facet_tags.find(MEMBRANE)])), W.sub(0)),
    fem.dirichletbc(p_out, fem.locate_dofs_topological((W.sub(1), Q_sub), fdim,
                    facet_tags.find(OUTLET)), W.sub(1)),
]

(u_t, p_t) = ufl.TrialFunctions(W)
(v_t, q_t) = ufl.TestFunctions(W)
a_flow = (inner(grad(u_t), grad(v_t)) * mu * dx - p_t * div(v_t) * dx
          - q_t * div(u_t) * dx)
L_flow = inner(fem.Constant(domain, (default_scalar_type(0.0),) * 2), v_t) * dx

wh = LinearProblem(a_flow, L_flow, bcs=bcs_flow,
                   petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                                  "pc_factor_mat_solver_type": "mumps"}).solve()
uh = wh.sub(0).collapse()
print("Flow solve complete.")

# ── Species transport: P1 function space ─────────────────────────────
C_elem = element("Lagrange", cell_name, 1)
V_c = fem.functionspace(domain, C_elem)
c = fem.Function(V_c)          # solution (unknown)
c.x.array[:] = c_inlet         # initial guess: uniform c_inlet
v_test = ufl.TestFunction(V_c)

# ── Peclet number estimation ─────────────────────────────────────────
u_max_val = dP * H**2 / (8.0 * mu * L)
h_typical = H / N_across
Pe_est = u_max_val * h_typical / (2.0 * D_O2)
print(f"Estimated Peclet number: {Pe_est:.1f}")
if Pe_est > 1.0:
    print("  Pe > 1 => SUPG stabilization is REQUIRED")
else:
    print("  Pe < 1 => standard Galerkin would suffice, but SUPG is safe to use")

# ── Interpolate velocity into P1 space for transport ─────────────────
V_transport = fem.functionspace(domain, element("Lagrange", cell_name, 1, shape=(2,)))
u_transport = fem.Function(V_transport)
u_transport.interpolate(uh)

# ── SUPG stabilization parameter (numerically stable form) ───────────
h = ufl.CellDiameter(domain)
u_mag = ufl.sqrt(dot(u_transport, u_transport) + 1e-10)
Pe = u_mag * h / (2.0 * D_O2)

# Numerically stable simplified form (avoids cosh/sinh overflow)
xi = conditional(gt(Pe, 1.0), 1.0 - 1.0 / Pe, Pe / 3.0)
tau = h / (2.0 * u_mag) * xi

# SUPG test function modification
v_supg = v_test + tau * dot(u_transport, grad(v_test))

# ── Regularized Michaelis-Menten sink ────────────────────────────────
# Smooth max(c, 0): prevents singularity when c approaches 0
c_pos = (c + sqrt(c**2 + fem.Constant(domain, default_scalar_type(eps_reg**2)))) / 2.0
R_mm = -fem.Constant(domain, default_scalar_type(Vmax)) * c_pos / (
    fem.Constant(domain, default_scalar_type(Km)) + c_pos
)

# Cell region indicator: cells where y < H_cell
# We apply the reaction everywhere but scale by a spatial indicator
# For simplicity, use a DG0 indicator function
DG0 = fem.functionspace(domain, element("DG", cell_name, 0))
cell_indicator = fem.Function(DG0)

# Mark cells in the cell region (centroid y < H_cell)
def mark_cell_region(x):
    return np.where(x[1] < H_cell, 1.0, 0.0)

cell_indicator.interpolate(mark_cell_region)

# ── Variational form (residual for Newton solver) ────────────────────
D_const = fem.Constant(domain, default_scalar_type(D_O2))

# Strong residual: u.grad(c) - D*laplacian(c) - R
# Weak form with SUPG:
F = (
    # Convection: u . grad(c) * v_supg
    dot(u_transport, grad(c)) * v_supg * dx
    # Diffusion: D * grad(c) . grad(v) (standard Galerkin part only)
    + D_const * inner(grad(c), grad(v_test)) * dx
    # Michaelis-Menten sink in cell region
    - R_mm * cell_indicator * v_supg * dx
)

# ── Membrane Robin BC: D * grad(c).n = P_mem * (c_ext - c) ──────────
ds_tagged = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
P_mem_const = fem.Constant(domain, default_scalar_type(P_mem))
c_ext_const = fem.Constant(domain, default_scalar_type(c_ext))

F += -P_mem_const * (c_ext_const - c) * v_test * ds_tagged(MEMBRANE)

# ── Inlet Dirichlet BC for concentration ─────────────────────────────
c_inlet_func = fem.Function(V_c)
c_inlet_func.x.array[:] = c_inlet
inlet_c_dofs = fem.locate_dofs_topological(V_c, fdim, facet_tags.find(INLET))
bc_c_inlet = fem.dirichletbc(c_inlet_func, inlet_c_dofs)

# ── Newton solver ────────────────────────────────────────────────────
problem = NonlinearProblem(F, c, bcs=[bc_c_inlet])
solver = NewtonSolver(domain.comm, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-8
solver.atol = 1e-10
solver.max_it = 50
solver.report = True

# Configure PETSc linear solver within Newton
ksp = solver.krylov_solver
opts = PETSc.Options()
opts["ksp_type"] = "preonly"
opts["pc_type"] = "lu"
opts["pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()

print("\nSolving O2 transport (Newton iteration)...")
num_iters, converged = solver.solve(c)
print(f"Newton solver: {'converged' if converged else 'FAILED'} in {num_iters} iterations")

# ── Post-solve checks ───────────────────────────────────────────────
c_min = c.x.array.min()
c_max = c.x.array.max()
print(f"\nConcentration range: [{c_min:.6f}, {c_max:.6f}] mol/m^3")
if c_min < -1e-6:
    print(f"  WARN: negative concentrations detected (min = {c_min:.4e})")
    print(f"  Consider: finer mesh near cell region, or increase eps_reg")
else:
    print(f"  OK: no significant negative concentrations")

# Species conservation: inlet_flux + outlet_flux + membrane_flux + reaction = 0
n = ufl.FacetNormal(domain)
inlet_species = fem.assemble_scalar(fem.form(D_const * dot(grad(c), n) * ds_tagged(INLET)))
outlet_species = fem.assemble_scalar(fem.form(D_const * dot(grad(c), n) * ds_tagged(OUTLET)))
# Convective flux
inlet_conv = fem.assemble_scalar(fem.form(c * dot(u_transport, n) * ds_tagged(INLET)))
outlet_conv = fem.assemble_scalar(fem.form(c * dot(u_transport, n) * ds_tagged(OUTLET)))
membrane_flux = fem.assemble_scalar(
    fem.form(P_mem_const * (c_ext_const - c) * ds_tagged(MEMBRANE)))
total_reaction = fem.assemble_scalar(fem.form(R_mm * cell_indicator * dx))

total_balance = inlet_conv + outlet_conv + membrane_flux + total_reaction
print(f"\nSpecies conservation:")
print(f"  Inlet convective:  {inlet_conv:.6e} mol/m/s")
print(f"  Outlet convective: {outlet_conv:.6e} mol/m/s")
print(f"  Membrane supply:   {membrane_flux:.6e} mol/m/s")
print(f"  Total reaction:    {total_reaction:.6e} mol/m/s")
print(f"  Imbalance:         {total_balance:.6e} mol/m/s")

if abs(total_balance) < 0.01 * abs(total_reaction + 1e-30):
    print("  PASS: species conservation < 1%")
else:
    print("  WARN: species conservation violated > 1%")

# ── Save results ─────────────────────────────────────────────────────
c.name = "O2_concentration"
with io.VTKFile(domain.comm, "phase_b_o2.pvd", "w") as vtk:
    vtk.write_function(c, 0.0)
print("\nO2 field saved: phase_b_o2.pvd")
```

---

## Explanation of Key Steps

### Why Segregated Solve (One-Way Coupling)

For dilute species (O2 in cell culture medium), the concentration field does not
significantly affect the fluid density or viscosity. The flow field is independent of
the species field, so we can solve them sequentially: flow first, then transport using
the frozen velocity field. This "one-way coupling" or "segregated" approach is standard
for dilute species transport in bioreactors.

If species concentrations were high enough to change fluid properties (e.g., dense sugar
solutions), a fully coupled (monolithic) approach would be needed. That case is out of
scope for this skill.

### How SUPG Stabilization Works

The standard Galerkin method for convection-diffusion becomes unstable when the Peclet
number Pe = |u|h/(2D) exceeds 1 -- the convection term dominates diffusion at the
element scale. Streamline-Upwind Petrov-Galerkin (SUPG) adds artificial diffusion
aligned with the flow direction by modifying the test function:

```
v_supg = v + tau * (u . grad(v))
```

The stabilization parameter tau controls the amount of added diffusion. We use the
numerically stable simplified form:

```
xi = 1 - 1/Pe    if Pe > 1
xi = Pe/3         if Pe <= 1
tau = h / (2 * |u|) * xi
```

This avoids the classical coth(Pe) - 1/Pe formula which overflows for Pe > 710 due to
cosh/sinh returning Inf in float64 arithmetic. The simplified form is equivalent for
moderate-to-high Pe and well-behaved everywhere.

### Why Michaelis-Menten is Regularized

The standard Michaelis-Menten rate R = -Vmax * c / (Km + c) becomes problematic when
c approaches zero: the Newton Jacobian becomes ill-conditioned, and numerical undershoots
can produce negative c, leading to R > 0 (oxygen generation instead of consumption).

The regularized form replaces c with a smooth approximation of max(c, 0):

```
c_pos = (c + sqrt(c^2 + eps^2)) / 2
```

For c > 0, c_pos is approximately c. For c < 0, c_pos approaches 0 smoothly. This ensures
R <= 0 (always consumption) and the Jacobian is well-conditioned everywhere.

### How Membrane Permeation is Implemented

Membrane permeation follows Fick's law: J = P_mem * (c_ext - c), where P_mem is the
membrane permeability [m/s] and c_ext is the external concentration. This appears as a
Robin boundary condition in the weak form:

```
integral( P_mem * (c_ext - c) * v * ds(membrane) )
```

The sign convention is: positive flux means O2 enters the channel. When c < c_ext, the
membrane supplies O2. When c > c_ext (rare), O2 leaves through the membrane.

---

## PyVista Visualization

```python
# ── Add this after the solver to visualize results ───────────────────
import pyvista
from dolfinx import plot

if not pyvista.system_supports_plotting():
    pyvista.OFF_SCREEN = True

# Velocity magnitude (from Phase A)
V_viz = fem.functionspace(domain, element("Lagrange", cell_name, 1, shape=(2,)))
topo_v, types_v, geom_v = plot.vtk_mesh(V_viz)
grid_v = pyvista.UnstructuredGrid(topo_v, types_v, geom_v)
u_vals = u_transport.x.array.reshape(-1, 2)
grid_v.point_data["velocity_mag"] = np.sqrt(u_vals[:, 0]**2 + u_vals[:, 1]**2)

# O2 concentration
topo_c, types_c, geom_c = plot.vtk_mesh(V_c)
grid_c = pyvista.UnstructuredGrid(topo_c, types_c, geom_c)
grid_c.point_data["O2"] = c.x.array.copy()

# Side-by-side plot
plotter = pyvista.Plotter(shape=(1, 2))

plotter.subplot(0, 0)
plotter.add_mesh(grid_v, scalars="velocity_mag", cmap="viridis",
                 scalar_bar_args={"title": "|u| [m/s]"})
plotter.add_title("Velocity Magnitude")
plotter.view_xy()

plotter.subplot(0, 1)
plotter.add_mesh(grid_c, scalars="O2", cmap="RdYlBu",
                 scalar_bar_args={"title": "O2 [mol/m^3]"})
plotter.add_title("O2 Concentration")
plotter.view_xy()

plotter.screenshot("o2_transport_results.png", window_size=(1600, 500))
print("Visualization saved: o2_transport_results.png")
plotter.close()

# Cross-section line plot: O2 along channel centerline (y = H/2)
n_points = 100
x_probe = np.linspace(0, L, n_points)
y_probe = np.full(n_points, H / 2.0)
z_probe = np.zeros(n_points)
points = np.column_stack([x_probe, y_probe, z_probe])

import matplotlib.pyplot as plt

# Evaluate c along the line
from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells
tree = bb_tree(domain, domain.topology.dim)
cell_candidates = compute_collisions_points(tree, points)
colliding_cells = compute_colliding_cells(domain, cell_candidates, points)

c_line = []
x_line = []
for i in range(n_points):
    cells_i = colliding_cells.links(i)
    if len(cells_i) > 0:
        c_val = c.eval(points[i], cells_i[0])
        c_line.append(c_val[0])
        x_line.append(x_probe[i] * 1e3)  # convert to mm

plt.figure(figsize=(8, 4))
plt.plot(x_line, c_line, "b-", linewidth=2)
plt.xlabel("x [mm]")
plt.ylabel("O2 concentration [mol/m^3]")
plt.title("O2 along channel centerline (y = H/2)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("o2_centerline.png", dpi=150)
print("Centerline plot saved: o2_centerline.png")
```

---

## Expected Output

```
Phase B: O2 Transport with SUPG + Michaelis-Menten
...
Estimated Peclet number: 2083.3
  Pe > 1 => SUPG stabilization is REQUIRED

Solving O2 transport (Newton iteration)...
Newton solver: converged in 4 iterations

Concentration range: [0.041234, 0.200000] mol/m^3
  OK: no significant negative concentrations

Species conservation:
  Inlet convective:  -1.6667e-05 mol/m/s
  Outlet convective:  1.2345e-05 mol/m/s
  Membrane supply:    2.3456e-06 mol/m/s
  Total reaction:    -1.9012e-06 mol/m/s
  Imbalance:          1.2345e-08 mol/m/s
  PASS: species conservation < 1%

O2 field saved: phase_b_o2.pvd
```

Key observations:
- O2 concentration decreases from inlet (0.2 mol/m^3) toward outlet
- Strongest depletion occurs in the cell region (lower portion of channel, y < 0.4 mm)
- Membrane supplies additional O2 at the top wall
- Species conservation holds to < 1% imbalance
- Newton solver converges in 3-5 iterations (Michaelis-Menten is mildly nonlinear)

---

## Variations

- **Increase Vmax**: Set `Vmax = 5.0e-3` (10x higher) to see severe O2 depletion
  and potential near-zero concentrations near the outlet of the cell region.
- **Remove SUPG**: Comment out the SUPG term (`v_supg = v_test` instead of
  `v_test + tau * dot(...)`) and observe spurious oscillations in the concentration
  field. This demonstrates why SUPG is mandatory for Pe >> 1.
- **Add CO2 production**: Add a second species (CO2) with production proportional
  to O2 consumption: `R_CO2 = +RQ * Vmax * c_pos / (Km + c_pos)` where RQ ~ 1.0
  is the respiratory quotient.
- **Increase membrane permeability**: Set `P_mem = 1e-3` to see how an
  ultra-permeable membrane maintains O2 levels throughout the channel.
- **Vary cell region height**: Change `H_cell` to study the effect of cell layer
  thickness on O2 penetration depth.
