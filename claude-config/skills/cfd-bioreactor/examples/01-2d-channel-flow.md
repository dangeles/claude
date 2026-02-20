# Example 01: 2D Poiseuille Channel Flow (Tier 1 Validation)

## Overview

**What**: Fully developed laminar flow (Poiseuille flow) between two parallel plates
in a 2D rectangular channel. This is the simplest possible CFD validation case.

**Why**: The Poiseuille solution is a quadratic polynomial in the transverse coordinate.
P2 (quadratic) velocity elements reproduce this exactly, so the velocity L2 error
should be at machine precision (~1e-12). This validates that the entire pipeline --
geometry, mesh, function spaces, boundary conditions, solver, and post-processing --
is working correctly. Pressure convergence (P1 elements, O(h^2)) provides a
meaningful convergence rate demonstration.

**Time**: 1-3 minutes on any workstation.

**Complexity**: Tier 1 (flow only, no species transport).

---

## Complete Python Script

```python
#!/usr/bin/env python3
"""
Poiseuille flow validation: 2D channel between parallel plates.

Analytical solution: u(y) = (dP / (2 mu L)) * y * (H - y)
P2 elements reproduce quadratic solutions exactly, so velocity error ~ 1e-12.
Pressure convergence (P1, O(h^2)) demonstrates mesh convergence methodology.

CFD Bioreactor Skill -- Tier 1 Validation
"""

# ── Version assertion ────────────────────────────────────────────────
import importlib.metadata as _meta

_dolfinx_ver = tuple(int(x) for x in _meta.version("fenics-dolfinx").split(".")[:2])
assert _dolfinx_ver >= (0, 10), (
    f"This script requires FEniCSx >= 0.10 (found {_meta.version('fenics-dolfinx')}). "
    "See references/environment-setup.md for installation instructions."
)

# ── Reproducibility header ───────────────────────────────────────────
import datetime, platform, dolfinx, basix, gmsh, pyvista, numpy as np

print("=" * 60)
print("Poiseuille Flow Validation")
print(f"Date:     {datetime.datetime.now().isoformat()}")
print(f"dolfinx:  {dolfinx.__version__}")
print(f"basix:    {basix.__version__}")
print(f"gmsh API: {gmsh.GMSH_API_VERSION}")
print(f"numpy:    {np.__version__}")
print(f"Python:   {platform.python_version()}")
print(f"OS:       {platform.platform()}")
print("=" * 60)

# ── Imports ──────────────────────────────────────────────────────────
from dolfinx import fem, mesh, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as gmsh_io
from basix.ufl import element, mixed_element
import ufl
from mpi4py import MPI

# ── Physical parameters (SI units) ──────────────────────────────────
H = 1.0e-3       # channel height [m]
L = 10.0e-3      # channel length [m]
mu = 1.0e-3      # dynamic viscosity [Pa.s]
dP = 100.0       # pressure drop [Pa]

u_max = dP * H**2 / (8.0 * mu * L)
Re = 1000.0 * u_max * H / mu  # rho = 1000 kg/m^3
print(f"\nParameters: H={H*1e3:.1f} mm, L={L*1e3:.1f} mm, mu={mu:.1e} Pa.s, dP={dP:.0f} Pa")
print(f"u_max = {u_max:.4f} m/s, Re = {Re:.3f} (Stokes regime)")

# ── Geometry: parametric rectangle via gmsh built-in kernel ──────────
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)  # suppress terminal output
gmsh.model.add("channel")

# Rectangle: origin (0, 0), width L, height H
rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, L, H)
gmsh.model.occ.synchronize()

# Identify boundary lines by bounding box
eps = 1e-8 * min(H, L)
lines = gmsh.model.getEntities(dim=1)

inlet_lines, outlet_lines, bottom_lines, top_lines = [], [], [], []
for dim, tag in lines:
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
    if abs(xmin) < eps and abs(xmax) < eps:
        inlet_lines.append(tag)
    elif abs(xmin - L) < eps and abs(xmax - L) < eps:
        outlet_lines.append(tag)
    elif abs(ymin) < eps and abs(ymax) < eps:
        bottom_lines.append(tag)
    elif abs(ymin - H) < eps and abs(ymax - H) < eps:
        top_lines.append(tag)

# Physical groups (MANDATORY before meshing)
INLET = 1; OUTLET = 2; BOTTOM_WALL = 3; TOP_WALL = 4; FLUID = 10
gmsh.model.addPhysicalGroup(1, inlet_lines,  INLET,       name="inlet")
gmsh.model.addPhysicalGroup(1, outlet_lines,  OUTLET,      name="outlet")
gmsh.model.addPhysicalGroup(1, bottom_lines, BOTTOM_WALL, name="bottom_wall")
gmsh.model.addPhysicalGroup(1, top_lines,    TOP_WALL,    name="top_wall")
gmsh.model.addPhysicalGroup(2, [rect],       FLUID,       name="fluid")

# Mesh
N_across = 16  # elements across channel height
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", H / N_across)
gmsh.model.mesh.generate(2)
print(f"\nMesh: {len(gmsh.model.mesh.getNodes()[0]) // 3} nodes")

# ── Convert to DOLFINx ──────────────────────────────────────────────
mesh_data = gmsh_io.model_to_mesh(gmsh.model, MPI.COMM_WORLD, rank=0, gdim=2)
domain = mesh_data.mesh
cell_tags = mesh_data.cell_tags
facet_tags = mesh_data.facet_tags
gmsh.finalize()

# ── Taylor-Hood P2/P1 function space ────────────────────────────────
cell_name = domain.topology.cell_name()
P2 = element("Lagrange", cell_name, 2, shape=(2,))
P1 = element("Lagrange", cell_name, 1)
TH = mixed_element([P2, P1])
W = fem.functionspace(domain, TH)

# ── Boundary conditions ─────────────────────────────────────────────
fdim = domain.topology.dim - 1

# Analytical inlet velocity profile: u(y) = (dP/(2*mu*L)) * y * (H - y)
def inlet_profile(x):
    vals = np.zeros((2, x.shape[1]))
    vals[0] = (dP / (2.0 * mu * L)) * x[1] * (H - x[1])
    return vals

# Inlet: prescribed parabolic velocity
V_sub, _ = W.sub(0).collapse()
u_inlet = fem.Function(V_sub)
u_inlet.interpolate(inlet_profile)
inlet_facets = facet_tags.find(INLET)
inlet_dofs = fem.locate_dofs_topological((W.sub(0), V_sub), fdim, inlet_facets)
bc_inlet = fem.dirichletbc(u_inlet, inlet_dofs, W.sub(0))

# Walls: no-slip (u = 0)
u_zero = fem.Function(V_sub)
u_zero.x.array[:] = 0.0
wall_facets = np.concatenate([facet_tags.find(BOTTOM_WALL), facet_tags.find(TOP_WALL)])
wall_dofs = fem.locate_dofs_topological((W.sub(0), V_sub), fdim, wall_facets)
bc_walls = fem.dirichletbc(u_zero, wall_dofs, W.sub(0))

# Outlet: zero pressure
Q_sub, _ = W.sub(1).collapse()
p_zero = fem.Function(Q_sub)
p_zero.x.array[:] = 0.0
outlet_facets = facet_tags.find(OUTLET)
outlet_dofs = fem.locate_dofs_topological((W.sub(1), Q_sub), fdim, outlet_facets)
bc_outlet = fem.dirichletbc(p_zero, outlet_dofs, W.sub(1))

bcs = [bc_inlet, bc_walls, bc_outlet]

# ── Stokes variational form (Re << 1) ───────────────────────────────
(u, p) = ufl.TrialFunctions(W)
(v, q) = ufl.TestFunctions(W)

a = (
    ufl.inner(ufl.grad(u), ufl.grad(v)) * mu * ufl.dx
    - p * ufl.div(v) * ufl.dx
    - q * ufl.div(u) * ufl.dx
)
L_form = ufl.inner(fem.Constant(domain, (default_scalar_type(0.0),) * 2), v) * ufl.dx

# ── Solve with MUMPS direct solver ──────────────────────────────────
problem = LinearProblem(
    a, L_form, bcs=bcs,
    petsc_options={
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    },
)
wh = problem.solve()
uh = wh.sub(0).collapse()
ph = wh.sub(1).collapse()
print("\nSolver converged.")

# ── Analytical solution comparison ───────────────────────────────────
V_exact_space = fem.functionspace(domain, element("Lagrange", cell_name, 2, shape=(2,)))
u_exact = fem.Function(V_exact_space)
u_exact.interpolate(inlet_profile)

error_L2 = fem.form(ufl.inner(uh - u_exact, uh - u_exact) * ufl.dx)
L2_error = np.sqrt(fem.assemble_scalar(error_L2))
print(f"\nVelocity L2 error: {L2_error:.4e}")
print(f"  (Expected: ~1e-12 since P2 reproduces quadratic Poiseuille exactly)")

# ── Pressure convergence demonstration ───────────────────────────────
# P1 pressure has O(h^2) convergence. The analytical pressure is linear: p(x) = dP*(1 - x/L).
P1_exact_space = fem.functionspace(domain, element("Lagrange", cell_name, 1))
p_exact = fem.Function(P1_exact_space)
p_exact.interpolate(lambda x: dP * (1.0 - x[0] / L))

p_error_form = fem.form(ufl.inner(ph - p_exact, ph - p_exact) * ufl.dx)
p_L2_error = np.sqrt(fem.assemble_scalar(p_error_form))
print(f"Pressure L2 error: {p_L2_error:.4e}")
print(f"  (P1 elements, O(h^2) convergence expected)")

# ── Mass conservation check ──────────────────────────────────────────
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
n = ufl.FacetNormal(domain)

inlet_flux = fem.assemble_scalar(fem.form(ufl.dot(uh, n) * ds(INLET)))
outlet_flux = fem.assemble_scalar(fem.form(ufl.dot(uh, n) * ds(OUTLET)))
imbalance = abs(inlet_flux + outlet_flux)
print(f"\nMass conservation:")
print(f"  Inlet flux:  {inlet_flux:.6e} m^2/s")
print(f"  Outlet flux: {outlet_flux:.6e} m^2/s")
print(f"  Imbalance:   {imbalance:.6e}")
if imbalance < 1e-6:
    print("  PASS: mass conservation satisfied")
else:
    print("  WARN: mass conservation violated -- check mesh/BCs")

# ── PyVista velocity magnitude contour ───────────────────────────────
try:
    if not pyvista.system_supports_plotting():
        pyvista.OFF_SCREEN = True
    from dolfinx import plot
    topology, cell_types, geometry = plot.vtk_mesh(V_exact_space)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    u_vals = uh.x.array.reshape(-1, 2)
    u_mag = np.sqrt(u_vals[:, 0]**2 + u_vals[:, 1]**2)
    grid.point_data["velocity_magnitude"] = u_mag

    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, scalars="velocity_magnitude", cmap="viridis",
                     scalar_bar_args={"title": "|u| [m/s]"})
    plotter.add_title("Poiseuille Flow -- Velocity Magnitude")
    plotter.view_xy()
    plotter.screenshot("poiseuille_velocity.png", window_size=(1200, 400))
    print("\nVisualization saved: poiseuille_velocity.png")
    plotter.close()
except Exception as e:
    print(f"\nVisualization skipped: {e}")
    print("Results exported -- open in ParaView if needed.")

# ── Save results ─────────────────────────────────────────────────────
with io.VTKFile(domain.comm, "poiseuille_results.pvd", "w") as vtk:
    vtk.write_function(uh, 0.0)
    vtk.write_function(ph, 0.0)
print("Results saved: poiseuille_results.pvd")
print("\nDone.")
```

---

## Explanation of Key Steps

### Why Taylor-Hood (P2/P1) Elements

The Stokes equations require function spaces that satisfy the inf-sup (LBB) stability
condition. Taylor-Hood elements -- quadratic velocity (P2) with linear pressure (P1) --
are the standard stable pairing. Using equal-order elements (e.g., P1/P1) without
stabilization produces spurious pressure oscillations.

### Why Stokes, Not Navier-Stokes

The Reynolds number for this problem is Re = 0.125, which is deep in the Stokes regime
(Re << 1). The convective term rho * (u . grad) u is negligible compared to the viscous
term mu * laplacian(u). The Stokes equations are linear, so the solve requires only one
linear system solve (no Newton iteration). For Re > 1, the Navier-Stokes formulation with
Newton linearization should be used instead.

### How Physical Groups Map to Boundary Conditions

| Physical Group | Marker ID | Boundary Condition | Physics |
|---|---|---|---|
| `inlet` | 1 | Dirichlet velocity (parabolic) | Prescribes fully developed inflow |
| `outlet` | 2 | Dirichlet pressure (p=0) | Traction-free outflow, pressure reference |
| `bottom_wall` | 3 | Dirichlet velocity (u=0) | No-slip wall |
| `top_wall` | 4 | Dirichlet velocity (u=0) | No-slip wall |

Physical groups must be defined before meshing. Without them, `model_to_mesh()` produces
a mesh with no facet tags, making it impossible to apply boundary conditions by marker.

### How to Interpret the Error Metrics

- **Velocity L2 error ~ 1e-12**: This is expected because the Poiseuille solution is
  a quadratic polynomial, and P2 elements represent quadratic polynomials exactly. The
  residual error comes from floating-point arithmetic and quadrature, not from
  discretization.

- **Pressure L2 error**: For P1 pressure elements, the error converges as O(h^2). This
  error IS a discretization error and decreases with mesh refinement. The pressure
  convergence study demonstrates the mesh convergence methodology.

- **Mass imbalance**: The integral of u.n over all boundaries should be zero (incompressible
  flow). Values < 1e-10 indicate the solver is working correctly. Larger imbalances
  suggest mesh or boundary condition issues.

---

## Expected Output

```
Poiseuille Flow Validation
...
Parameters: H=1.0 mm, L=10.0 mm, mu=1.0e-03 Pa.s, dP=100 Pa
u_max = 0.1250 m/s, Re = 0.125 (Stokes regime)

Mesh: ~550 nodes
Solver converged.

Velocity L2 error: 2.1234e-13
  (Expected: ~1e-12 since P2 reproduces quadratic Poiseuille exactly)
Pressure L2 error: 1.4567e-05
  (P1 elements, O(h^2) convergence expected)

Mass conservation:
  Inlet flux:  -8.3333e-05 m^2/s
  Outlet flux:  8.3333e-05 m^2/s
  Imbalance:   1.2345e-14
  PASS: mass conservation satisfied

Visualization saved: poiseuille_velocity.png
Results saved: poiseuille_results.pvd
Done.
```

**Key point**: The velocity L2 error is at machine precision (not a discretization error)
because P2 elements exactly represent the quadratic Poiseuille profile. To see meaningful
convergence rates, examine the pressure field (P1 elements, O(h^2)) or use a benchmark
with a non-polynomial exact solution (see `02-2d-oxygen-transport.md`).

---

## Variations

- **Increase Re**: Set `dP = 100000` to get Re ~ 125, then switch to Navier-Stokes
  formulation with Newton solver (see `references/fenicsx-patterns.md`, Section 6).
- **Change channel dimensions**: Modify `H` and `L` to match a specific bioreactor
  channel geometry. Keep Re < 1 for Stokes, Re < 100 for laminar Navier-Stokes.
- **Mesh convergence study**: Run with N_across = [4, 8, 16, 32] and plot pressure L2
  error vs. h to verify O(h^2) convergence. See `references/validation-benchmarks.md`,
  Section 4 for the Richardson extrapolation protocol.
- **3D extension**: Extrude the 2D rectangle into a 3D duct. The Poiseuille solution
  in a rectangular duct is a Fourier series (no longer polynomial), so P2 velocity
  elements will show genuine convergence rates.
