# Validation Benchmarks

Reference benchmark cases with analytical solutions, error metrics, and expected convergence
rates. Use these to validate generated simulation code before running production cases.

---

## Benchmark 1: Poiseuille Flow (2D Channel)

### Problem Description

Fully developed laminar flow between two infinite parallel plates (2D channel flow).
The exact solution is a parabolic velocity profile. This benchmark validates the
Stokes/Navier-Stokes solver, boundary condition application, and pressure gradient handling.

### Analytical Solution

```
u(y) = (dP / (2 * mu * L)) * y * (H - y)
v(y) = 0
p(x) = p_outlet + dP * (L - x) / L
```

where y is the cross-channel coordinate (y=0 at bottom wall, y=H at top wall),
x is the streamwise coordinate, and dP is the total pressure drop over length L.

### Reference Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Channel height | H | 1e-3 | m |
| Channel length | L | 10e-3 | m |
| Dynamic viscosity | mu | 1e-3 | Pa.s |
| Pressure drop | dP | 100 | Pa |
| Density | rho | 1000 | kg/m^3 |

### Derived Quantities

- Maximum velocity: u_max = dP * H^2 / (8 * mu * L) = 0.125 m/s
- Mean velocity: u_mean = (2/3) * u_max = 0.08333 m/s
- Reynolds number: Re = rho * u_max * H / mu = 0.125 (Stokes regime, Re << 1)
- Volume flow rate per unit depth: Q = u_mean * H = 8.333e-5 m^2/s

### Boundary Conditions

- Inlet (x=0): prescribed parabolic velocity profile u(y) (Dirichlet)
- Outlet (x=L): zero pressure, p=0 (Dirichlet on pressure)
- Walls (y=0, y=H): no-slip, u=v=0 (Dirichlet)

### Expected Accuracy with Taylor-Hood (P2/P1) Elements

**IMPORTANT**: The Poiseuille velocity profile is a quadratic polynomial in y.
P2 finite elements represent quadratic polynomials exactly. Therefore:

- **Velocity L2 error**: ~1e-12 (machine precision) on ANY mesh, even the coarsest.
  There is no meaningful velocity convergence to observe -- the solution is exact
  by construction regardless of mesh size.
- **Pressure L2 error**: The pressure field is linear in x. P1 elements represent
  linear functions exactly, so pressure error is also at machine precision.
- The theoretical convergence rates (O(h^3) for P2 velocity, O(h^2) for P1 pressure)
  apply to GENERAL smooth solutions, NOT to this specific benchmark where the exact
  solution is a low-order polynomial.

**For meaningful convergence rate demonstration**, use either:
1. The **pressure field convergence** on a problem with non-linear pressure (not this one).
2. The **1D diffusion-reaction benchmark** below (cosh solution is non-polynomial).

### Pass/Fail Criteria

- Velocity L2 error < 1e-10 (should be at machine precision ~1e-12 to 1e-14)
- Maximum pointwise velocity error < 1e-10
- Mass conservation: integral of u.n over all boundaries < 1e-12

### FEniCSx Code Template

```python
"""Benchmark 1: Poiseuille flow validation (2D channel)."""
import numpy as np
from mpi4py import MPI
import dolfinx
from dolfinx import fem, mesh, io
import ufl
from dolfinx.fem.petsc import LinearProblem

assert dolfinx.__version__.startswith("0.10"), (
    f"Requires FEniCSx v0.10, got {dolfinx.__version__}"
)

# --- Parameters ---
H = 1e-3       # channel height [m]
L = 10e-3      # channel length [m]
mu_val = 1e-3  # dynamic viscosity [Pa.s]
dP = 100.0     # pressure drop [Pa]
nx, ny = 80, 16  # mesh resolution

# --- Analytical solution ---
def u_exact_expr(x):
    """Exact parabolic velocity profile."""
    return np.stack([
        (dP / (2.0 * mu_val * L)) * x[1] * (H - x[1]),
        np.zeros_like(x[0])
    ])

def p_exact_expr(x):
    """Exact linear pressure field."""
    return dP * (L - x[0]) / L

# --- Mesh ---
domain = mesh.create_rectangle(
    MPI.COMM_WORLD,
    [np.array([0.0, 0.0]), np.array([L, H])],
    [nx, ny],
    cell_type=mesh.CellType.triangle,
)

# --- Function spaces (Taylor-Hood P2/P1) ---
P2 = ufl.VectorElement("Lagrange", domain.ufl_cell(), 2)
P1 = ufl.FiniteElement("Lagrange", domain.ufl_cell(), 1)
TH = ufl.MixedElement([P2, P1])
W = fem.functionspace(domain, TH)

# --- Boundary conditions ---
V_sub, _ = W.sub(0).collapse()
Q_sub, _ = W.sub(1).collapse()

# Inlet: parabolic profile
def inlet_marker(x):
    return np.isclose(x[0], 0.0)

inlet_facets = mesh.locate_entities_boundary(domain, 1, inlet_marker)
u_inlet = fem.Function(V_sub)
u_inlet.interpolate(u_exact_expr)
inlet_dofs = fem.locate_dofs_topological((W.sub(0), V_sub), 1, inlet_facets)
bc_inlet = fem.dirichletbc(u_inlet, inlet_dofs, W.sub(0))

# Walls: no-slip
def wall_marker(x):
    return np.logical_or(np.isclose(x[1], 0.0), np.isclose(x[1], H))

wall_facets = mesh.locate_entities_boundary(domain, 1, wall_marker)
u_noslip = fem.Function(V_sub)
u_noslip.x.array[:] = 0.0
wall_dofs = fem.locate_dofs_topological((W.sub(0), V_sub), 1, wall_facets)
bc_walls = fem.dirichletbc(u_noslip, wall_dofs, W.sub(0))

# Outlet: zero pressure
def outlet_marker(x):
    return np.isclose(x[0], L)

outlet_facets = mesh.locate_entities_boundary(domain, 1, outlet_marker)
p_outlet = fem.Function(Q_sub)
p_outlet.x.array[:] = 0.0
outlet_dofs = fem.locate_dofs_topological((W.sub(1), Q_sub), 1, outlet_facets)
bc_outlet = fem.dirichletbc(p_outlet, outlet_dofs, W.sub(1))

bcs = [bc_inlet, bc_walls, bc_outlet]

# --- Variational form (Stokes) ---
(u, p) = ufl.TrialFunctions(W)
(v, q) = ufl.TestFunctions(W)
mu = fem.Constant(domain, mu_val)

a = (ufl.inner(mu * ufl.grad(u), ufl.grad(v)) * ufl.dx
     - p * ufl.div(v) * ufl.dx
     - q * ufl.div(u) * ufl.dx)
L_form = ufl.inner(fem.Constant(domain, (0.0, 0.0)), v) * ufl.dx

# --- Solve ---
problem = LinearProblem(a, L_form, bcs=bcs,
                        petsc_options={"ksp_type": "preonly",
                                       "pc_type": "lu",
                                       "pc_factor_mat_solver_type": "mumps"})
wh = problem.solve()
uh = wh.sub(0).collapse()
ph = wh.sub(1).collapse()

# --- Compute errors ---
u_ex = fem.Function(V_sub)
u_ex.interpolate(u_exact_expr)

error_u = fem.form(ufl.inner(uh - u_ex, uh - u_ex) * ufl.dx)
L2_u = np.sqrt(MPI.COMM_WORLD.allreduce(fem.assemble_scalar(error_u), op=MPI.SUM))

p_ex = fem.Function(Q_sub)
p_ex.interpolate(p_exact_expr)
error_p = fem.form((ph - p_ex)**2 * ufl.dx)
L2_p = np.sqrt(MPI.COMM_WORLD.allreduce(fem.assemble_scalar(error_p), op=MPI.SUM))

print(f"Velocity L2 error: {L2_u:.2e}")
print(f"Pressure L2 error: {L2_p:.2e}")

# P2 elements reproduce quadratic Poiseuille exactly -> machine precision
assert L2_u < 1e-10, f"FAIL: velocity L2 error {L2_u:.2e} > 1e-10"
assert L2_p < 1e-10, f"FAIL: pressure L2 error {L2_p:.2e} > 1e-10"
print("PASS: Poiseuille benchmark within machine precision")
```

---

## Benchmark 2: 1D Diffusion-Reaction

### Problem Description

Steady-state 1D diffusion with first-order reaction: -D * d^2c/dx^2 + k * c = 0.
This benchmark validates the species transport solver and reaction term handling.
The cosh analytical solution is non-polynomial, making this benchmark suitable for
convergence rate studies with P1 and P2 elements.

### Analytical Solution

```
c(x) = c0 * cosh(sqrt(k/D) * (L - x)) / cosh(sqrt(k/D) * L)
```

### Reference Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Diffusivity | D | 3e-9 | m^2/s |
| Reaction rate | k | 1e-3 | 1/s |
| Domain length | L | 1e-3 | m |
| Inlet concentration | c0 | 0.2 | mol/m^3 |

### Derived Quantities

- Thiele modulus: phi = L * sqrt(k/D) = 1e-3 * sqrt(1e-3 / 3e-9) = 577.4
- At this Thiele modulus, the reaction layer is very thin near x=0.
  The solution decays rapidly from c0 at x=0 to near zero in the interior.
- Damkoehler number: Da = k * L^2 / D = 333.3

### Boundary Conditions

- x=0: c = c0 (Dirichlet)
- x=L: dc/dx = 0 (Neumann, zero-flux)

### Expected Convergence Rates

Unlike the Poiseuille benchmark, the cosh solution is transcendental (non-polynomial).
Finite elements cannot represent it exactly, so convergence rates are meaningful:

- **P1 elements**: O(h^2) convergence in L2 norm
- **P2 elements**: O(h^3) convergence in L2 norm

This is the recommended benchmark for demonstrating mesh convergence methodology.

### Pass/Fail Criteria

- L2 error < 1e-4 on 100-element mesh with P1
- L2 error < 1e-6 on 100-element mesh with P2
- Observed convergence rate within 10% of theoretical (1.8-2.2 for P1, 2.7-3.3 for P2)

### FEniCSx Code Template

```python
"""Benchmark 2: 1D diffusion-reaction validation."""
import numpy as np
from mpi4py import MPI
import dolfinx
from dolfinx import fem, mesh
import ufl
from dolfinx.fem.petsc import LinearProblem

assert dolfinx.__version__.startswith("0.10"), (
    f"Requires FEniCSx v0.10, got {dolfinx.__version__}"
)

# --- Parameters ---
D_val = 3e-9   # diffusivity [m^2/s]
k_val = 1e-3   # reaction rate [1/s]
L = 1e-3       # domain length [m]
c0 = 0.2       # inlet concentration [mol/m^3]
nel = 100      # number of elements
degree = 1     # polynomial degree (1 or 2)

lam = np.sqrt(k_val / D_val)  # sqrt(k/D)

# --- Analytical solution ---
def c_exact_expr(x):
    return c0 * np.cosh(lam * (L - x[0])) / np.cosh(lam * L)

# --- Mesh (1D interval) ---
domain = mesh.create_interval(MPI.COMM_WORLD, nel, [0.0, L])

# --- Function space ---
V = fem.functionspace(domain, ("Lagrange", degree))

# --- Boundary conditions ---
def left_marker(x):
    return np.isclose(x[0], 0.0)

left_facets = mesh.locate_entities_boundary(domain, 0, left_marker)
left_dofs = fem.locate_dofs_topological(V, 0, left_facets)
bc = fem.dirichletbc(fem.Constant(domain, c0), left_dofs, V)

# --- Variational form ---
c = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
D = fem.Constant(domain, D_val)
k = fem.Constant(domain, k_val)

a = D * ufl.inner(ufl.grad(c), ufl.grad(v)) * ufl.dx + k * c * v * ufl.dx
L_form = fem.Constant(domain, 0.0) * v * ufl.dx

# --- Solve ---
problem = LinearProblem(a, L_form, bcs=[bc],
                        petsc_options={"ksp_type": "preonly",
                                       "pc_type": "lu"})
ch = problem.solve()

# --- Compute error ---
c_ex = fem.Function(V)
c_ex.interpolate(c_exact_expr)

error_form = fem.form((ch - c_ex)**2 * ufl.dx)
L2_err = np.sqrt(MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(error_form), op=MPI.SUM))

# Normalize by exact solution norm
norm_form = fem.form(c_ex**2 * ufl.dx)
L2_exact = np.sqrt(MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(norm_form), op=MPI.SUM))

print(f"L2 error (absolute): {L2_err:.2e}")
print(f"L2 error (relative): {L2_err / L2_exact:.2e}")

if degree == 1:
    assert L2_err < 1e-4, f"FAIL: P1 L2 error {L2_err:.2e} > 1e-4"
elif degree == 2:
    assert L2_err < 1e-6, f"FAIL: P2 L2 error {L2_err:.2e} > 1e-6"

print(f"PASS: 1D diffusion-reaction benchmark (P{degree})")
```

### Convergence Study (Use This for Rate Demonstration)

```python
"""Mesh convergence study using Benchmark 2 (non-polynomial solution)."""
errors = []
h_values = []
nel_list = [25, 50, 100, 200]

for nel in nel_list:
    domain = mesh.create_interval(MPI.COMM_WORLD, nel, [0.0, L])
    V = fem.functionspace(domain, ("Lagrange", degree))

    # ... (solve as above) ...

    h = L / nel
    h_values.append(h)
    errors.append(L2_err)

# Compute observed convergence rates
for i in range(1, len(errors)):
    rate = np.log(errors[i-1] / errors[i]) / np.log(h_values[i-1] / h_values[i])
    print(f"h={h_values[i]:.1e}  L2={errors[i]:.2e}  rate={rate:.2f}")

# Expected: rate ~ 2.0 for P1, ~ 3.0 for P2
```

---

## Benchmark 3: 2D Convection-Diffusion (SUPG Validation)

### Problem Description

Uniform flow carrying a species through a 2D domain with diffusion. This benchmark
validates the SUPG (Streamline Upwind Petrov-Galerkin) stabilization. Without SUPG,
standard Galerkin FEM produces spurious oscillations when the mesh Peclet number > 1.

### Setup

- Domain: rectangle [0, L] x [0, H], L=10e-3 m, H=1e-3 m
- Velocity: uniform u=(U, 0) with U=1e-3 m/s (prescribed, not solved)
- Diffusivity: D=1e-9 m^2/s
- Mesh Peclet number: Pe_h = U * h / (2*D) (must be > 10 for test to be meaningful)
- Inlet (x=0): c = c_in = 0.2 mol/m^3
- Outlet (x=L): zero-flux (natural BC)
- Walls: zero-flux (natural BC)

### Analytical Solution (1D Reduction)

For uniform flow with constant D and no reaction, the steady-state solution
in the flow direction is c(x) = c_in (constant). Oscillations appear in the
discrete solution near the inlet boundary layer when Pe_h >> 1 and SUPG is absent.

### Pass/Fail Criteria

- With SUPG active: max(c) <= c_in * 1.01 (no overshoots > 1%)
- Without SUPG (for comparison): oscillations visible, max(c) >> c_in
- SUPG stabilization parameter: tau = h / (2 * |u|) * (1 - 1/Pe_h) for Pe_h > 1

### Diagnostic Code Pattern

```python
# After solving with and without SUPG:
c_max_supg = ch_supg.x.array.max()
c_max_gal = ch_galerkin.x.array.max()

print(f"Galerkin max(c):      {c_max_gal:.4e} (expect overshoots)")
print(f"SUPG max(c):          {c_max_supg:.4e} (expect <= {c_in * 1.01:.4e})")

assert c_max_supg <= c_in * 1.01, (
    f"FAIL: SUPG overshoot {c_max_supg:.4e} > {c_in * 1.01:.4e}"
)
print("PASS: SUPG stabilization prevents spurious oscillations")
```

---

## Mesh Convergence Study Protocol

Use this protocol with Benchmark 2 (diffusion-reaction) for convergence rate
demonstration. Do NOT use Poiseuille flow for velocity convergence -- P2 elements
give machine-precision results on any mesh for that quadratic solution.

### Procedure

1. **Select at least 4 refinement levels**: h, h/2, h/4, h/8
   (e.g., nel = 25, 50, 100, 200 for 1D; or nx = 10, 20, 40, 80 for 2D)

2. **For each level**: solve the problem, compute the error norm against the
   analytical solution (or a QoI such as outlet average concentration).

3. **Compute observed convergence rate** between successive levels:
   ```
   rate_i = log(e_{i-1} / e_i) / log(h_{i-1} / h_i)
   ```

4. **Apply Richardson extrapolation** (using the two finest meshes):
   ```
   f_exact ~ f_fine + (f_fine - f_coarse) / (r^p - 1)
   ```
   where r = h_coarse / h_fine and p = observed convergence order.

5. **Report table**:
   | Level | h | DOFs | Error | Rate |
   |-------|---|------|-------|------|

### Expected Rates

| Element | Error norm | Expected rate |
|---------|------------|---------------|
| P1 | L2 | ~2.0 |
| P1 | H1 | ~1.0 |
| P2 | L2 | ~3.0 |
| P2 | H1 | ~2.0 |

### Non-Monotone Convergence Diagnostics

If the convergence rate is not monotonically approaching the theoretical value:

- **Rate < expected**: Check for singularities (sharp corners, re-entrant corners),
  boundary layer resolution, or insufficient asymptotic regime (h not small enough).
- **Rate > expected**: Possible super-convergence at special points, or pre-asymptotic
  regime where higher-order terms dominate.
- **Rate oscillates**: Check mesh quality, verify boundary conditions are applied
  correctly, check for bugs in error computation.
- **Rate ~ 0 (flat error)**: Either the solution is in the FE space (like Poiseuille
  with P2), or there is a bug in the error computation (common: forgetting to
  interpolate the exact solution onto the FE space).

---

## Mass and Species Conservation Check Protocol

Conservation checks are essential for validating that the variational form and
boundary conditions are correctly implemented. Run these after every solve.

### Mass Conservation (Fluid)

The integral of the normal velocity over all boundaries must be zero for
incompressible flow:

```python
"""Mass conservation check for fluid solver."""
n = ufl.FacetNormal(domain)
ds = ufl.Measure("ds", domain=domain)

# Total mass flux over all boundaries
mass_flux_form = fem.form(ufl.dot(uh, n) * ds)
mass_flux = MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(mass_flux_form), op=MPI.SUM)

# Relative to inlet flux
inlet_flux_form = fem.form(ufl.dot(uh, n) * ds(1))  # tag 1 = inlet
inlet_flux = abs(MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(inlet_flux_form), op=MPI.SUM))

rel_error = abs(mass_flux) / max(inlet_flux, 1e-30)
print(f"Mass conservation: net flux = {mass_flux:.2e}, relative = {rel_error:.2e}")
assert rel_error < 1e-6, f"FAIL: mass conservation violation {rel_error:.2e}"
```

### Species Conservation

For steady-state species transport with reaction:
inlet_flux + outlet_flux + total_reaction = 0

```python
"""Species conservation check for transport solver."""
n = ufl.FacetNormal(domain)
ds = ufl.Measure("ds", domain=domain)
dx = ufl.Measure("dx", domain=domain)

# Diffusive + convective flux at inlet
inlet_species = fem.form(
    (ch * ufl.dot(uh, n) - D_val * ufl.dot(ufl.grad(ch), n)) * ds(1))
inlet_flux = MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(inlet_species), op=MPI.SUM)

# Diffusive + convective flux at outlet
outlet_species = fem.form(
    (ch * ufl.dot(uh, n) - D_val * ufl.dot(ufl.grad(ch), n)) * ds(2))
outlet_flux = MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(outlet_species), op=MPI.SUM)

# Total reaction (volumetric sink)
reaction_total = fem.form(k_val * ch * dx)
total_rxn = MPI.COMM_WORLD.allreduce(
    fem.assemble_scalar(reaction_total), op=MPI.SUM)

balance = inlet_flux + outlet_flux + total_rxn
ref = max(abs(inlet_flux), 1e-30)
rel_err = abs(balance) / ref
print(f"Species balance: in={inlet_flux:.2e}, out={outlet_flux:.2e}, "
      f"rxn={total_rxn:.2e}, balance={balance:.2e}, rel={rel_err:.2e}")
assert rel_err < 0.01, f"FAIL: species conservation violation {rel_err:.2e}"
```

### Tolerance Thresholds

| Check | Absolute tolerance | Relative tolerance | Action if violated |
|-------|-------------------|-------------------|-------------------|
| Mass conservation | 1e-10 m^3/s | 1e-6 | Check mesh, BCs, variational form |
| Species conservation | -- | 0.01 (1%) | Refine mesh, check reaction term sign |
| Divergence-free | 1e-8 | -- | Check mixed element pairing |

---

## Error Quantification Methods

### When to Use Which Metric

| Metric | Formula | Use when |
|--------|---------|----------|
| L2 norm | sqrt(integral((u_h - u_exact)^2 dx)) | Overall solution accuracy; convergence studies |
| H1 semi-norm | sqrt(integral(|grad(u_h - u_exact)|^2 dx)) | Gradient accuracy matters (e.g., wall shear stress) |
| Pointwise error | |u_h(x_probe) - u_exact(x_probe)| | Validating at sensor locations; comparing to experiments |
| QoI error | |Q(u_h) - Q(u_exact)| / |Q(u_exact)| | Engineering design decisions (outlet concentration, max velocity) |

### Code Patterns

```python
# L2 error
error_L2 = fem.form(ufl.inner(uh - u_ex, uh - u_ex) * ufl.dx)
L2 = np.sqrt(comm.allreduce(fem.assemble_scalar(error_L2), op=MPI.SUM))

# H1 semi-norm error
error_H1 = fem.form(ufl.inner(ufl.grad(uh - u_ex), ufl.grad(uh - u_ex)) * ufl.dx)
H1 = np.sqrt(comm.allreduce(fem.assemble_scalar(error_H1), op=MPI.SUM))

# Pointwise error (at a probe location)
from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells
tree = bb_tree(domain, domain.topology.dim)
probe = np.array([[5e-3, 0.5e-3, 0.0]])  # midpoint
cell_candidates = compute_collisions_points(tree, probe)
cells = compute_colliding_cells(domain, cell_candidates, probe)
if len(cells.links(0)) > 0:
    val_h = uh.eval(probe, cells.links(0)[0])
    val_ex = u_exact_expr(probe.T)
    print(f"Pointwise error at probe: {abs(val_h - val_ex[0]):.2e}")

# QoI: outlet average concentration
outlet_area_form = fem.form(fem.Constant(domain, 1.0) * ds(2))
outlet_area = comm.allreduce(fem.assemble_scalar(outlet_area_form), op=MPI.SUM)
c_avg_form = fem.form(ch * ds(2))
c_avg = comm.allreduce(fem.assemble_scalar(c_avg_form), op=MPI.SUM) / outlet_area
print(f"Outlet average concentration: {c_avg:.4e} mol/m^3")
```

### Reporting Standard

Always report errors in this format for consistency:

```
=== Validation Report ===
Benchmark: [name]
Mesh: [nel] elements, h = [value] m
Element: P[degree]

L2 error (absolute): X.XXe-XX
L2 error (relative): X.XXe-XX
H1 error (absolute): X.XXe-XX  [if applicable]
Pointwise at probe:  X.XXe-XX  [if applicable]
QoI error:           X.XX%     [if applicable]

Mass conservation:    X.XXe-XX (relative)
Species conservation: X.XXe-XX (relative) [if applicable]

Convergence rate: X.XX (expected: X.XX)
Status: PASS / FAIL
=========================
```
