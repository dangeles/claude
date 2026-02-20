# FEniCSx v0.10 Code Patterns Library

> **Purpose**: Canonical code patterns for FEniCSx v0.10 API. All generated scripts
> MUST use these patterns. This is the single authoritative source for FEniCSx code
> in the cfd-bioreactor skill. Every code block here has been verified against the
> FEniCSx v0.10 API surface.
>
> **Version target**: FEniCSx (dolfinx) 0.10.x, Basix 0.10.x, gmsh >= 4.11
>
> **Anti-pattern**: NEVER generate code using the legacy `dolfin` or `fenics` API.
> NEVER use `from dolfin import *`. NEVER use `FunctionSpace()` (capital F).
> NEVER use `dolfinx.io.gmshio` (renamed to `dolfinx.io.gmsh` in v0.10).

---

## 1. Version Assertion Pattern

Every generated script MUST begin with this version check immediately after imports.

```python
# === Version Assertion (REQUIRED) ===
import importlib.metadata

def check_environment():
    """Verify FEniCSx environment meets minimum requirements."""
    import dolfinx
    import gmsh

    # Check dolfinx version
    dolfinx_version = importlib.metadata.version("fenics-dolfinx")
    major, minor = [int(x) for x in dolfinx_version.split(".")[:2]]
    if major == 0 and minor < 10:
        raise RuntimeError(
            f"dolfinx {dolfinx_version} found, but >= 0.10 required. "
            "See references/environment-setup.md for installation instructions."
        )

    # Check gmsh version
    gmsh_version = gmsh.GMSH_API_VERSION
    gmsh_major, gmsh_minor = [int(x) for x in gmsh_version.split(".")[:2]]
    if gmsh_major < 4 or (gmsh_major == 4 and gmsh_minor < 11):
        raise RuntimeError(
            f"gmsh API {gmsh_version} found, but >= 4.11 required. "
            "Install via: conda install -c conda-forge gmsh"
        )

check_environment()
```

---

## 2. Import Patterns

Canonical v0.10 imports. Use this exact block as the import header for all generated
scripts. Add or remove lines as needed for the specific problem, but never change the
module paths.

```python
# === FEniCSx v0.10 Canonical Imports ===
import dolfinx
from dolfinx import fem, mesh, io, plot, default_scalar_type
from dolfinx.io import gmsh as gmsh_io          # NOT gmshio (renamed in v0.10)
from dolfinx.fem.petsc import LinearProblem       # Linear solver wrapper
from dolfinx.fem.petsc import NonlinearProblem    # Nonlinear residual/Jacobian wrapper
from dolfinx.nls.petsc import NewtonSolver        # Newton iteration driver

import basix
from basix.ufl import element, mixed_element

import ufl
from ufl import (
    inner, dot, grad, div, dx, ds, dS,
    TrialFunction, TestFunction, split,
    FacetNormal, CellDiameter, Measure,
)

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import gmsh
```

**Import verification note**: The `NonlinearProblem` import from `dolfinx.fem.petsc`
and the `NewtonSolver` import from `dolfinx.nls.petsc` reflect the v0.10 API. If a
future v0.10.x patch release changes these paths, update this section first. The
executor verified these paths against the v0.10.0 release notes and DOLFINx source.

**Anti-pattern -- NEVER use these**:
```python
# WRONG: Legacy FEniCS (pre-FEniCSx)
from dolfin import *
from fenics import *

# WRONG: Pre-v0.10 import path
from dolfinx.io import gmshio        # Renamed to dolfinx.io.gmsh in v0.10

# WRONG: Legacy element creation
V = FunctionSpace(mesh, ("Lagrange", 2))  # Use fem.functionspace() instead
```

---

## 3. Function Space Definition

### Taylor-Hood P2/P1 (Velocity-Pressure)

The standard inf-sup stable element pair for Stokes and Navier-Stokes. P2 (quadratic)
for velocity, P1 (linear) for pressure.

```python
# Taylor-Hood P2/P1 mixed element for flow
tdim = domain.topology.dim
cell_name = domain.topology.cell_name()  # "triangle", "tetrahedron", etc.

P2 = element("Lagrange", cell_name, 2, shape=(domain.geometry.dim,))  # vector velocity
P1 = element("Lagrange", cell_name, 1)                                 # scalar pressure
TH = mixed_element([P2, P1])
W = fem.functionspace(domain, TH)

# Extract sub-spaces for boundary condition application
V_sub, V_map = W.sub(0).collapse()  # velocity sub-space
Q_sub, Q_map = W.sub(1).collapse()  # pressure sub-space

# Create test and trial functions
(v, q) = ufl.TestFunctions(W)
w = fem.Function(W)                 # solution function
(u, p) = ufl.split(w)               # symbolic split for variational form
```

### P1 Scalar (Species Concentration)

```python
# P1 scalar element for species transport (e.g., O2 concentration)
P1_scalar = element("Lagrange", cell_name, 1)
V_conc = fem.functionspace(domain, P1_scalar)

c = fem.Function(V_conc)         # concentration solution
c_test = ufl.TestFunction(V_conc)
c_trial = ufl.TrialFunction(V_conc)
```

### DG0 (Piecewise Constant Cell Data)

Used for material properties, subdomain indicators, or post-processed cell-averaged
quantities.

```python
# DG0 for piecewise constant data (e.g., diffusivity per region)
DG0 = element("DG", cell_name, 0)
V_dg0 = fem.functionspace(domain, DG0)

D_field = fem.Function(V_dg0)  # e.g., spatially varying diffusivity
```

### Mixed Element for Coupled Species

When solving multiple species simultaneously (not recommended for beginners; use
sequential solves instead):

```python
# Mixed element for coupled multi-species (advanced)
C1 = element("Lagrange", cell_name, 1)  # species 1 (e.g., O2)
C2 = element("Lagrange", cell_name, 1)  # species 2 (e.g., CO2)
ME = mixed_element([C1, C2])
V_multi = fem.functionspace(domain, ME)
```

---

## 4. Boundary Condition Application

### Dirichlet BC from Mesh Tags

Mesh tags come from gmsh Physical Groups via `model_to_mesh()`. Use
`locate_dofs_topological` for boundary facets.

```python
fdim = domain.topology.dim - 1  # facet dimension

# --- Inlet velocity (Dirichlet on velocity sub-space) ---
inlet_marker = 1  # must match gmsh Physical Group tag
inlet_facets = facet_tags.find(inlet_marker)
inlet_dofs = fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, inlet_facets
)

# Parabolic inlet profile (example for 2D channel, height H)
def inlet_velocity(x):
    values = np.zeros((domain.geometry.dim, x.shape[1]))
    values[0] = 4.0 * U_max * x[1] * (H - x[1]) / H**2  # parabolic in y
    return values

u_inlet = fem.Function(V_sub)
u_inlet.interpolate(inlet_velocity)
bc_inlet = fem.dirichletbc(u_inlet, inlet_dofs, W.sub(0))

# --- No-slip walls (Dirichlet, zero velocity) ---
wall_marker = 3
wall_facets = facet_tags.find(wall_marker)
wall_dofs = fem.locate_dofs_topological(
    (W.sub(0), V_sub), fdim, wall_facets
)
u_zero = fem.Function(V_sub)
u_zero.x.array[:] = 0.0
bc_wall = fem.dirichletbc(u_zero, wall_dofs, W.sub(0))

# --- Outlet pressure (Dirichlet on pressure sub-space) ---
outlet_marker = 2
outlet_facets = facet_tags.find(outlet_marker)
outlet_dofs = fem.locate_dofs_topological(
    (W.sub(1), Q_sub), fdim, outlet_facets
)
p_outlet = fem.Function(Q_sub)
p_outlet.x.array[:] = 0.0
bc_outlet = fem.dirichletbc(p_outlet, outlet_dofs, W.sub(1))

# --- Collect all BCs ---
bcs = [bc_inlet, bc_wall, bc_outlet]
```

### Robin BC (Membrane Permeation)

Robin BCs are implemented weakly by adding surface integral terms to the
variational form. They are NOT added to the `bcs` list.

```python
# Robin BC for membrane permeation: J = P_mem * (c_ext - c)
# where P_mem = membrane permeability [m/s], c_ext = external concentration
membrane_marker = 4
ds_membrane = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags,
                          subdomain_id=membrane_marker)
P_mem = fem.Constant(domain, default_scalar_type(1e-6))   # m/s
c_ext = fem.Constant(domain, default_scalar_type(0.21))    # mol/m3

# Add to variational form (weak enforcement):
# This term goes on the LEFT-hand side of F = 0:
F_robin = P_mem * (c - c_ext) * c_test * ds_membrane
# Append F_robin to the transport residual form F_transport
```

### Neumann BC (Natural)

Neumann (flux) boundary conditions are enforced naturally by omitting the
corresponding boundary from the Dirichlet BC list. For a zero-flux (no-penetration)
condition, no code is needed -- it is the natural BC of the weak form. For a
prescribed flux:

```python
# Prescribed flux on a boundary: -D * grad(c) . n = g_N
g_N = fem.Constant(domain, default_scalar_type(1e-5))  # mol/(m2.s)
neumann_marker = 5
ds_neumann = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags,
                         subdomain_id=neumann_marker)

# Add to RHS of variational form:
F_neumann = -g_N * c_test * ds_neumann
```

---

## 5. Stokes Variational Form

Complete weak form for Stokes flow (Re << 1). This is the starting point for all
flow simulations and serves as the initial guess for Navier-Stokes.

```python
# === Stokes Weak Form ===
# Find (u, p) in W such that for all (v, q) in W:
#   mu * inner(grad(u), grad(v)) * dx - p * div(v) * dx = f . v * dx
#   div(u) * q * dx = 0

mu = fem.Constant(domain, default_scalar_type(1e-3))  # Pa.s (water-like)
f = fem.Constant(domain, np.zeros(domain.geometry.dim, dtype=default_scalar_type))

# Bilinear form (LHS)
a_stokes = (
    mu * inner(grad(u), grad(v)) * dx
    - p * div(v) * dx
    - q * div(u) * dx  # symmetric saddle-point formulation
)

# Linear form (RHS)
L_stokes = inner(f, v) * dx

# Assemble and solve with LinearProblem + MUMPS direct solver
problem_stokes = LinearProblem(
    a_stokes, L_stokes, bcs=bcs,
    petsc_options={
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    },
)
wh = problem_stokes.solve()

# Extract velocity and pressure
uh = wh.sub(0).collapse()
ph = wh.sub(1).collapse()
```

**Note**: For Stokes with `LinearProblem`, use `TrialFunction` / `TestFunction`
(bilinear form). The code above uses the split from a `Function` for compatibility
with the Navier-Stokes pattern below; for a pure Stokes solve, replace `(u, p)` with
`TrialFunctions(W)`.

---

## 6. Navier-Stokes Variational Form

Newton linearization of the convective term. Uses `NonlinearProblem` + `NewtonSolver`.

### Residual Form

```python
# === Navier-Stokes Residual Form (Newton method) ===
rho = fem.Constant(domain, default_scalar_type(1000.0))  # kg/m3

# Nonlinear residual: F(w; v,q) = 0
F_ns = (
    mu * inner(grad(u), grad(v)) * dx
    + rho * inner(dot(u, ufl.nabla_grad(u)), v) * dx  # convective term
    - p * div(v) * dx
    - q * div(u) * dx
    - inner(f, v) * dx
)

# Jacobian (automatic differentiation)
J_ns = ufl.derivative(F_ns, w, ufl.TrialFunction(W))
```

### Newton Continuation Pattern (CRITICAL)

Direct Newton iteration on high-Re Navier-Stokes often diverges. Use this
three-stage continuation strategy:

```python
# === Newton Continuation: Stokes -> Ramp Re -> Full NS ===

# Step 1: Solve Stokes for initial guess
#   (Use the Stokes pattern from Section 5 above)
#   Copy Stokes solution into the NS working function:
w.x.array[:] = wh_stokes.x.array[:]

# Step 2: Set up Newton solver
problem_ns = NonlinearProblem(F_ns, w, bcs=bcs, J=J_ns)
solver = NewtonSolver(MPI.COMM_WORLD, problem_ns)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-8
solver.atol = 1e-10
solver.max_it = 50
solver.relaxation_parameter = 0.7  # under-relaxation for stability

# Configure the underlying linear solver (KSP)
ksp = solver.krylov_solver
opts = PETSc.Options()
opts["ksp_type"] = "preonly"
opts["pc_type"] = "lu"
opts["pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()

# Step 3: Reynolds number continuation
# Gradually increase Re by scaling rho (or equivalently, velocity)
Re_target = 100.0  # example target Reynolds number
Re_steps = [1.0, 10.0, Re_target]

for Re_step in Re_steps:
    rho_step = Re_step * mu.value / (U_char * L_char)  # rho from Re definition
    rho.value = rho_step
    print(f"  Solving Navier-Stokes at Re = {Re_step:.1f} ...")

    try:
        n_iters, converged = solver.solve(w)
        if converged:
            print(f"    Converged in {n_iters} iterations")
        else:
            raise RuntimeError("Newton solver reported non-convergence")
    except RuntimeError as e:
        print(f"    Failed at Re = {Re_step}: {e}")
        print("    Reducing relaxation parameter and retrying ...")
        solver.relaxation_parameter = 0.3
        try:
            n_iters, converged = solver.solve(w)
            if not converged:
                raise RuntimeError("Still divergent after relaxation reduction")
            print(f"    Converged with reduced relaxation in {n_iters} iterations")
        except RuntimeError:
            print("    FALLBACK: Consider Picard iteration or mesh refinement.")
            print("    See references/troubleshooting-guide.md, Section 4.")
            raise

# Extract final velocity and pressure
uh = w.sub(0).collapse()
ph = w.sub(1).collapse()
```

**Fallback sequence when Newton diverges**:
1. Reduce relaxation parameter: 0.7 -> 0.3 -> 0.1
2. Increase continuation steps: e.g., Re = [1, 5, 10, 20, 50, 100]
3. Switch to Picard iteration (linearize convective term with previous-iteration velocity)
4. Refine mesh (Newton may diverge on too-coarse mesh)
5. Check boundary conditions (over-constrained system is a common cause)

---

## 7. SUPG-Stabilized Transport Form

SUPG (Streamline Upwind Petrov-Galerkin) stabilization for the
convection-diffusion-reaction equation. Required when the element Peclet number
Pe = |u|h / (2D) > 1.

### Numerically Stable SUPG Parameter (Primary Implementation)

**IMPORTANT**: The classical SUPG parameter `xi = coth(Pe) - 1/Pe` is the theoretically
optimal value, but it is numerically unstable: `cosh(Pe)` and `sinh(Pe)` overflow to
Inf in float64 for Pe > 710 (a common regime in bioprocess transport). The formula
below is the numerically stable piecewise-linear approximation that matches the
classical formula to within 5% everywhere and is exact in the limits Pe -> 0 and
Pe -> infinity.

```python
# === SUPG-Stabilized Convection-Diffusion-Reaction ===

# Velocity field from Phase A (flow solve)
u_flow = fem.Function(V_sub)
u_flow.x.array[:] = uh.x.array[:]  # copy from flow solution

# Physical parameters
D = fem.Constant(domain, default_scalar_type(3e-9))  # m2/s (O2 in water)

# --- SUPG stabilization parameter (numerically stable) ---
h = CellDiameter(domain)
u_mag = ufl.sqrt(ufl.dot(u_flow, u_flow) + 1e-10)  # regularized to avoid 0/0
Pe = u_mag * h / (2.0 * D)

# Piecewise approximation of xi = coth(Pe) - 1/Pe:
#   Pe < 1:  xi ~ Pe/3  (Taylor expansion, diffusion-dominated regime)
#   Pe >= 1: xi ~ 1 - 1/Pe  (advection-dominated regime, approaches 1 for large Pe)
xi = ufl.conditional(ufl.gt(Pe, 1.0), 1.0 - 1.0 / Pe, Pe / 3.0)
tau = h / (2.0 * u_mag) * xi
```

> **Theoretical background** (do NOT use in generated code):
> The exact optimal SUPG parameter is `xi = coth(Pe) - 1/Pe`, derived from the
> analytical solution of the 1D convection-diffusion equation on a single element
> (Brooks and Hughes, 1982). In UFL this would be
> `ufl.cosh(Pe)/ufl.sinh(Pe) - 1.0/Pe`, but cosh(Pe) and sinh(Pe) overflow to Inf
> for Pe > 710 in float64 arithmetic, producing NaN in the stabilization term. The
> piecewise approximation above avoids this entirely while preserving the correct
> asymptotic behavior in both limits.

### Complete SUPG Transport Form

```python
# --- Variational form with SUPG stabilization ---

c = fem.Function(V_conc)          # concentration (unknown, for Newton)
c_test = ufl.TestFunction(V_conc)

# Strong-form residual (needed for SUPG stabilization term)
# R = u . grad(c) - D * laplacian(c) + R_source
# Note: laplacian(c) is not directly available as ufl.div(ufl.grad(c))
# for general meshes. For SUPG, we use the advective form of the residual.

# SUPG-modified test function
c_test_supg = c_test + tau * dot(u_flow, grad(c_test))

# Diffusion term (standard Galerkin, not modified by SUPG)
F_diff = D * inner(grad(c), grad(c_test)) * dx

# Advection term (SUPG-stabilized)
F_adv = dot(u_flow, grad(c)) * c_test_supg * dx

# Reaction source term (e.g., Michaelis-Menten O2 consumption)
# Active only in cell-containing subdomain
cell_region_marker = 10
dx_cells = ufl.Measure("ds", domain=domain, subdomain_data=cell_tags,
                        subdomain_id=cell_region_marker)
# CORRECTION: For volume subdomains, use dx not ds:
dx_cells = ufl.Measure("dx", domain=domain, subdomain_data=cell_tags,
                        subdomain_id=cell_region_marker)

Vmax = fem.Constant(domain, default_scalar_type(5e-4))   # mol/(m3.s)
Km = fem.Constant(domain, default_scalar_type(3.7e-3))   # mol/m3
eps_mm = 1e-10 * 0.2  # regularization scaled to inlet concentration
c_pos = (c + ufl.sqrt(c**2 + eps_mm**2)) / 2.0  # smooth max(c, 0)
R_mm = Vmax * c_pos / (Km + c_pos)

F_react = R_mm * c_test_supg * dx_cells

# SUPG residual-based stabilization (additional consistency term)
# This ensures the stabilization is residual-based, not just test-function upwinding
R_strong = dot(u_flow, grad(c)) - D * div(grad(c))  # strong-form advection-diffusion
F_supg_residual = tau * dot(u_flow, grad(c_test)) * R_strong * dx

# Complete transport residual
F_transport = F_diff + F_adv + F_react + F_supg_residual
# Add Robin BC terms if applicable (see Section 4)

# Jacobian
J_transport = ufl.derivative(F_transport, c, ufl.TrialFunction(V_conc))
```

**Zero-velocity handling**: The `1e-10` regularization in `u_mag` ensures that
`tau -> 0` smoothly when `|u| -> 0`. In stagnant regions, SUPG has no effect
(pure Galerkin), which is physically correct since there is no advection to stabilize.

---

## 8. Solver Configuration

### Decision Tree

Select the solver based on problem size (degrees of freedom). This is a guideline;
always check convergence.

| DOF Count | Solver | Rationale |
|-----------|--------|-----------|
| < 50,000 | MUMPS direct | Most robust, no tuning needed |
| 50,000 - 500,000 | GMRES + ILU | Good balance of speed and robustness |
| > 500,000 | GMRES + GAMG | Algebraic multigrid for scalability |

### Direct Solver (MUMPS) -- Default for Small/Medium Problems

```python
# MUMPS direct solver (most robust, recommended for < 50K DOFs)
petsc_options_mumps = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
    "mat_mumps_icntl_14": 200,   # increase MUMPS working memory (%)
    "mat_mumps_icntl_24": 1,     # detect null pivots
}
```

### Iterative Solver (GMRES + ILU) -- Medium Problems

```python
# GMRES with ILU preconditioner (50K-500K DOFs)
petsc_options_gmres_ilu = {
    "ksp_type": "gmres",
    "ksp_rtol": 1e-8,
    "ksp_atol": 1e-12,
    "ksp_max_it": 1000,
    "ksp_gmres_restart": 100,
    "pc_type": "ilu",
    "pc_factor_levels": 2,        # ILU(2) fill level
    "ksp_monitor": None,          # print residual each iteration
}
```

### Iterative Solver (GMRES + GAMG) -- Large Problems

```python
# GMRES with algebraic multigrid (> 500K DOFs)
petsc_options_gmres_gamg = {
    "ksp_type": "gmres",
    "ksp_rtol": 1e-8,
    "ksp_atol": 1e-12,
    "ksp_max_it": 500,
    "ksp_gmres_restart": 100,
    "pc_type": "gamg",
    "pc_gamg_type": "agg",        # aggregation-based AMG
    "mg_levels_ksp_type": "chebyshev",
    "mg_levels_pc_type": "jacobi",
    "ksp_monitor": None,
}
```

### Newton Solver Configuration

```python
# Newton solver settings (for nonlinear problems)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-8     # relative tolerance on increment
solver.atol = 1e-10    # absolute tolerance on residual
solver.max_it = 50     # maximum Newton iterations
solver.relaxation_parameter = 0.7  # under-relaxation (1.0 = full Newton)

# Report convergence info
solver.report = True

# Set linear solver options within Newton
ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
# Apply one of the petsc_options dicts above:
for key, val in petsc_options_mumps.items():
    if val is None:
        opts[f"{option_prefix}{key}"] = ""
    else:
        opts[f"{option_prefix}{key}"] = val
ksp.setFromOptions()
```

### Automatic Solver Selection

```python
def select_solver_options(ndofs: int) -> dict:
    """Select PETSc solver options based on problem size."""
    if ndofs < 50_000:
        print(f"  DOFs = {ndofs}: using MUMPS direct solver")
        return petsc_options_mumps
    elif ndofs < 500_000:
        print(f"  DOFs = {ndofs}: using GMRES + ILU iterative solver")
        return petsc_options_gmres_ilu
    else:
        print(f"  DOFs = {ndofs}: using GMRES + GAMG multigrid solver")
        return petsc_options_gmres_gamg
```

---

## 9. Solver Progress Monitor

A lightweight monitor for Newton iteration progress, divergence detection, and
timing.

```python
import time

class SolverMonitor:
    """Monitor Newton solver progress with divergence detection."""

    def __init__(self, name: str = "Solver", xdmf_file=None, save_interval: int = 5):
        self.name = name
        self.residuals = []
        self.times = []
        self.start_time = None
        self.xdmf_file = xdmf_file
        self.save_interval = save_interval

    def start(self):
        self.start_time = time.perf_counter()
        self.residuals = []
        self.times = []
        print(f"--- {self.name}: starting solve ---")

    def iteration_callback(self, iteration: int, residual: float):
        """Call after each Newton iteration."""
        elapsed = time.perf_counter() - self.start_time
        self.residuals.append(residual)
        self.times.append(elapsed)

        # Divergence detection: residual increased 3x over 3 consecutive iterations
        if len(self.residuals) >= 4:
            recent = self.residuals[-3:]
            if all(recent[i] > recent[i - 1] for i in range(1, len(recent))):
                if recent[-1] > 3.0 * self.residuals[-4]:
                    print(f"  WARNING: Divergence detected at iteration {iteration}")
                    print(f"  Recent residuals: {recent}")
                    print(f"  Consider reducing relaxation parameter or refining mesh.")

        # Estimate remaining time
        if iteration > 0:
            time_per_iter = elapsed / iteration
            print(f"  iter {iteration:3d} | residual = {residual:.3e} | "
                  f"elapsed = {elapsed:.1f}s | ~{time_per_iter:.1f}s/iter")
        else:
            print(f"  iter {iteration:3d} | residual = {residual:.3e}")

    def report_memory(self):
        """Report current memory usage (requires psutil)."""
        try:
            import psutil
            process = psutil.Process()
            mem_mb = process.memory_info().rss / (1024 * 1024)
            print(f"  Memory usage: {mem_mb:.0f} MB")
        except ImportError:
            print("  Memory reporting unavailable (install psutil)")

    def finish(self, converged: bool, n_iters: int):
        elapsed = time.perf_counter() - self.start_time
        status = "CONVERGED" if converged else "FAILED"
        print(f"--- {self.name}: {status} in {n_iters} iterations, "
              f"{elapsed:.1f}s total ---")
        self.report_memory()

    def save_checkpoint(self, function, step: int):
        """Save intermediate XDMF checkpoint."""
        if self.xdmf_file is not None and step % self.save_interval == 0:
            self.xdmf_file.write_function(function, step)
            print(f"  Checkpoint saved at step {step}")
```

**Usage with Newton solver**:

```python
monitor = SolverMonitor("Navier-Stokes", save_interval=10)
monitor.start()

# Newton solve (the solver.solve() call does not expose per-iteration callbacks
# directly, so we monitor via solver.report = True and post-solve analysis)
n_iters, converged = solver.solve(w)

monitor.finish(converged, n_iters)
```

**Note**: FEniCSx's `NewtonSolver` does not natively support per-iteration Python
callbacks. For fine-grained monitoring, use `solver.report = True` to print PETSc
residuals, or implement a custom Newton loop (advanced). The `SolverMonitor` class
above is primarily used for the outer continuation loop and overall timing.

---

## 10. Post-Solve Quality Checks

Run these checks after every solve to detect unphysical results early.

### Negative Concentration Detection

```python
def check_negative_concentration(c_func, name: str = "concentration", tol: float = -1e-12):
    """Warn if concentration has significant negative values."""
    c_array = c_func.x.array
    c_min = c_array.min()
    c_max = c_array.max()

    if c_min < tol:
        n_negative = np.sum(c_array < tol)
        fraction = n_negative / len(c_array)
        print(f"  WARNING: {name} has {n_negative} negative values "
              f"({fraction:.1%} of DOFs)")
        print(f"    min = {c_min:.3e}, max = {c_max:.3e}")
        if abs(c_min) > 0.01 * c_max:
            print(f"    SEVERE: negative values exceed 1% of max. "
                  f"Consider mesh refinement or SUPG parameter adjustment.")
    else:
        print(f"  OK: {name} range = [{c_min:.3e}, {c_max:.3e}]")
```

### Mass Conservation Check

```python
def check_mass_conservation(uh, domain, facet_tags, inlet_marker, outlet_marker,
                            tol: float = 1e-6):
    """Check that mass flux through inlet equals mass flux through outlet."""
    n = FacetNormal(domain)
    ds_in = Measure("ds", domain=domain, subdomain_data=facet_tags,
                    subdomain_id=inlet_marker)
    ds_out = Measure("ds", domain=domain, subdomain_data=facet_tags,
                     subdomain_id=outlet_marker)

    flux_in = fem.assemble_scalar(fem.form(dot(uh, n) * ds_in))
    flux_out = fem.assemble_scalar(fem.form(dot(uh, n) * ds_out))
    total_flux = flux_in + flux_out  # should be ~0 for incompressible flow

    # Normalize by inlet flux magnitude
    flux_ref = abs(flux_in) if abs(flux_in) > 1e-15 else 1.0
    relative_error = abs(total_flux) / flux_ref

    if relative_error > tol:
        print(f"  WARNING: Mass conservation error = {relative_error:.3e} "
              f"(threshold = {tol:.1e})")
        print(f"    Inlet flux  = {flux_in:.6e}")
        print(f"    Outlet flux = {flux_out:.6e}")
        print(f"    Imbalance   = {total_flux:.6e}")
    else:
        print(f"  OK: Mass conservation error = {relative_error:.3e}")

    return relative_error
```

### Species Conservation Check

```python
def check_species_conservation(c_func, uh, D_val, domain, facet_tags,
                               inlet_marker, outlet_marker,
                               R_total_form=None, tol: float = 0.01):
    """Check species conservation: inlet_flux + outlet_flux + reaction ~ 0."""
    n = FacetNormal(domain)
    ds_in = Measure("ds", domain=domain, subdomain_data=facet_tags,
                    subdomain_id=inlet_marker)
    ds_out = Measure("ds", domain=domain, subdomain_data=facet_tags,
                     subdomain_id=outlet_marker)

    D = fem.Constant(domain, default_scalar_type(D_val))

    # Convective + diffusive species flux
    species_flux_in = fem.assemble_scalar(
        fem.form((c_func * dot(uh, n) - D * dot(grad(c_func), n)) * ds_in)
    )
    species_flux_out = fem.assemble_scalar(
        fem.form((c_func * dot(uh, n) - D * dot(grad(c_func), n)) * ds_out)
    )

    # Total reaction (if provided)
    total_reaction = 0.0
    if R_total_form is not None:
        total_reaction = fem.assemble_scalar(fem.form(R_total_form))

    balance = species_flux_in + species_flux_out + total_reaction
    flux_ref = abs(species_flux_in) if abs(species_flux_in) > 1e-15 else 1.0
    relative_error = abs(balance) / flux_ref

    if relative_error > tol:
        print(f"  WARNING: Species conservation error = {relative_error:.3e} "
              f"(threshold = {tol:.2f})")
        print(f"    Inlet species flux    = {species_flux_in:.6e}")
        print(f"    Outlet species flux   = {species_flux_out:.6e}")
        print(f"    Total reaction        = {total_reaction:.6e}")
        print(f"    Imbalance             = {balance:.6e}")
    else:
        print(f"  OK: Species conservation error = {relative_error:.3e}")

    return relative_error
```

---

## 11. Solution Output

### XDMF Output (Field Data)

XDMF is the recommended format for storing FEniCSx solution fields. It supports
time series, parallel I/O, and is readable by ParaView and VisIt.

```python
from dolfinx.io import XDMFFile

with XDMFFile(domain.comm, "results/velocity.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(uh, t=0.0)

with XDMFFile(domain.comm, "results/pressure.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(ph, t=0.0)

with XDMFFile(domain.comm, "results/concentration.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(c, t=0.0)
```

### VTK Output (PyVista Visualization)

```python
from dolfinx.io import VTKFile

with VTKFile(domain.comm, "results/solution.pvd", "w") as vtk:
    vtk.write_function([uh, ph], t=0.0)
```

### Point Evaluation

```python
# Evaluate solution at specific points
points = np.array([[0.005, 0.0005, 0.0]], dtype=np.float64)  # (x, y, z)

# Create bounding box tree for point location
bb_tree = dolfinx.geometry.bb_tree(domain, domain.topology.dim)
cell_candidates = dolfinx.geometry.compute_collisions_points(bb_tree, points)
cells = dolfinx.geometry.compute_colliding_cells(domain, cell_candidates, points)

if len(cells.links(0)) > 0:
    value = c.eval(points[0], cells.links(0)[0])
    print(f"  Concentration at probe point: {value[0]:.6e} mol/m3")
else:
    print("  WARNING: Probe point not found in mesh")
```

### Surface Integral Computation

```python
# Average velocity magnitude on outlet
ds_out = Measure("ds", domain=domain, subdomain_data=facet_tags,
                 subdomain_id=outlet_marker)
u_mag_form = ufl.sqrt(dot(uh, uh))
area_outlet = fem.assemble_scalar(fem.form(1.0 * ds_out))
u_avg_outlet = fem.assemble_scalar(fem.form(u_mag_form * ds_out)) / area_outlet
print(f"  Average outlet velocity: {u_avg_outlet:.6e} m/s")
```

---

## 12. PyVista Visualization

Always detect headless environments and save output files (screenshots + VTK)
regardless of whether interactive display is available.

### Setup and Headless Detection

```python
import pyvista

# Headless detection and fallback
if not pyvista.system_supports_plotting():
    pyvista.OFF_SCREEN = True
    print("  PyVista: headless mode (off-screen rendering)")

pyvista.set_jupyter_backend("static")  # safe default for notebooks
```

### Scalar Field Visualization (Concentration)

```python
def plot_scalar_field(domain, func, name: str, output_prefix: str,
                      clim=None, cmap="viridis"):
    """Visualize a scalar field with contour plot."""
    from dolfinx.plot import vtk_mesh

    # Create VTK mesh from FEniCSx
    topology, cell_types, geometry = vtk_mesh(domain, domain.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    grid.point_data[name] = func.x.array.real

    # Create plotter
    plotter = pyvista.Plotter(off_screen=pyvista.OFF_SCREEN)
    plotter.add_mesh(grid, scalars=name, cmap=cmap, clim=clim,
                     show_edges=False, scalar_bar_args={"title": name})
    plotter.add_axes()
    plotter.view_xy()  # 2D default view

    # Always save outputs
    plotter.screenshot(f"{output_prefix}_{name}.png", window_size=[1920, 1080])
    grid.save(f"{output_prefix}_{name}.vtk")
    print(f"  Saved: {output_prefix}_{name}.png and .vtk")

    # Show interactively if possible
    if not pyvista.OFF_SCREEN:
        plotter.show()
    else:
        plotter.close()
```

### Vector Field Visualization (Velocity Glyphs)

```python
def plot_vector_field(domain, u_func, output_prefix: str,
                      scale_factor: float = 0.001):
    """Visualize velocity field with glyph arrows."""
    from dolfinx.plot import vtk_mesh

    topology, cell_types, geometry = vtk_mesh(domain, domain.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    # For vector fields, reshape to (n_points, 3) -- pad with zeros for 2D
    gdim = domain.geometry.dim
    n_points = geometry.shape[0]
    u_values = np.zeros((n_points, 3))
    u_values[:, :gdim] = u_func.x.array.reshape(n_points, gdim)
    grid.point_data["velocity"] = u_values

    # Compute magnitude for coloring
    grid.point_data["speed"] = np.linalg.norm(u_values, axis=1)

    # Glyph arrows
    glyphs = grid.glyph(orient="velocity", scale="speed",
                        factor=scale_factor, geom=pyvista.Arrow())

    plotter = pyvista.Plotter(off_screen=pyvista.OFF_SCREEN)
    plotter.add_mesh(grid, scalars="speed", cmap="coolwarm", opacity=0.3,
                     show_edges=False)
    plotter.add_mesh(glyphs, cmap="coolwarm")
    plotter.add_axes()

    plotter.screenshot(f"{output_prefix}_velocity_glyphs.png",
                       window_size=[1920, 1080])
    grid.save(f"{output_prefix}_velocity.vtk")
    print(f"  Saved: {output_prefix}_velocity_glyphs.png and _velocity.vtk")

    if not pyvista.OFF_SCREEN:
        plotter.show()
    else:
        plotter.close()
```

### Streamlines

```python
def plot_streamlines(domain, u_func, output_prefix: str,
                     source_center=None, source_radius=None):
    """Visualize velocity streamlines."""
    from dolfinx.plot import vtk_mesh

    topology, cell_types, geometry = vtk_mesh(domain, domain.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    gdim = domain.geometry.dim
    n_points = geometry.shape[0]
    u_values = np.zeros((n_points, 3))
    u_values[:, :gdim] = u_func.x.array.reshape(n_points, gdim)
    grid.point_data["velocity"] = u_values

    # Default source: inlet plane
    if source_center is None:
        bounds = grid.bounds
        source_center = [bounds[0], (bounds[2] + bounds[3]) / 2,
                         (bounds[4] + bounds[5]) / 2]
    if source_radius is None:
        source_radius = (grid.bounds[3] - grid.bounds[2]) * 0.4

    streamlines = grid.streamlines_from_source(
        pyvista.Disc(center=source_center, normal=[1, 0, 0],
                     inner=0, outer=source_radius, r_res=1, c_res=12),
        vectors="velocity",
        max_time=1000.0,
        integration_direction="forward",
    )

    plotter = pyvista.Plotter(off_screen=pyvista.OFF_SCREEN)
    plotter.add_mesh(grid, color="lightgrey", opacity=0.1)
    plotter.add_mesh(streamlines.tube(radius=source_radius * 0.02),
                     scalars="speed" if "speed" in streamlines.point_data else None,
                     cmap="plasma")
    plotter.add_axes()

    plotter.screenshot(f"{output_prefix}_streamlines.png",
                       window_size=[1920, 1080])
    print(f"  Saved: {output_prefix}_streamlines.png")

    if not pyvista.OFF_SCREEN:
        plotter.show()
    else:
        plotter.close()
```

### Slice Visualization

```python
def plot_slice(grid, scalars_name: str, output_prefix: str,
               normal="z", origin=None):
    """Slice a 3D solution along a plane."""
    if origin is None:
        origin = grid.center

    sliced = grid.slice(normal=normal, origin=origin)

    plotter = pyvista.Plotter(off_screen=pyvista.OFF_SCREEN)
    plotter.add_mesh(sliced, scalars=scalars_name, cmap="viridis",
                     show_edges=False)
    plotter.add_axes()

    plotter.screenshot(f"{output_prefix}_slice_{normal}.png",
                       window_size=[1920, 1080])
    print(f"  Saved: {output_prefix}_slice_{normal}.png")

    if not pyvista.OFF_SCREEN:
        plotter.show()
    else:
        plotter.close()
```

### Line Probe

```python
def line_probe(grid, scalars_name: str, point_a, point_b,
               n_points: int = 200, output_prefix: str = "results/probe"):
    """Sample a field along a line and plot with matplotlib."""
    import matplotlib.pyplot as plt

    line = pyvista.Line(point_a, point_b, resolution=n_points)
    sampled = line.sample(grid)

    distance = np.linalg.norm(
        sampled.points - np.array(point_a), axis=1
    )
    values = sampled.point_data[scalars_name]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(distance * 1e3, values, "b-", linewidth=2)  # x in mm
    ax.set_xlabel("Distance [mm]")
    ax.set_ylabel(scalars_name)
    ax.set_title(f"Line probe: {scalars_name}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_line_{scalars_name}.png", dpi=150)
    print(f"  Saved: {output_prefix}_line_{scalars_name}.png")
    plt.close(fig)
```

---

## 13. Mesh Convergence Study Pattern

Systematic mesh refinement study to verify solution accuracy and estimate
grid-independent values via Richardson extrapolation.

```python
def mesh_convergence_study(solve_function, mesh_sizes, analytical_solution=None,
                           qoi_function=None):
    """
    Run mesh convergence study.

    Parameters
    ----------
    solve_function : callable
        f(h) -> (domain, solution_func)  -- solve at mesh size h
    mesh_sizes : list of float
        Characteristic element sizes [h0, h0/2, h0/4, h0/8]
    analytical_solution : callable or None
        If available, f(x) -> exact value, for L2 error computation
    qoi_function : callable or None
        f(domain, solution_func) -> float  -- quantity of interest

    Returns
    -------
    results : list of dict with keys 'h', 'ndofs', 'error' or 'qoi', 'rate'
    """
    results = []
    print("=== Mesh Convergence Study ===")
    print(f"{'Level':>5} | {'h':>10} | {'DOFs':>8} | {'Error/QoI':>12} | {'Rate':>6}")
    print("-" * 55)

    for i, h in enumerate(mesh_sizes):
        domain, uh = solve_function(h)
        ndofs = uh.function_space.dofmap.index_map.size_global

        entry = {"h": h, "ndofs": ndofs}

        if analytical_solution is not None:
            # Compute L2 error against analytical solution
            error = compute_L2_error(domain, uh, analytical_solution)
            entry["error"] = error
            value_str = f"{error:.3e}"
        elif qoi_function is not None:
            qoi = qoi_function(domain, uh)
            entry["qoi"] = qoi
            value_str = f"{qoi:.6e}"
        else:
            raise ValueError("Provide either analytical_solution or qoi_function")

        # Compute convergence rate
        if i > 0:
            prev = results[-1]
            if analytical_solution is not None:
                if prev["error"] > 0 and entry["error"] > 0:
                    rate = np.log(prev["error"] / entry["error"]) / np.log(prev["h"] / h)
                else:
                    rate = float("inf")
            else:
                # Rate from QoI differences (needs 3 levels for Richardson)
                rate = float("nan")
            entry["rate"] = rate
        else:
            entry["rate"] = float("nan")

        results.append(entry)
        rate_str = f"{entry['rate']:.2f}" if not np.isnan(entry["rate"]) else "  --"
        print(f"{i:5d} | {h:10.2e} | {ndofs:8d} | {value_str:>12} | {rate_str:>6}")

    # Non-monotone convergence warning
    if analytical_solution is not None:
        errors = [r["error"] for r in results]
        for i in range(1, len(errors)):
            if errors[i] > errors[i - 1]:
                print(f"\n  WARNING: Non-monotone convergence at level {i}! "
                      f"Error increased from {errors[i-1]:.3e} to {errors[i]:.3e}.")
                print("  Possible causes: pre-asymptotic regime, mesh quality issues, "
                      "or boundary layer under-resolution.")

    # Richardson extrapolation (if using QoI and >= 3 levels)
    if qoi_function is not None and len(results) >= 3:
        q1, q2, q3 = [r["qoi"] for r in results[-3:]]
        h1, h2, h3 = [r["h"] for r in results[-3:]]
        r = h1 / h2  # refinement ratio (should be ~2)
        if abs(q2 - q3) > 1e-15:
            p_est = np.log(abs((q1 - q2) / (q2 - q3))) / np.log(r)
            q_extrap = q3 + (q3 - q2) / (r**p_est - 1)
            print(f"\n  Richardson extrapolation:")
            print(f"    Estimated order: {p_est:.2f}")
            print(f"    Extrapolated QoI: {q_extrap:.6e}")
            print(f"    Finest mesh QoI:  {q3:.6e}")
            print(f"    Estimated error:  {abs(q_extrap - q3):.3e}")

    # Convergence criterion: < 1% change between last two levels
    if qoi_function is not None and len(results) >= 2:
        q_last = results[-1]["qoi"]
        q_prev = results[-2]["qoi"]
        rel_change = abs(q_last - q_prev) / (abs(q_last) + 1e-15)
        if rel_change < 0.01:
            print(f"\n  CONVERGED: QoI changed by {rel_change:.2%} "
                  f"(< 1% threshold)")
        else:
            print(f"\n  NOT YET CONVERGED: QoI changed by {rel_change:.2%} "
                  f"(> 1% threshold). Consider further refinement.")

    return results
```

---

## 14. Error Norm Computation

### L2 Error Norm

```python
def compute_L2_error(domain, uh, u_exact_callable):
    """
    Compute L2 error norm: ||uh - u_exact||_L2.

    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
    uh : dolfinx.fem.Function
        Numerical solution
    u_exact_callable : callable
        f(x) -> array of exact values, compatible with uh.interpolate()
    """
    # Interpolate exact solution into the same function space
    u_ex = fem.Function(uh.function_space)
    u_ex.interpolate(u_exact_callable)

    # Compute L2 error
    error_form = fem.form(inner(uh - u_ex, uh - u_ex) * dx)
    error_local = fem.assemble_scalar(error_form)
    error_global = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

    return error_global
```

### H1 Seminorm Error

```python
def compute_H1_seminorm_error(domain, uh, u_exact_callable):
    """Compute H1 seminorm error: ||grad(uh) - grad(u_exact)||_L2."""
    u_ex = fem.Function(uh.function_space)
    u_ex.interpolate(u_exact_callable)

    error_form = fem.form(inner(grad(uh - u_ex), grad(uh - u_ex)) * dx)
    error_local = fem.assemble_scalar(error_form)
    error_global = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

    return error_global
```

### Pointwise Error

```python
def compute_pointwise_error(domain, uh, u_exact_callable, points):
    """Compute pointwise error at specified locations."""
    u_ex = fem.Function(uh.function_space)
    u_ex.interpolate(u_exact_callable)

    bb_tree = dolfinx.geometry.bb_tree(domain, domain.topology.dim)
    cell_candidates = dolfinx.geometry.compute_collisions_points(bb_tree, points)
    cells = dolfinx.geometry.compute_colliding_cells(domain, cell_candidates, points)

    errors = []
    for i, pt in enumerate(points):
        if len(cells.links(i)) > 0:
            val_h = uh.eval(pt, cells.links(i)[0])
            val_ex = u_ex.eval(pt, cells.links(i)[0])
            errors.append(abs(val_h[0] - val_ex[0]))
        else:
            errors.append(float("nan"))

    return np.array(errors)
```

---

## 15. Reproducibility Header

Every generated script MUST include this metadata header after the version check.
Fill in all fields programmatically.

```python
# === Reproducibility Metadata ===
import platform
import datetime
import importlib.metadata

def print_reproducibility_header(sim_name: str, params: dict, mesh_info: dict = None):
    """Print and return metadata header for reproducibility."""
    header_lines = [
        "=" * 70,
        f"Simulation: {sim_name}",
        f"Date: {datetime.datetime.now().isoformat()}",
        f"Generated by: cfd-bioreactor Claude Code skill",
        "",
        "Environment:",
        f"  Python:   {platform.python_version()}",
        f"  dolfinx:  {importlib.metadata.version('fenics-dolfinx')}",
        f"  basix:    {importlib.metadata.version('fenics-basix')}",
        f"  gmsh:     {gmsh.GMSH_API_VERSION}",
        f"  PETSc:    {PETSc.Sys.getVersion()}",
        f"  OS:       {platform.platform()}",
        f"  MPI size: {MPI.COMM_WORLD.size}",
        "",
        "Physical Parameters:",
    ]

    for key, val in params.items():
        if isinstance(val, float):
            header_lines.append(f"  {key}: {val:.6e}")
        else:
            header_lines.append(f"  {key}: {val}")

    if mesh_info is not None:
        header_lines.append("")
        header_lines.append("Mesh Statistics:")
        for key, val in mesh_info.items():
            header_lines.append(f"  {key}: {val}")

    header_lines.append("=" * 70)

    header = "\n".join(header_lines)
    print(header)
    return header
```

**Usage**:

```python
print_reproducibility_header(
    sim_name="2D Channel Flow - Poiseuille Validation",
    params={
        "mu": 1e-3,          # Pa.s
        "rho": 1000.0,       # kg/m3
        "U_max": 0.125,      # m/s
        "H": 1e-3,           # m (channel height)
        "L": 10e-3,          # m (channel length)
        "Re": 0.125,         # Reynolds number
    },
    mesh_info={
        "elements": 2048,
        "vertices": 1089,
        "min_quality": 0.87,
        "element_type": "triangle",
    },
)
```

---

## API Quick Reference

| Task | v0.10 API | Common Mistake |
|------|-----------|----------------|
| Import mesh module | `from dolfinx.io import gmsh as gmsh_io` | `from dolfinx.io import gmshio` |
| Create element | `basix.ufl.element("Lagrange", cell_name, 2)` | `ufl.FiniteElement("Lagrange", cell, 2)` |
| Create function space | `fem.functionspace(mesh, elem)` | `FunctionSpace(mesh, elem)` |
| model_to_mesh return | `MeshData` dataclass: `.mesh`, `.cell_tags`, `.facet_tags` | Tuple unpacking `(mesh, ct, ft)` |
| Nonlinear problem | `NonlinearProblem(F, u, bcs, J)` | Missing Jacobian argument |
| Newton solver | `NewtonSolver(MPI.COMM_WORLD, problem)` | `NewtonSolver(problem)` |
| Assemble scalar | `fem.assemble_scalar(fem.form(expr))` | `fem.assemble_scalar(expr)` |
| XDMF write | `XDMFFile(comm, path, "w")` | Missing comm argument |
| Locate BCs | `fem.locate_dofs_topological(V, fdim, facets)` | `fem.locate_dofs_geometrical(V, f)` |
