# Physics Models for Bioprocess Cartridge CFD

Reference file for the `cfd-bioreactor` skill. Contains governing equations,
variational formulations, O2 sink models, boundary condition patterns,
bioreactor parameter lookup tables, and coupling strategy.

All equations target FEniCSx v0.10 variational formulations. All parameter
values are in SI units unless explicitly stated otherwise.

---

## 1. Incompressible Navier-Stokes Equations

### Strong Form

For an incompressible Newtonian fluid with constant density rho and dynamic
viscosity mu:

```
Momentum:   rho (du/dt + u . grad u) = -grad p + mu laplacian u + f
Continuity: div u = 0
```

For steady-state (the primary regime for bioprocess cartridges at low flow
rates):

```
rho (u . grad u) = -grad p + mu laplacian u + f
div u = 0
```

### Reynolds Number Estimation

```
Re = rho * U * L / mu
```

Where:
- U = mean inlet velocity (m/s)
- L = characteristic length, typically channel hydraulic diameter (m)
- rho = fluid density (kg/m^3)
- mu = dynamic viscosity (Pa.s)

**Decision tree**:
- Re < 1: Use Stokes equations (drop convective term, linear problem)
- 1 <= Re < 2000: Use Navier-Stokes with Newton continuation
- Re >= 2000: Turbulent regime -- outside scope of this skill; warn user

Typical bioprocess cartridges operate at Re = 0.001 to 10. Most cases fall
in the Stokes regime.

### Stokes Simplification (Re << 1)

When Re < 1, drop the convective term u . grad u:

```
-grad p + mu laplacian u + f = 0
div u = 0
```

This yields a linear system -- no Newton iteration required.

### Weak Form for FEniCSx

Let V be the velocity test function space and Q the pressure test function space.
Trial functions (u, p), test functions (v, q):

**Stokes weak form**:
```
a((u,p), (v,q)) = mu * inner(grad(u), grad(v)) * dx
                 - inner(p, div(v)) * dx
                 - inner(div(u), q) * dx
L((v,q))        = inner(f, v) * dx
```

**Navier-Stokes weak form** (Newton linearization):
```
F = rho * inner(dot(u, nabla_grad(u)), v) * dx
  + mu * inner(grad(u), grad(v)) * dx
  - inner(p, div(v)) * dx
  - inner(div(u), q) * dx
  - inner(f, v) * dx
```

---

## 2. Species Transport Equations

### Convection-Diffusion-Reaction Equation

Full transient form:

```
dc/dt + u . grad c = D * laplacian c + R(c)
```

Where:
- c = species concentration (mol/m^3)
- u = velocity field (from flow solution)
- D = diffusion coefficient (m^2/s)
- R(c) = reaction source/sink term (mol/(m^3.s))

### Steady-State Form (Primary for This Skill)

For bioprocess cartridges at steady operation:

```
u . grad c = D * laplacian c + R(c)
```

This is the default equation. Transient formulation is only needed for
startup dynamics or periodic operation.

### Peclet Number

```
Pe = |u| * h / (2 * D)
```

Where h is the local cell diameter.

**Decision tree**:
- Pe < 1: Standard Galerkin is acceptable (diffusion-dominated)
- Pe >= 1: SUPG stabilization is REQUIRED (convection-dominated)

Since Pe >= 1 is common in practice and SUPG has negligible cost when Pe < 1,
ALWAYS use SUPG by default. This is not optional.

### SUPG Weak Form (Default -- Always Use)

The SUPG-stabilized weak form modifies the test function:

```
v_supg = v + tau * dot(u, grad(v))
```

Stabilization parameter tau:

```
u_mag = sqrt(dot(u, u) + 1e-10)   # regularized magnitude
Pe_local = u_mag * h / (2 * D)
tau = h / (2 * u_mag) * min(1, Pe_local / 3)
```

The complete SUPG weak form:

```
F = inner(dot(u, grad(c)), v) * dx
  + D * inner(grad(c), grad(v)) * dx
  - inner(R(c), v) * dx
  + inner(dot(u, grad(c)) - D * div(grad(c)) - R(c),
          tau * dot(u, grad(v))) * dx
```

The last line is the SUPG stabilization term (residual-weighted).

### Zero-Velocity Special Case

When the velocity field is zero (pure diffusion, or in stagnant zones):

```
if |u| < eps:
    tau = 0  (SUPG stabilization deactivates gracefully)
```

The regularization `sqrt(dot(u, u) + 1e-10)` in u_mag prevents division by
zero. When u = 0, the SUPG term contributes nothing because tau scales with
h / u_mag while the test function modification scales with dot(u, grad(v)),
which is also zero. The formulation is inherently safe.

---

## 3. O2 Sink Models

Three models are provided, from most physically realistic to simplest.

### 3.1 Regularized Michaelis-Menten (PRIMARY Model)

This is the default O2 consumption model for cell-containing regions.
Standard Michaelis-Menten kinetics with smooth regularization to prevent
numerical issues at c -> 0.

**Formulation**:

```
c_pos = (c + sqrt(c^2 + eps^2)) / 2     # smooth approximation of max(c, 0)
R = -Vmax * c_pos / (Km + c_pos)
```

**Parameters**:
- Vmax: maximum consumption rate (mol/(m^3.s)) -- depends on cell type and density
- Km: Michaelis constant (mol/m^3) -- half-saturation concentration
- eps: regularization parameter

**Regularization parameter selection**:

```
eps = 1e-10 * c_inlet
```

Where c_inlet is the inlet O2 concentration. This scales the regularization
with the problem's concentration magnitude, ensuring:
- eps is small enough to not affect the solution where c >> eps
- eps prevents sqrt(c^2) from evaluating to |c| with discontinuous derivative at c = 0
- The smooth max function c_pos approaches c when c > 0 and approaches 0 when c < 0

**Mathematical justification**: The function `(c + sqrt(c^2 + eps^2)) / 2`
is a C-infinity approximation to max(c, 0). As eps -> 0, it converges to the
exact max function. The derivative is continuous everywhere, which is required
for Newton convergence in the nonlinear solver.

**When to use**: Always, for any Michaelis-Menten kinetics problem. The
regularization has negligible computational cost and prevents solver failures.

**Implementation note for subdomain-specific reactions**:

```python
# Only apply reaction in cell-containing subdomain (marker = 5)
F_reaction = -inner(R_mm, v) * dx(cell_region_marker)
# No reaction elsewhere (R = 0 by omission)
```

### 3.2 Membrane Permeation (Robin Boundary Condition)

For O2 supply through a permeable membrane separating the culture medium from
an external O2 source.

**Formulation (Fick's law of permeation)**:

```
J = P * (c_ext - c)      # flux through membrane (mol/(m^2.s))
```

Where:
- P: membrane permeability (m/s) -- see lookup table below
- c_ext: external O2 concentration (mol/m^3) -- typically atmospheric saturation
- c: local O2 concentration at membrane surface (mol/m^3)

This enters the weak form as a Robin (mixed) boundary condition on the
membrane surface:

```python
# Robin BC term added to the variational form
F_robin = -P * (c_ext - c) * v * ds(membrane_marker)
```

**Biot number estimation**:

```
Bi = P * L / D
```

Where L is the characteristic length (channel height or membrane-to-center
distance) and D is the O2 diffusion coefficient.

**Interpretation**:
- Bi > 1000: membrane resistance is negligible; effectively a Dirichlet BC
  (c = c_ext on membrane). Consider replacing Robin BC with Dirichlet for
  computational simplicity.
- 0.001 < Bi < 1000: membrane resistance matters; use Robin BC as formulated.
- Bi < 0.001: membrane is effectively impermeable; O2 flux through membrane
  is negligible. Consider omitting the membrane BC entirely.

### 3.3 Constant Reaction Rate (Simplest Model)

For quick estimates or when detailed kinetic data is unavailable.

**Formulation**:

```
R = -k * c       # first-order consumption
```

Where k is a first-order rate constant (1/s).

**When to use**: Initial scoping calculations, parameter sensitivity studies,
or when Vmax/Km data is not available for the specific cell type.

**Implementation with subdomain markers**:

```python
# First-order consumption only in cell region
R_first_order = -k * c
F_reaction = -inner(R_first_order, v) * dx(cell_region_marker)
```

### Damkohler Number

The Damkohler number characterizes the ratio of reaction rate to transport rate:

```
Da = Vmax * L / (D * c_inlet)     # for Michaelis-Menten
Da = k * L^2 / D                   # for first-order kinetics
```

**Interpretation**:
- Da << 1: Transport is fast relative to reaction; O2 distribution is nearly uniform
- Da ~ 1: Transport and reaction are comparable; interesting physics regime
- Da >> 1: Reaction is fast relative to transport; sharp concentration gradients
  expected near cell regions. Mesh refinement near reaction zones is critical.

---

## 4. Boundary Condition Patterns

### Inlet Conditions

**Parabolic velocity profile** (fully developed flow in a channel of height H):

```
u_x(y) = 6 * U_mean * y * (H - y) / H^2
u_y = 0
```

Where U_mean is the mean inlet velocity.

For a circular pipe of radius R:

```
u_z(r) = 2 * U_mean * (1 - (r/R)^2)
```

**Uniform velocity** (plug flow):

```
u = (U_inlet, 0, 0)   # or appropriate direction
```

**Species inlet**: Dirichlet BC with specified concentration.

```
c = c_inlet    # typically O2 saturation at 37C
```

### Outlet Conditions

**Zero-pressure (traction-free)**: the standard outlet condition.

```
p = 0   on outlet
```

This is a natural (Neumann) condition for the pressure in the weak form --
no explicit BC code is needed; it is satisfied automatically when no Dirichlet
BC is applied on the outlet for pressure.

**Species outlet**: zero-flux (natural Neumann). No explicit BC needed.

### Wall Conditions

**No-slip** (flow):

```
u = 0   on walls
```

**No-flux** (species): the default Neumann condition. No explicit BC needed.
Physically, this means no O2 passes through solid walls.

### Symmetry Conditions

**Slip condition** (flow):

```
u . n = 0   on symmetry plane
```

Implementation: apply zero normal velocity as a Dirichlet BC on the normal
component.

### Membrane Conditions

**Robin BC** (species): see Section 3.2 above. Applied via weak form
modification, not as a Dirichlet BC.

### BC Validation Rules

Before solving, verify the following:

1. **Pressure reference**: At least one boundary must have a pressure
   condition (Dirichlet pressure or zero-traction). If all boundaries have
   prescribed velocity, the pressure is determined only up to a constant
   and the solver may fail or produce wrong results.

2. **At least one outlet**: The domain must have at least one outflow
   boundary. A closed domain with only inlets has no steady-state solution.

3. **Flow rate compatibility**: For steady-state incompressible flow,
   conservation of mass requires:
   ```
   integral(u . n, inlet) + integral(u . n, outlet) = 0
   ```
   The total volumetric flow rate entering must equal the total leaving.
   (Note: with zero-pressure outlet, this is automatically satisfied.)

4. **Species BC completeness**: At least one boundary or source must provide
   the species. If all boundaries are no-flux and there is no source term,
   the species equation has only the trivial solution c = 0.

---

## 5. Bioreactor Parameter Lookup Tables

All values are in SI units. Ranges reflect literature variability. Use the
"Typical" column as default when user does not specify.

### 5.1 Fluid Properties at 37C

| Medium | Density rho (kg/m^3) | Dynamic Viscosity mu (Pa.s) | Notes |
|--------|---------------------|-----------------------------|-------|
| Water | 993 | 6.92e-4 | Pure water at 37C |
| PBS | 1005 | 7.0e-4 | Phosphate buffered saline |
| DMEM | 1007 | 7.8e-4 | Dulbecco's Modified Eagle Medium |
| DMEM + 10% FBS | 1010 | 8.5e-4 | With fetal bovine serum |
| Hydrogel (0.5% agarose) | 1005 | ~1e-2 (effective) | Non-Newtonian; treat as porous medium |
| Hydrogel (collagen, 3 mg/mL) | 1005 | ~5e-3 (effective) | Use Brinkman model for flow in gel |

**Default**: If medium is unspecified, use DMEM properties (rho = 1007, mu = 7.8e-4).

### 5.2 O2 Transport Properties

| Property | Value | Units | Conditions |
|----------|-------|-------|------------|
| D_O2 in water | 2.8e-9 | m^2/s | 37C |
| D_O2 in DMEM | 2.5e-9 | m^2/s | 37C (Radisic et al., 2006) |
| D_O2 in PBS | 2.7e-9 | m^2/s | 37C |
| D_O2 in agarose (0.5%) | 2.1e-9 | m^2/s | 37C (Avgoustiniatos, 2002) |
| D_O2 in collagen gel | 1.8e-9 | m^2/s | 37C, 3 mg/mL (Demol et al., 2011) |
| D_O2 in tissue (generic) | 1.5e-9 | m^2/s | 37C, estimate for dense cell mass |
| Henry's constant H_O2 | 1.3e-3 | mol/(m^3.Pa) | O2 in water at 37C |
| c_sat (air, 1 atm) | 0.200 | mol/m^3 | O2 saturation at 37C, atmospheric |
| c_sat (pure O2, 1 atm) | 1.0 | mol/m^3 | 100% O2 gas phase at 37C |

**Default**: If D_O2 is unspecified, use D_O2 in DMEM = 2.5e-9 m^2/s.
**Default**: If c_inlet is unspecified, use c_sat(air) = 0.200 mol/m^3.

### 5.3 Michaelis-Menten Kinetics by Cell Type

| Cell Type | Vmax (mol/(m^3.s)) | Km (mol/m^3) | Cell Density (cells/m^3) | Source |
|-----------|-------------------|-------------|------------------------|--------|
| HEK-293 | 3.0e-4 | 5.0e-3 | 1e12 | Jorjani & Bhatt, 2005 |
| CHO | 5.5e-4 | 3.8e-3 | 2e12 | Goudar et al., 2011 |
| T cells (activated) | 8.0e-4 | 4.0e-3 | 5e12 | Wagner et al., 2019 |
| Hepatocytes (primary) | 1.2e-3 | 7.0e-3 | 1e12 | Allen & Bhatia, 2003 |
| MSCs | 2.0e-4 | 4.5e-3 | 5e11 | Zhao et al., 2005 |
| iPSC-derived cardiomyocytes | 1.5e-3 | 6.0e-3 | 1e12 | Correia et al., 2018 |
| Fibroblasts (NIH 3T3) | 1.5e-4 | 3.5e-3 | 1e12 | Casey & Arthur, 2000 |

**Vmax scaling**: Vmax is reported at the reference cell density shown. To
scale for a different cell density n_cells:

```
Vmax_scaled = Vmax_ref * (n_cells / n_ref)
```

**Default**: If cell type is unspecified, use HEK-293 parameters.

**Important**: These are representative literature values. Actual consumption
rates vary significantly with passage number, growth phase (log vs.
stationary), and culture conditions. Users should provide experimentally
measured values when available.

### 5.4 Membrane Permeability

| Material | Permeability P (m/s) | O2 Permeability coeff. (barrer) | Thickness range (um) | Notes |
|----------|---------------------|-------------------------------|---------------------|-------|
| PDMS (Sylgard 184) | 3.0e-3 | 600 | 50-500 | Most common; high O2 permeability |
| Silicone rubber (generic) | 2.5e-3 | 500 | 100-1000 | Similar to PDMS |
| PTFE (Teflon) | 5.0e-4 | 100 | 50-200 | Lower permeability |
| PES (polyethersulfone) | 1.0e-5 | ~2 | 100-300 | Low O2, used for filtration |
| Polycarbonate (track-etched) | 8.0e-5 | ~15 | 10-50 | Porous; permeability depends on pore size |

**Note**: Permeability P listed here is the overall transfer coefficient
(permeation coefficient / membrane thickness). If the user specifies membrane
thickness t and permeation coefficient k_perm separately:

```
P = k_perm / t
```

**Default**: If membrane is unspecified, use PDMS (P = 3.0e-3 m/s).

### 5.5 Dimensionless Numbers -- Estimation and Typical Ranges

| Number | Formula | Typical Range (Bioprocess) | Interpretation |
|--------|---------|---------------------------|----------------|
| Reynolds (Re) | rho * U * L / mu | 0.001 - 10 | Flow regime: Re < 1 -> Stokes |
| Peclet (Pe) | U * L / D | 1 - 10000 | Transport regime: Pe > 1 -> SUPG needed |
| Damkohler (Da) | Vmax * L / (D * c0) | 0.1 - 100 | Reaction vs. transport rate |
| Biot (Bi) | P * L / D | 0.1 - 10000 | Membrane resistance significance |
| Sherwood (Sh) | k_m * L / D | 1 - 100 | Mass transfer coefficient (output) |

**Quick estimation for a typical bioprocess cartridge**:

```
Given: L = 1 mm = 1e-3 m, U = 1e-3 m/s (DMEM), D = 2.5e-9 m^2/s

Re = 1007 * 1e-3 * 1e-3 / 7.8e-4 = 1.3    -> borderline Stokes/N-S
Pe = 1e-3 * 1e-3 / 2.5e-9 = 400            -> SUPG required
Da = 3e-4 * 1e-3 / (2.5e-9 * 0.2) = 600    -> reaction-dominated, refine mesh
Bi = 3e-3 * 1e-3 / 2.5e-9 = 1200           -> PDMS membrane is nearly Dirichlet
```

---

## 6. Non-Dimensionalization Guide

### Reference Scales

Choose characteristic scales from the problem:

| Quantity | Scale | Typical Choice |
|----------|-------|----------------|
| Length | L_ref | Channel height or hydraulic diameter |
| Velocity | U_ref | Mean inlet velocity |
| Pressure | P_ref = mu * U_ref / L_ref | Viscous pressure scale (Stokes) |
| Concentration | c_ref = c_inlet | Inlet species concentration |
| Time | t_ref = L_ref / U_ref | Advective time scale |

### Scaled Equations

Substituting dimensionless variables x* = x/L, u* = u/U, p* = p/P_ref,
c* = c/c_ref:

**Flow (Stokes)**:
```
-grad* p* + laplacian* u* = 0
div* u* = 0
```

**Transport**:
```
u* . grad* c* = (1/Pe) * laplacian* c* + Da * R*(c*)
```

Where R* is the dimensionless reaction term.

### When to Non-Dimensionalize

- **Recommended** for 3D problems where variables span many orders of magnitude
  (pressure in Pa vs. concentration in mol/m^3). Improves matrix conditioning.
- **Not needed** for 2D test problems or when using a direct solver (MUMPS),
  which handles poorly-conditioned systems robustly.
- **Required** when comparing results across different geometries or operating
  conditions for parameter studies.

For this skill, dimensional equations are the default. Non-dimensionalization
is offered as an option when the user encounters convergence issues or when
performing parameter sweeps.

---

## 7. Coupling Strategy

### Segregated Approach (Default)

Solve in two sequential phases:

```
Phase A: Solve flow equations (Stokes or Navier-Stokes) -> obtain u, p
Phase B: Solve species transport using u from Phase A -> obtain c
```

**Advantages**:
- Each phase uses different solver settings optimized for its equation type
- Flow solution can be reused across multiple transport scenarios
- Smaller linear systems (easier to solve, less memory)
- Debugging is straightforward: verify flow first, then transport

**When segregated is appropriate** (covers all cases in this skill's scope):
- Dilute species: O2 concentration does not affect fluid density or viscosity
- Constant fluid properties: rho and mu are independent of c
- No buoyancy coupling: density variations due to species are negligible
- No osmotic effects: species does not drive flow

### When Monolithic Coupling is Needed

**Outside this skill's scope**, but for reference:
- Non-dilute species that change fluid density (e.g., high-concentration sugar solutions)
- Buoyancy-driven flows where density depends on temperature or concentration
- Electrokinetic flows with charge transport coupling

If a user requests any of these, warn that the segregated approach is invalid
and a monolithic formulation or iterative coupling with convergence checks
would be needed.

### Iteration Between Phases

For the problems within this skill's scope, a single pass (Phase A then
Phase B) is sufficient. No iteration between phases is needed because the
species transport does not feed back into the flow equations.

If future versions add temperature-dependent viscosity or density-driven
convection, add a Picard outer loop:

```
while not converged:
    Solve Phase A with current c -> update u
    Solve Phase B with current u -> update c
    Check ||c_new - c_old|| / ||c_new|| < tol
```

This is noted for extensibility but is NOT implemented in v1.0.
