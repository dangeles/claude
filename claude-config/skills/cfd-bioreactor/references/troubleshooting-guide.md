# Troubleshooting Guide

Error catalog organized by workflow stage, with causes, diagnostic commands, and fixes.
When a simulation fails, identify the stage and look up the error message below.

---

## Stage 0: Environment and Installation Errors

### ImportError: cannot import name 'mesh' from 'dolfinx'

**Cause**: FEniCSx version mismatch. The v0.10 API differs significantly from v0.7/v0.8.
**Diagnostic**:
```python
import dolfinx
print(dolfinx.__version__)
```
**Fix**: Recreate the environment from the pinned `environment.yml`:
```bash
conda deactivate
conda env remove -n cfd-bioreactor
conda env create -f environment.yml
conda activate cfd-bioreactor
```

### AttributeError: module 'gmsh.model' has no attribute 'occ'

**Cause**: The gmsh package was installed without OpenCASCADE (OCC) kernel support.
This happens when installing gmsh via pip instead of conda-forge.
**Fix**: Install gmsh from conda-forge (which includes OCC):
```bash
conda install -c conda-forge gmsh python-gmsh
```
**Verify**:
```python
import gmsh
gmsh.initialize()
gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
print("OCC kernel: OK")
gmsh.finalize()
```

### ModuleNotFoundError: No module named 'petsc4py'

**Cause**: PETSc Python bindings not installed or version mismatch with dolfinx.
**Fix**:
```bash
conda install -c conda-forge petsc4py
```
**Note**: petsc4py version must match the PETSc version used to build dolfinx.
Using conda-forge for both ensures compatibility.

### RuntimeError: Cannot create display / No module named 'vtkmodules'

**Cause**: Running on a headless server (no X11 display) or VTK not installed.
**Fix** (headless):
```python
import pyvista
pyvista.OFF_SCREEN = True
pyvista.start_xvfb()  # Linux only; requires xvfb package
```
**Fix** (missing VTK):
```bash
conda install -c conda-forge pyvista vtk
```

### Segmentation fault on import

**Cause**: Version mismatch between PETSc, petsc4py, and dolfinx compiled libraries.
This is the hardest error to debug because there is no Python traceback.
**Fix**: The only reliable fix is to recreate the conda environment from scratch:
```bash
conda deactivate
conda env remove -n cfd-bioreactor
conda env create -f environment.yml
conda activate cfd-bioreactor
```
**Prevention**: Never mix conda-forge and pip for core packages (dolfinx, petsc4py, gmsh).

### ImportError: cannot import name 'gmshio' from 'dolfinx.io'

**Cause**: In FEniCSx v0.10, `dolfinx.io.gmshio` was renamed to `dolfinx.io.gmsh`.
**Fix**: Use the v0.10 import path:
```python
# v0.10 correct:
from dolfinx.io import gmsh as gmsh_io
mesh_data = gmsh_io.model_to_mesh(
    gmsh.model, MPI.COMM_WORLD, 0, gdim=3
)
domain = mesh_data.mesh
cell_tags = mesh_data.cell_tags
facet_tags = mesh_data.facet_tags
```
If `dolfinx.io.gmsh` is not found, verify dolfinx version is 0.10.x.

---

## Stage 1: Geometry and Mesh Errors

### OCC exception: BRep_Tool / Standard_NullObject

**Cause**: Corrupted or incompatible STEP/IGES file. The OCC kernel cannot parse the
boundary representation data.
**Diagnostic**:
```python
try:
    gmsh.model.occ.importShapes("geometry.step")
    gmsh.model.occ.synchronize()
except Exception as e:
    print(f"STEP import failed: {e}")
```
**Fix**:
1. Re-export the file as STEP AP214 from the CAD tool (not IGES, not STEP AP203).
2. If re-export is not possible, try healing the geometry:
   ```python
   gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
   gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
   gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)
   gmsh.model.occ.importShapes("geometry.step")
   gmsh.model.occ.synchronize()
   ```
3. As last resort, simplify the geometry in CAD and re-export.

### Meshing error with no detail / gmsh crashes silently

**Cause**: Small geometric features (fillets, chamfers, thin walls) cause the mesher
to fail when it cannot place elements.
**Diagnostic**:
```python
gmsh.option.setNumber("General.Verbosity", 99)  # maximum verbosity
gmsh.model.mesh.generate(3)
```
**Fix**:
1. Defeature: remove fillets < 0.5 mm, small holes, thin features.
   ```python
   gmsh.option.setNumber("Geometry.Tolerance", 1e-4)
   gmsh.option.setNumber("Mesh.MeshSizeMin", 5e-5)  # set floor on element size
   ```
2. Increase minimum element size until meshing succeeds.
3. If the geometry is too complex, create a simplified fluid domain manually.

### 0 volumes after STEP import

**Cause**: The STEP file contains only 2D geometry (surfaces, no solids).
**Diagnostic**:
```python
gmsh.model.occ.importShapes("geometry.step")
gmsh.model.occ.synchronize()
volumes = gmsh.model.getEntities(3)
surfaces = gmsh.model.getEntities(2)
print(f"Volumes: {len(volumes)}, Surfaces: {len(surfaces)}")
```
**Fix**:
1. If 2D simulation is acceptable, mesh the surfaces instead:
   ```python
   gmsh.model.mesh.generate(2)  # 2D mesh instead of 3D
   ```
2. If 3D is required, extrude a 2D cross-section:
   ```python
   gmsh.model.occ.extrude(surfaces, 0, 0, depth)
   gmsh.model.occ.synchronize()
   ```
3. Re-export from CAD as a solid body (not sheet body).

### Multiple volumes in STEP file

**Cause**: The STEP file is a multi-body assembly (housing + cartridge + seals).
Only the fluid domain volume should be meshed.
**Diagnostic**:
```python
volumes = gmsh.model.getEntities(3)
for v in volumes:
    mass = gmsh.model.occ.getMass(v[0], v[1])
    print(f"Volume {v[1]}: mass proxy = {mass:.6e}")
```
**Fix**: Select only the fluid domain volume by tag:
```python
fluid_volume_tag = 2  # identify by inspection or mass proxy
# Remove non-fluid volumes
for v in volumes:
    if v[1] != fluid_volume_tag:
        gmsh.model.occ.remove([v])
gmsh.model.occ.synchronize()
```

### Mesh quality warning (min quality < 0.01)

**Cause**: Degenerate elements near sharp features or thin regions.
**Diagnostic**:
```python
gmsh.plugin.setNumber("AnalyseMeshQuality", "JacobianDeterminant", 1)
gmsh.plugin.run("AnalyseMeshQuality")
```
**Fix**:
1. Increase `Mesh.MeshSizeMin` near the offending region.
2. Defeature sharp corners or narrow channels.
3. Use `Mesh.OptimizeNetgen = 1` for 3D mesh optimization.

### Out of memory during meshing

**Cause**: Element count is too high for available RAM. A 3D mesh with 10M elements
requires approximately 8-12 GB of RAM.
**Diagnostic**: Estimate element count before meshing:
```python
# Rough estimate: volume / (h^3 * 6) for tetrahedra
import psutil
mem_gb = psutil.virtual_memory().total / 1e9
max_elements = mem_gb * 1e6  # ~1M elements per GB as rule of thumb
```
**Fix**:
1. Increase `Mesh.MeshSizeMin` and `Mesh.MeshSizeMax`.
2. Use mesh size fields to refine only near walls and regions of interest.
3. For very large geometries, consider 2D axisymmetric simulation instead of full 3D.

---

## Stage 2-3: Physics Setup and Solver Configuration Errors

### TypeError in boundary condition application

**Cause**: Wrong function space sub-component index when applying BCs to mixed spaces.
**Example error**: `TypeError: Cannot apply boundary condition to component 2 of space with 2 sub-components`
**Fix**: For Taylor-Hood (P2/P1) mixed space:
```python
W = fem.functionspace(domain, TH)  # mixed space
# W.sub(0) = velocity (vector, P2)
# W.sub(1) = pressure (scalar, P1)
# W.sub(0).sub(0) = velocity x-component
# W.sub(0).sub(1) = velocity y-component
# W.sub(0).sub(2) = velocity z-component (3D only)
```
Check that the sub-component index matches the physical dimension.

### ValueError: incompatible shapes

**Cause**: Element dimension mismatch between trial/test functions and forms.
Common when mixing 2D and 3D elements or scalar and vector quantities.
**Fix**:
1. Verify that the mesh geometric dimension (`gdim`) matches the vector element dimension.
2. For 2D simulations, use `VectorElement("Lagrange", cell, 2, dim=2)`.
3. Check that `ufl.inner`, `ufl.dot`, and `ufl.grad` are used with compatible operands.

### Over-constrained system (all Dirichlet, no pressure reference)

**Cause**: Applying Dirichlet velocity BCs on all boundaries without a pressure
reference point makes the pressure indeterminate (unique only up to a constant).
The solver may fail or produce an arbitrary pressure offset.
**Symptom**: MUMPS reports "numerical pivot too small" or solver diverges.
**Fix**: Add a pressure pin at one node:
```python
# Pin pressure at one point (e.g., outlet corner)
zero_pressure = fem.Constant(domain, 0.0)
dof = fem.locate_dofs_geometrical(Q_sub, lambda x: np.logical_and(
    np.isclose(x[0], L), np.isclose(x[1], 0.0)))
bc_pressure = fem.dirichletbc(zero_pressure, dof[:1], W.sub(1))
```

### No outlet boundary defined

**Cause**: All boundaries have velocity (Dirichlet) conditions but no outlet.
Incompressible flow requires at least one boundary with a pressure/traction condition
to allow fluid to exit.
**Fix**: Designate an outlet with zero-pressure or zero-traction BC. In the mixed
formulation, a zero-pressure outlet is a Dirichlet condition on the pressure sub-space.

---

## Stage 4: Solver Execution Errors

### Newton solver did not converge

This is the most common runtime error. Use the following diagnostic procedure.

**Step 1: Check residual history**
```python
from dolfinx.nls.petsc import NewtonSolver
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-6
solver.atol = 1e-8
solver.max_it = 50
solver.report = True  # prints residual at each iteration

log = dolfinx.log
log.set_log_level(log.LogLevel.INFO)

n_iters, converged = solver.solve(uh)
if not converged:
    print(f"DIVERGED after {n_iters} iterations")
```

**Step 2: Identify the cause**

| Residual pattern | Likely cause | Fix |
|-----------------|--------------|-----|
| Increasing monotonically | Initial guess too far from solution | Use Stokes solution as initial guess |
| Oscillating without decreasing | High Re without continuation | Ramp Re gradually (Re=1, 10, 100, ...) |
| Decreasing then stalling | Stiff reaction terms | Reduce Vmax or increase Km in Michaelis-Menten |
| First iteration already huge | Mesh too coarse for BCs | Refine mesh near walls and inlet |

**Step 3: Apply fixes in sequence**

1. **Reduce relaxation factor** (damped Newton):
   ```python
   solver.relaxation_parameter = 0.5  # try 0.5, 0.3, 0.1
   ```

2. **Use Stokes as initial guess** (removes nonlinearity from Re):
   ```python
   # Solve Stokes (linear) first
   # ... (see Benchmark 1 code) ...
   # Copy Stokes solution as initial guess for Navier-Stokes
   uh_ns.x.array[:] = uh_stokes.x.array[:]
   ```

3. **Reynolds continuation** (gradually increase Re):
   ```python
   for Re_target in [1, 10, 50, 100, Re_final]:
       mu_current = rho * U * L_char / Re_target
       mu_const.value = mu_current
       solver.solve(uh)
       print(f"Re = {Re_target}: converged")
   ```

4. **Picard iteration** (linearize convective term):
   ```python
   # Instead of Newton on full Navier-Stokes:
   # Iterate: solve linearized system using previous velocity for convection
   for picard_iter in range(20):
       # Use uh_prev for the convective velocity
       # Solve linear system for uh_new
       # Check ||uh_new - uh_prev|| < tol
       pass
   ```

5. **Mesh refinement**: If all above fail, the mesh is likely too coarse near
   boundary layers or geometric features. Refine and retry.

### MUMPS: numerical pivot too small / zero pivot

**Cause**: The assembled matrix is singular or nearly singular.
**Common reasons**:
1. Missing pressure reference (all-Dirichlet velocity BCs). See Stage 2-3 fix.
2. Degenerate mesh elements (zero-volume tetrahedra).
3. Incompatible boundary conditions (conflicting Dirichlet values at corners).
**Diagnostic**:
```python
# Check for degenerate elements
from dolfinx.mesh import compute_midpoints
cells = domain.topology.index_map(domain.topology.dim)
# If the matrix condition number is needed:
# Use SLEPc or estimate via PETSc KSP
```
**Fix**: Address the root cause (BC, mesh, or formulation).

### PETSc out of memory (malloc failed)

**Cause**: Direct solver (MUMPS/SuperLU) requires O(N^2) memory for the factored matrix.
For problems > 500K DOFs, direct solvers may exhaust available RAM.
**Diagnostic**:
```python
import psutil
mem = psutil.virtual_memory()
print(f"Available: {mem.available / 1e9:.1f} GB")
ndofs = W.dofmap.index_map.size_global
print(f"DOFs: {ndofs:,}")
# Rule of thumb: direct solver needs ~40 bytes * ndofs^1.5 for 3D
```
**Fix**: Switch to an iterative solver with preconditioning:
```python
petsc_options = {
    "ksp_type": "gmres",
    "ksp_rtol": 1e-8,
    "pc_type": "hypre",
    "pc_hypre_type": "boomeramg",
    "ksp_max_it": 500,
}
```

### KSP did not converge (iterative solver stall)

**Cause**: Poor conditioning of the system matrix, inadequate preconditioner.
**Diagnostic**:
```python
# Add monitoring
petsc_options["ksp_monitor"] = ""  # prints residual per iteration
petsc_options["ksp_converged_reason"] = ""
```
**Fix** (in order of preference):
1. Increase `ksp_max_it` (maybe it just needs more iterations).
2. Try a stronger preconditioner: `"pc_type": "ilu"` or `"pc_type": "lu"` (if memory allows).
3. Scale the equations (non-dimensionalize).
4. Check mesh quality -- poorly shaped elements degrade conditioning.

### NaN in solution

**Cause**: Numerical blow-up, typically from one of:
1. Degenerate mesh elements producing singular element matrices.
2. Boundary conditions not properly applied (missing BC on a boundary).
3. Division by zero in a source term or material property.
**Diagnostic**:
```python
import numpy as np
if np.any(np.isnan(uh.x.array)):
    nan_dofs = np.where(np.isnan(uh.x.array))[0]
    print(f"NaN at {len(nan_dofs)} DOFs (first 5: {nan_dofs[:5]})")
    # Map DOF to spatial location for debugging
```
**Fix**:
1. Check all boundary conditions are applied.
2. Verify mesh quality (no inverted or zero-volume elements).
3. Check source terms and material properties for division by zero.
4. If using Michaelis-Menten kinetics, verify regularization is active:
   ```python
   # Bad: Vmax * c / (Km + c)  -- can blow up if Km + c ~ 0
   # Good: Vmax * c / (Km + c + eps)  where eps = 1e-15
   ```

### Negative concentrations in species transport

**Cause**: Standard Galerkin FEM does not guarantee a discrete maximum principle.
Convection-dominated transport (high Peclet) or stiff reactions can produce
undershoots below zero.
**Fix** (in order of preference):
1. **Verify SUPG is active**: Check that the SUPG stabilization term is included
   in the variational form and that the stabilization parameter tau is positive.
2. **Increase mesh resolution**: Especially near inlet boundaries and reaction zones.
3. **Apply clipping** (post-processing, not physically rigorous but practical):
   ```python
   ch.x.array[ch.x.array < 0] = 0.0
   ```
4. **Use DG (Discontinuous Galerkin) elements**: They provide better local conservation
   and are less prone to oscillations, but increase DOF count and implementation complexity.

---

## Stage 5-6: Validation and Visualization Errors

### Conservation violation > 1%

**Cause**: The variational form or boundary conditions have an error, or the mesh
is too coarse to resolve the solution features.
**Diagnostic**: Run the conservation check protocol from validation-benchmarks.md.
**Fix**:
1. Double-check the variational form: verify the sign of each term.
2. Verify boundary condition tags match the correct physical boundaries.
3. Refine the mesh (especially near boundaries with large gradients).
4. Check units: a common error is mixing SI and non-SI units.

### Mesh convergence study shows no convergence

**Cause**: Usually a geometry singularity (sharp corner, re-entrant corner) limiting
the convergence rate, or a bug in the error computation.
**Diagnostic**:
1. Check that the error actually changes between refinement levels.
   If error is flat at machine precision, the solution is in the FE space (see
   Poiseuille note in validation-benchmarks.md).
2. If error decreases but rate is < 1, check for corner singularities.
**Fix**:
1. Use a benchmark with known convergence (e.g., diffusion-reaction from
   validation-benchmarks.md) to rule out code bugs.
2. For singularities: use adaptive mesh refinement or graded meshes near corners.
3. Verify the exact solution is correctly interpolated onto the FE space.

### PyVista: empty plot or all-white rendering

**Cause**: The data array is not attached to the VTK grid, or the data range is
degenerate (min == max).
**Diagnostic**:
```python
import pyvista
grid = pyvista.read("output.vtk")
print(grid.array_names)       # should list "velocity", "pressure", etc.
print(grid["velocity"].min(), grid["velocity"].max())
```
**Fix**:
1. Verify the XDMF/VTK file was written correctly:
   ```python
   with io.XDMFFile(domain.comm, "output.xdmf", "w") as xdmf:
       xdmf.write_mesh(domain)
       xdmf.write_function(uh, t=0.0)  # make sure function is passed
   ```
2. If min == max, the solution is likely zero everywhere (solver did not run
   or BCs are wrong).

### PyVista: slow rendering or freezing

**Cause**: Trying to render a very large mesh (> 1M cells) interactively.
**Fix**:
1. Subsample for visualization:
   ```python
   grid_decimated = grid.decimate(0.9)  # keep 10% of triangles
   ```
2. Use off-screen rendering and save to file:
   ```python
   pyvista.OFF_SCREEN = True
   plotter = pyvista.Plotter()
   plotter.add_mesh(grid, scalars="velocity")
   plotter.screenshot("velocity.png")
   ```
3. Extract a 2D slice instead of rendering the full 3D mesh:
   ```python
   slice_plane = grid.slice(normal="z", origin=(0, 0, 0.5e-3))
   plotter.add_mesh(slice_plane, scalars="velocity")
   ```

### No streamlines generated

**Cause**: The velocity field is zero, very small, or the seed points are outside
the flow domain.
**Diagnostic**:
```python
vel = grid["velocity"]
print(f"Velocity magnitude: min={np.linalg.norm(vel, axis=1).min():.2e}, "
      f"max={np.linalg.norm(vel, axis=1).max():.2e}")
```
**Fix**:
1. If velocity is zero: the flow solver did not converge or BCs are wrong.
2. If velocity is very small (< 1e-15): check units and pressure drop magnitude.
3. Ensure seed points are inside the domain:
   ```python
   streamlines = grid.streamlines(
       vectors="velocity",
       source_center=grid.center,  # safe default
       source_radius=H/4,
       n_points=20,
       max_time=L / u_mean,  # integration time ~ domain transit
   )
   ```

---

## Quick Diagnostic Checklist

When a simulation fails, work through this checklist:

### 1. Which stage failed?

| Symptom | Stage | Go to section |
|---------|-------|---------------|
| Cannot import modules | Stage 0 | Environment errors |
| STEP file won't load | Stage 1 | Geometry errors |
| Meshing fails or crashes | Stage 1 | Mesh errors |
| TypeError / ValueError during setup | Stage 2-3 | Physics setup |
| Solver diverges or gives NaN | Stage 4 | Solver execution |
| Results look wrong | Stage 5-6 | Validation |
| Visualization blank or slow | Stage 5-6 | Visualization |

### 2. Is the environment correct?

```python
import dolfinx; print(dolfinx.__version__)   # expect 0.10.x
import gmsh; print(gmsh.__version__)          # expect 4.x
import pyvista; print(pyvista.__version__)    # expect 0.42+
import petsc4py; print(petsc4py.__version__)  # expect 3.x
```

### 3. Is the mesh valid?

```python
# Check basic mesh properties
print(f"Cells: {domain.topology.index_map(domain.topology.dim).size_global}")
print(f"Vertices: {domain.topology.index_map(0).size_global}")
print(f"Geometric dimension: {domain.geometry.dim}")
```

### 4. Are boundary conditions correctly applied?

```python
# Print DOF counts for each BC
for i, bc in enumerate(bcs):
    print(f"BC {i}: {len(bc.dof_indices()[0])} DOFs constrained")
# If 0 DOFs: the boundary marker function is not selecting any entities
```

### 5. Are units consistent?

All quantities must be in SI base units:

| Quantity | Unit | Common error |
|----------|------|--------------|
| Length | m | Using mm (factor 1e-3) |
| Velocity | m/s | Using cm/s |
| Pressure | Pa | Using bar (factor 1e5) or mmHg |
| Viscosity | Pa.s | Using cP (factor 1e-3) |
| Diffusivity | m^2/s | Using cm^2/s (factor 1e-4) |
| Concentration | mol/m^3 | Using mM (= mol/m^3, no conversion needed) |
| Reaction rate | 1/s or mol/(m^3.s) | Mixing volumetric and surface rates |
