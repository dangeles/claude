# Mesh Generation Guide for Bioprocess Cartridges

Reference file for the `cfd-bioreactor` skill. Covers gmsh Python API patterns
for STEP file import, physical group definition, mesh refinement, quality
assessment, memory estimation, and conversion to DOLFINx.

All code targets gmsh >= 4.11 with the Python API and DOLFINx v0.10.

---

## 1. gmsh Python API Fundamentals

### Geometry Kernel Selection

gmsh provides two geometry kernels:

| Kernel | Use Case | API Prefix |
|--------|----------|------------|
| OpenCASCADE (OCC) | STEP/IGES import, boolean operations, CAD interop | `gmsh.model.occ` |
| Built-in | Parametric geometry from scratch, simple shapes | `gmsh.model.geo` |

**Rule**: If importing a STEP file, always use OCC. If constructing geometry
programmatically from dimensions, prefer built-in for simpler shapes or OCC
for boolean operations.

**Never mix kernels** in a single model. Pick one and use it throughout.

### Initialization and Finalization Pattern

Every gmsh script must follow this pattern:

```python
import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)  # Enable terminal output
gmsh.option.setNumber("General.Verbosity", 3)  # 0=silent, 3=normal, 99=debug
gmsh.model.add("bioreactor")

# ... all geometry, meshing, and export operations ...

gmsh.finalize()  # MUST call to free resources
```

Use a try/finally block to guarantee finalization:

```python
import gmsh

gmsh.initialize()
try:
    gmsh.model.add("bioreactor")
    # ... operations ...
finally:
    gmsh.finalize()
```

---

## 2. STEP File Import with Error Handling

### Complete Import Pattern

```python
import gmsh
import os

def import_step_file(step_path: str) -> dict:
    """Import a STEP file and return diagnostic information.

    Returns a dict with keys: volumes, surfaces, curves, points,
    bounding_box, dimensionality.
    """
    if not os.path.isfile(step_path):
        raise FileNotFoundError(f"STEP file not found: {step_path}")

    file_size_mb = os.path.getsize(step_path) / 1e6
    if file_size_mb > 200:
        print(f"WARNING: STEP file is {file_size_mb:.0f} MB. "
              "Import may be very slow and memory-intensive. "
              "Consider simplifying the CAD model first.")
    elif file_size_mb > 50:
        print(f"NOTE: STEP file is {file_size_mb:.0f} MB. "
              "Import may take a few minutes.")

    try:
        entities = gmsh.model.occ.importShapes(step_path)
    except Exception as e:
        raise RuntimeError(
            f"Failed to import STEP file: {e}\n"
            "Troubleshooting:\n"
            "  1. Re-export from CAD as STEP AP214 (most compatible)\n"
            "  2. Remove small features (fillets < 0.1 mm, chamfers)\n"
            "  3. Export only the fluid domain (not solid parts)\n"
            "  4. Check file is not corrupted (open in FreeCAD first)\n"
            "  5. Try IGES format as fallback"
        )

    gmsh.model.occ.synchronize()

    # Post-import diagnostics
    volumes = gmsh.model.getEntities(dim=3)
    surfaces = gmsh.model.getEntities(dim=2)
    curves = gmsh.model.getEntities(dim=1)
    points = gmsh.model.getEntities(dim=0)

    print(f"Imported: {len(volumes)} volumes, {len(surfaces)} surfaces, "
          f"{len(curves)} curves, {len(points)} points")

    # Bounding box
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin
    print(f"Bounding box: [{xmin:.4e}, {ymin:.4e}, {zmin:.4e}] "
          f"to [{xmax:.4e}, {ymax:.4e}, {zmax:.4e}]")
    print(f"Dimensions: Lx={Lx:.4e}, Ly={Ly:.4e}, Lz={Lz:.4e} m")

    # Unit mismatch warning
    max_dim = max(Lx, Ly, Lz)
    if max_dim > 1.0:
        print(f"WARNING: Maximum dimension is {max_dim:.2f} m ({max_dim*1e3:.1f} mm). "
              "If the geometry was exported in millimeters, you must scale "
              "to meters. Bioprocess cartridges are typically 1-100 mm. "
              "Use gmsh.model.occ.dilate() to rescale if needed:")
        print("  gmsh.model.occ.dilate(entities, 0, 0, 0, 1e-3, 1e-3, 1e-3)")
    elif max_dim < 1e-4:
        print(f"WARNING: Maximum dimension is {max_dim:.2e} m ({max_dim*1e6:.1f} um). "
              "This is unusually small. Verify that the STEP file uses "
              "meter units. If in micrometers, scale up accordingly.")

    # Dimensionality check
    dimensionality = 3 if len(volumes) > 0 else (2 if len(surfaces) > 0 else 1)
    if dimensionality < 3:
        print(f"NOTE: No 3D volumes found. Geometry is {dimensionality}D. "
              "Options: (a) mesh as 2D, (b) extrude to 3D, "
              "(c) re-export STEP with solid bodies.")

    # Multi-body detection
    if len(volumes) > 1:
        print(f"NOTE: Found {len(volumes)} separate volumes. "
              "This may indicate a multi-body assembly. "
              "Options: (a) boolean union to merge, "
              "(b) assign separate physical groups, "
              "(c) select only the fluid domain volume.")

    return {
        "volumes": volumes,
        "surfaces": surfaces,
        "curves": curves,
        "points": points,
        "bounding_box": (xmin, ymin, zmin, xmax, ymax, zmax),
        "dimensionality": dimensionality,
    }
```

### Fallback: Parametric Geometry

When STEP import fails or is unavailable, construct geometry programmatically:

```python
# Simple 2D channel (built-in kernel)
gmsh.model.geo.addPoint(0, 0, 0, meshSize=0.1, tag=1)
gmsh.model.geo.addPoint(L, 0, 0, meshSize=0.1, tag=2)
gmsh.model.geo.addPoint(L, H, 0, meshSize=0.1, tag=3)
gmsh.model.geo.addPoint(0, H, 0, meshSize=0.1, tag=4)
gmsh.model.geo.addLine(1, 2, tag=1)  # bottom wall
gmsh.model.geo.addLine(2, 3, tag=2)  # outlet
gmsh.model.geo.addLine(3, 4, tag=3)  # top wall (or membrane)
gmsh.model.geo.addLine(4, 1, tag=4)  # inlet
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], tag=1)
gmsh.model.geo.addPlaneSurface([1], tag=1)
gmsh.model.geo.synchronize()
```

---

## 3. Physical Group Definition

### Physical Groups Are MANDATORY

Physical groups map geometric entities to named boundary/subdomain markers.
Without them, DOLFINx has no way to identify where to apply boundary
conditions or subdomain-specific source terms. Defining physical groups is
MANDATORY before meshing.

**Anti-pattern** (NEVER do this):
```python
# BAD: meshing without physical groups
gmsh.model.mesh.generate(3)
# Result: DOLFINx mesh has no tags, BCs cannot be applied
```

**Correct pattern**: always define physical groups BEFORE meshing.

### Naming Convention

Use consistent integer markers across all generated scripts
(inlet=1, outlet=2, walls=3, membrane=4, cell_region=5, fluid_region=6):

| Physical Group | Marker ID | Dimension | Description |
|---------------|-----------|-----------|-------------|
| inlet | 1 | facet (dim-1) | Inlet boundary |
| outlet | 2 | facet (dim-1) | Outlet boundary |
| walls | 3 | facet (dim-1) | No-slip / no-flux walls |
| membrane | 4 | facet (dim-1) | Permeable membrane surface |
| cell_region | 5 | cell (dim) | Subdomain containing cells |
| fluid_region | 6 | cell (dim) | Bulk fluid subdomain |

Additional markers (7+) can be assigned for problem-specific features.

### Complete Code Pattern

```python
def define_physical_groups_2d():
    """Define physical groups for a 2D channel geometry.

    Assumes geometry has 4 boundary curves:
    - curve 4 = inlet (left)
    - curve 2 = outlet (right)
    - curve 1 = bottom wall
    - curve 3 = top wall (or membrane)
    And 1 surface = fluid domain.
    """
    # Facet groups (1D curves in 2D)
    gmsh.model.addPhysicalGroup(1, [4], tag=1)
    gmsh.model.setPhysicalName(1, 1, "inlet")

    gmsh.model.addPhysicalGroup(1, [2], tag=2)
    gmsh.model.setPhysicalName(1, 2, "outlet")

    gmsh.model.addPhysicalGroup(1, [1], tag=3)
    gmsh.model.setPhysicalName(1, 3, "walls")

    gmsh.model.addPhysicalGroup(1, [3], tag=4)
    gmsh.model.setPhysicalName(1, 4, "membrane")

    # Cell group (2D surface)
    gmsh.model.addPhysicalGroup(2, [1], tag=6)
    gmsh.model.setPhysicalName(2, 6, "fluid_region")
```

### Identifying Surfaces by Location

For imported STEP geometry where entity tags are not known a priori, use
bounding box queries to identify surfaces:

```python
def find_surfaces_at_position(coord: str, value: float, tol: float = 1e-6):
    """Find surfaces whose bounding box center is at a given coordinate value.

    Args:
        coord: "x", "y", or "z"
        value: coordinate value to match
        tol: tolerance for matching

    Returns:
        List of surface tags matching the criterion.
    """
    idx = {"x": 0, "y": 1, "z": 2}[coord]
    matching = []
    for dim, tag in gmsh.model.getEntities(dim=2):
        bbox = gmsh.model.getBoundingBox(dim, tag)
        # bbox = (xmin, ymin, zmin, xmax, ymax, zmax)
        center = (bbox[idx] + bbox[idx + 3]) / 2.0
        if abs(center - value) < tol:
            matching.append(tag)
    return matching

# Example: find inlet surfaces at x = 0
inlet_surfaces = find_surfaces_at_position("x", 0.0, tol=1e-4)
gmsh.model.addPhysicalGroup(2, inlet_surfaces, tag=1)
gmsh.model.setPhysicalName(2, 1, "inlet")
```

### Volume Groups for Subdomain Markers

```python
# For multi-region problems: assign cell_region to specific volumes
cell_volume_tags = [2, 3]  # volumes containing cell scaffold
fluid_volume_tags = [1]     # bulk fluid volume

gmsh.model.addPhysicalGroup(3, cell_volume_tags, tag=5)
gmsh.model.setPhysicalName(3, 5, "cell_region")

gmsh.model.addPhysicalGroup(3, fluid_volume_tags, tag=6)
gmsh.model.setPhysicalName(3, 6, "fluid_region")
```

---

## 4. Mesh Refinement Strategies

### 4.1 Distance + Threshold Field (Beginner-Friendly)

This approach refines the mesh near specified geometric entities using a
smooth transition from fine to coarse elements.

```python
def add_refinement_near_surfaces(surface_tags: list,
                                  h_min: float,
                                  h_max: float,
                                  dist_min: float,
                                  dist_max: float):
    """Add mesh refinement near specified surfaces.

    Args:
        surface_tags: gmsh surface tags to refine near
        h_min: minimum element size (at the surface)
        h_max: maximum element size (far from the surface)
        dist_min: distance below which h = h_min
        dist_max: distance above which h = h_max
    """
    # Distance field: computes distance from surfaces
    dist_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist_field, "SurfacesList", surface_tags)
    gmsh.model.mesh.field.setNumber(dist_field, "Sampling", 100)

    # Threshold field: maps distance to element size
    thresh_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thresh_field, "InField", dist_field)
    gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin", h_min)
    gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax", h_max)
    gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", dist_min)
    gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", dist_max)

    return thresh_field
```

**Typical usage for bioprocess cartridges**:

```python
L_char = 1e-3  # characteristic length (1 mm channel height)

# Refine near membrane interface
f1 = add_refinement_near_surfaces(
    membrane_surfaces,
    h_min=L_char / 50,     # fine: 20 um
    h_max=L_char / 5,      # coarse: 200 um
    dist_min=L_char / 20,  # 50 um
    dist_max=L_char / 2,   # 500 um
)

# Refine near inlet/outlet
f2 = add_refinement_near_surfaces(
    inlet_surfaces + outlet_surfaces,
    h_min=L_char / 20,
    h_max=L_char / 5,
    dist_min=L_char / 10,
    dist_max=L_char,
)

# Combine fields: take minimum size
min_field = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [f1, f2])
gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
```

### 4.2 BoundaryLayer Field (Production Quality)

For resolving velocity and concentration boundary layers near walls and
membranes. Produces structured prismatic layers in 3D (or quadrilateral
layers in 2D) adjacent to surfaces.

```python
def add_boundary_layer(surface_tags: list,
                        hwall_n: float,
                        ratio: float = 1.3,
                        thickness: float = None,
                        nb_layers: int = 8):
    """Add boundary layer mesh near surfaces.

    Args:
        surface_tags: surfaces to grow layers from
        hwall_n: first layer height (normal to surface)
        ratio: growth ratio between successive layers
        thickness: total boundary layer thickness (auto-computed if None)
        nb_layers: number of layers
    """
    bl_field = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(bl_field, "SurfacesList", surface_tags)
    gmsh.model.mesh.field.setNumber(bl_field, "Size", hwall_n)
    gmsh.model.mesh.field.setNumber(bl_field, "Ratio", ratio)
    gmsh.model.mesh.field.setNumber(bl_field, "NbLayers", nb_layers)
    if thickness is not None:
        gmsh.model.mesh.field.setNumber(bl_field, "Thickness", thickness)
    gmsh.model.mesh.field.setAsBoundaryLayer(bl_field)
    return bl_field
```

**Typical values for bioprocess cartridges**:

| Parameter | 2D Channel | 3D Cartridge | Units |
|-----------|-----------|-------------|-------|
| hwall_n | L/100 | L/50 | m |
| ratio | 1.2 | 1.3 | - |
| nb_layers | 10 | 8 | - |
| thickness | L/5 | L/4 | m |

Where L is the channel height or characteristic dimension.

### 4.3 Mesh Size Bounds

Always set global mesh size bounds to prevent pathological meshes:

```python
L_char = 1e-3  # characteristic length in meters

# Prevent astronomically fine meshes from small CAD features
gmsh.option.setNumber("Mesh.MeshSizeMin", L_char / 200)

# Prevent excessively coarse meshes
gmsh.option.setNumber("Mesh.MeshSizeMax", L_char / 2)

# Extend mesh size from boundary (helps smooth transition)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
```

**Rationale for MeshSizeMin**: Small CAD features (fillets, chamfers, holes)
can cause gmsh to generate extremely fine elements locally, leading to
millions of unexpected elements. Setting a floor prevents this.

---

## 5. Memory Estimation

### Formula

For P2/P1 (Taylor-Hood) elements on tetrahedra:

```
DOFs per P2 tetrahedron ~ 10 (velocity) + 4 (pressure) = 14
DOFs per P1 tetrahedron (transport) ~ 4
Total DOFs per element ~ 18 (coupled flow + transport)

Matrix memory ~ num_elements * DOFs_per_element^2 * sparsity_factor * 8 bytes
Sparsity factor ~ 50 (typical for FEM with ~50 nonzeros per row)

RAM_GB ~ num_elements * 18 * 50 * 8 / 1e9
       ~ num_elements * 7.2e-6 GB
       ~ num_elements / 140000 GB
```

**Rule of thumb**: approximately 1 GB per 100K P2 tetrahedra for coupled
flow + transport with a direct solver.

### Memory Table by Problem Size

| Dimension | Elements | Approx DOFs | Est. RAM (direct solver) | Est. RAM (iterative) |
|-----------|----------|-------------|------------------------|---------------------|
| 2D | 1,000 | 10K | < 0.1 GB | < 0.1 GB |
| 2D | 10,000 | 100K | 0.2 GB | 0.1 GB |
| 2D | 100,000 | 1M | 3 GB | 1 GB |
| 3D | 10,000 | 100K | 0.5 GB | 0.2 GB |
| 3D | 50,000 | 500K | 4 GB | 1.5 GB |
| 3D | 100,000 | 1M | 10 GB | 3 GB |
| 3D | 500,000 | 5M | 60 GB | 12 GB |
| 3D | 1,000,000 | 10M | > 100 GB | 25 GB |

### System RAM Check

```python
import psutil

def check_memory_for_mesh(num_elements: int, gdim: int = 3) -> bool:
    """Estimate memory requirement and check against available RAM.

    Args:
        num_elements: number of mesh elements
        gdim: geometric dimension (2 or 3)

    Returns:
        True if estimated memory fits in available RAM.
    """
    if gdim == 2:
        gb_per_element = 5e-6  # triangles have fewer DOFs
    else:
        gb_per_element = 7.2e-6  # tetrahedra

    estimated_gb = num_elements * gb_per_element
    available_gb = psutil.virtual_memory().available / 1e9
    total_gb = psutil.virtual_memory().total / 1e9

    print(f"Estimated memory: {estimated_gb:.1f} GB")
    print(f"Available RAM:    {available_gb:.1f} GB / {total_gb:.1f} GB total")

    if estimated_gb > 0.7 * available_gb:
        print(f"WARNING: Estimated memory ({estimated_gb:.1f} GB) exceeds 70% "
              f"of available RAM ({available_gb:.1f} GB). "
              "Consider coarsening the mesh or using an iterative solver.")
        return False

    if estimated_gb > available_gb:
        print(f"ERROR: Estimated memory ({estimated_gb:.1f} GB) exceeds "
              f"available RAM ({available_gb:.1f} GB). "
              "The solver will likely crash. Reduce mesh element count.")
        return False

    print("Memory check: OK")
    return True
```

### Maximum Element Recommendations

| Configuration | Max Elements | Rationale |
|--------------|-------------|-----------|
| Single core, laptop (8 GB) | 50K-100K | Direct solver fits in memory |
| Single core, workstation (32 GB) | 200K-300K | Direct solver comfortable |
| MPI (4 cores, 64 GB) | 500K-1M | Iterative solver with domain decomposition |
| HPC cluster | 1M+ | GAMG preconditioner recommended |

**For this skill's scope**: target 200K elements max for single-core. If the
mesh exceeds this, suggest coarsening or switching to iterative solver.

---

## 6. Mesh Quality Metrics

### Quality Check Code

```python
def check_mesh_quality(element_type: str = "tetrahedra") -> dict:
    """Compute and report mesh quality metrics.

    Returns dict with min, max, mean quality and element count.
    """
    # Get all element types present
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

    total_elements = sum(len(tags) for tags in elem_tags)
    print(f"Total elements: {total_elements}")

    if total_elements == 0:
        raise RuntimeError("No elements found. Did you call gmsh.model.mesh.generate()?")

    # Compute quality (gamma metric: 0 = degenerate, 1 = ideal)
    qualities = gmsh.model.mesh.getElementQualities(qualityType="minSJ")

    if len(qualities) == 0:
        raise RuntimeError("No quality data returned.")

    q_min = min(qualities)
    q_max = max(qualities)
    q_mean = sum(qualities) / len(qualities)
    q_below_01 = sum(1 for q in qualities if q < 0.1) / len(qualities) * 100

    print(f"Quality (min scaled Jacobian):")
    print(f"  Min:  {q_min:.4f}")
    print(f"  Max:  {q_max:.4f}")
    print(f"  Mean: {q_mean:.4f}")
    print(f"  Elements with quality < 0.1: {q_below_01:.1f}%")

    return {
        "total_elements": total_elements,
        "quality_min": q_min,
        "quality_max": q_max,
        "quality_mean": q_mean,
        "pct_below_01": q_below_01,
    }
```

### Quality Thresholds

| Metric | Poor | Acceptable | Good | Excellent |
|--------|------|-----------|------|-----------|
| Min scaled Jacobian | < 0.01 | 0.01 - 0.1 | 0.1 - 0.3 | > 0.3 |
| Aspect ratio | > 20 | 10 - 20 | 5 - 10 | < 5 |
| Mean quality | < 0.3 | 0.3 - 0.5 | 0.5 - 0.7 | > 0.7 |
| % elements < 0.1 quality | > 5% | 1-5% | < 1% | 0% |

**Action thresholds**:
- Min quality < 0.01: Solver will likely fail. Defeature geometry (remove
  small fillets, chamfers) or increase MeshSizeMin.
- More than 5% of elements below 0.1: Expect slow convergence. Consider
  mesh smoothing (`gmsh.model.mesh.optimize("Laplace2D")` or
  `gmsh.model.mesh.optimize("Netgen")`).
- Total elements < 100 (2D) or < 1000 (3D): Mesh is too coarse for
  meaningful results. Reduce MeshSizeMax.

### Mesh Optimization

```python
# Smooth mesh after generation (improves quality without changing topology)
gmsh.model.mesh.optimize("Laplace2D")   # for 2D
gmsh.model.mesh.optimize("Netgen")       # for 3D (if available)
gmsh.model.mesh.optimize("HighOrder")    # for P2 elements
```

---

## 7. gmsh-to-DOLFINx Conversion

### DOLFINx v0.10 API

The correct conversion function in DOLFINx v0.10:

```python
from dolfinx.io import gmsh as gmsh_io
from mpi4py import MPI

# Convert gmsh model to DOLFINx mesh (while gmsh is still initialized)
mesh_data = gmsh_io.model_to_mesh(
    gmsh.model,
    comm=MPI.COMM_WORLD,
    rank=0,
    gdim=3,  # geometric dimension: 2 or 3
)

# Access mesh components
domain = mesh_data.mesh          # dolfinx.mesh.Mesh
cell_tags = mesh_data.cell_tags  # dolfinx.mesh.MeshTags (subdomain markers)
facet_tags = mesh_data.facet_tags  # dolfinx.mesh.MeshTags (boundary markers)
```

**Important**: Call `model_to_mesh` BEFORE `gmsh.finalize()`. The gmsh model
must still be in memory.

### Using MeshTags for Boundary Conditions

```python
import dolfinx.fem as fem

# Identify facets with a given marker
inlet_marker = 1
inlet_facets = facet_tags.find(inlet_marker)

# Locate DOFs on those facets
fdim = domain.topology.dim - 1
inlet_dofs = fem.locate_dofs_topological(V, fdim, inlet_facets)

# Apply Dirichlet BC
bc_inlet = fem.dirichletbc(inlet_value, inlet_dofs, V)
```

### Using MeshTags for Subdomain Integration

```python
import ufl

# Create integration measures with subdomain data
dx = ufl.Measure("dx", domain=domain, subdomain_data=cell_tags)
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

# Integrate over specific subdomains
cell_region_marker = 5
membrane_marker = 4

F_reaction = -inner(R, v) * dx(cell_region_marker)  # reaction in cell region only
F_robin = -P * (c_ext - c) * v * ds(membrane_marker)  # Robin BC on membrane
```

---

## 8. Common Pitfalls

### Pitfall 1: Forgetting Physical Groups

**Symptom**: DOLFINx mesh loads but has no tags; `facet_tags.find(1)` returns
empty array.

**Cause**: Physical groups were not defined before meshing.

**Fix**: Always define physical groups before calling `gmsh.model.mesh.generate()`.

### Pitfall 2: Using `gmshio` Instead of `gmsh`

**Symptom**: `ImportError: cannot import name 'gmshio' from 'dolfinx.io'`

**Cause**: In DOLFINx v0.10, the module was renamed from `dolfinx.io.gmshio`
to `dolfinx.io.gmsh`.

**Fix**: Use `from dolfinx.io import gmsh as gmsh_io` (aliased to avoid
conflict with the `gmsh` package itself).

### Pitfall 3: Not Synchronizing After OCC Operations

**Symptom**: Meshing produces unexpected results or crashes. Entities are
not found when defining physical groups.

**Cause**: OCC operations (import, boolean, extrude) modify an internal
representation that must be synchronized to the gmsh model.

**Fix**: Always call `gmsh.model.occ.synchronize()` after any OCC operation
and before meshing or physical group definition.

### Pitfall 4: Mesh Size Field Conflicts

**Symptom**: Mesh is unexpectedly fine or coarse in regions where a size field
was set.

**Cause**: Multiple active size fields override each other. The
`MeshSizeFromPoints` and `MeshSizeFromCurvature` options compete with
background fields.

**Fix**: Disable default sizing when using custom fields:

```python
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
```

### Pitfall 5: Non-Manifold Geometry from STEP Import

**Symptom**: Meshing fails with "non-manifold" or "self-intersecting" errors.

**Cause**: The STEP file contains touching-but-not-merged surfaces, duplicate
entities, or zero-thickness features.

**Fix**:
1. In CAD software: use "heal" or "repair" tools before exporting
2. In gmsh: try `gmsh.model.occ.removeAllDuplicates()` after import
3. Increase OCC tolerance: `gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)`
4. As last resort: reconstruct geometry parametrically using gmsh built-in kernel

### Pitfall 6: Calling model_to_mesh After gmsh.finalize()

**Symptom**: Segmentation fault or empty mesh.

**Cause**: `gmsh.finalize()` frees the gmsh model from memory. The DOLFINx
converter needs the model to still be active.

**Fix**: Convert to DOLFINx mesh BEFORE calling `gmsh.finalize()`:

```python
gmsh.model.mesh.generate(3)

# Convert FIRST
mesh_data = gmsh_io.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)

# THEN finalize
gmsh.finalize()
```
