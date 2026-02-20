# Environment Setup: FEniCSx for Bioprocess CFD

This document covers installation, platform-specific configuration, and pre-flight
validation for the FEniCSx-based bioreactor simulation environment. Every generated
script depends on the packages described here. If the environment is not correctly
set up, no simulation will run.

---

## 1. Conda Installation (Primary Method)

Conda is the recommended installation path. The provided `environment.yml` pins
FEniCSx v0.10 and all transitive dependencies to tested, compatible versions.

### 1.1 One-Command Install

```bash
# From the skill examples directory (where environment.yml lives):
conda env create -f environment.yml

# Activate the environment:
conda activate cfd-bioreactor
```

If you already have the environment and need to update it:

```bash
conda env update -f environment.yml --prune
```

### 1.2 Verification Commands

After installation, verify each critical package individually:

```bash
# FEniCSx core
python -c "import dolfinx; print(f'dolfinx {dolfinx.__version__}')"

# Basix (finite element definitions)
python -c "import basix; print(f'basix {basix.__version__}')"

# UFL (unified form language)
python -c "import ufl; print(f'ufl {ufl.__version__}')"

# gmsh (mesh generation)
python -c "import gmsh; gmsh.initialize(); print(f'gmsh {gmsh.GMSH_API_VERSION}'); gmsh.finalize()"

# PETSc (linear algebra backend)
python -c "from petsc4py import PETSc; print(f'PETSc {PETSc.Sys.getVersion()}')"

# PyVista (visualization)
python -c "import pyvista; print(f'pyvista {pyvista.__version__}')"

# MPI
python -c "from mpi4py import MPI; print(f'mpi4py comm_size={MPI.COMM_WORLD.Get_size()}')"
```

All commands must succeed without error. If any fail, see Section 6 (Graceful
Degradation Modes) for partial-installation guidance.

### 1.3 Why conda-forge Exclusively

FEniCSx, PETSc, and gmsh are compiled C/C++ libraries with Python bindings. Mixing
channels (e.g., defaults + conda-forge) causes ABI incompatibilities that manifest
as segfaults during mesh import or solver execution. Always use `conda-forge` as
the sole channel.

---

## 2. Docker Alternative

Docker provides a self-contained environment with zero dependency management. Prefer
Docker when conda installation fails (especially on Windows or HPC clusters with
restricted permissions).

### 2.1 Quick Start

```bash
# Pull the official FEniCSx v0.10 image
docker pull dolfinx/dolfinx:v0.10.0.0

# Run interactively with the current directory mounted
docker run --init -ti \
  -v "$(pwd)":/work \
  -w /work \
  -p 8888:8888 \
  dolfinx/dolfinx:v0.10.0.0
```

Inside the container, all FEniCSx packages are pre-installed. You may need to
install PyVista and psutil additionally:

```bash
pip install pyvista psutil meshio
```

### 2.2 Docker Compose for Persistent Workspace

For repeated use, create a `docker-compose.yml`:

```yaml
version: "3.8"
services:
  cfd:
    image: dolfinx/dolfinx:v0.10.0.0
    volumes:
      - ./simulations:/work
    working_dir: /work
    ports:
      - "8888:8888"
    init: true
    stdin_open: true
    tty: true
```

Then run:

```bash
docker compose up -d
docker compose exec cfd bash
```

### 2.3 When to Use Docker Over Conda

- Windows (WSL2 Docker Desktop is simpler than WSL2 + conda)
- HPC clusters where conda-forge packages conflict with system MPI
- CI/CD pipelines (reproducible container builds)
- When conda environment creation fails due to solver conflicts

---

## 3. Platform-Specific Notes

### 3.1 macOS (Apple Silicon -- M1/M2/M3/M4)

Conda-forge provides ARM64 (osx-arm64) builds for FEniCSx and most dependencies.
However, availability can change between releases.

**Before attempting conda install**, verify ARM64 availability:

```bash
conda search -c conda-forge fenics-dolfinx --platform osx-arm64
```

If ARM64 builds are not yet available for v0.10, use Docker instead (see Section 2).

Known considerations:
- The `mpich` ARM64 build works correctly on Apple Silicon.
- OpenCASCADE (used by gmsh for STEP import) is available on ARM64 via conda-forge
  but should be tested with the pre-flight script (Section 4) after installation.
- PyVista rendering uses the native Metal backend; set `PYVISTA_OFF_SCREEN=true`
  if running headless.

### 3.2 macOS (Intel)

Fully supported via conda-forge. No special considerations. The environment.yml
works without modification.

### 3.3 Linux (Ubuntu 22.04+)

Two options:

**Option A: conda-forge (recommended for version control)**

Use the environment.yml as described in Section 1. This provides exact version
pinning independent of the system package manager.

**Option B: System packages (Ubuntu 24.04+)**

Ubuntu 24.04 provides FEniCSx packages in the official repositories:

```bash
sudo apt install python3-dolfinx python3-gmsh python3-pyvista
```

Note: System packages may lag behind the latest release. Verify the installed
version matches v0.10 before proceeding. Mixing system FEniCSx with conda
packages will cause conflicts.

### 3.4 Windows (WSL2 Required)

Native Windows is not supported. FEniCSx requires a POSIX environment.

**Setup steps**:

1. Install WSL2 with Ubuntu 22.04+:
   ```powershell
   wsl --install -d Ubuntu-24.04
   ```

2. Inside WSL2, install Miniconda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

3. Proceed with conda installation (Section 1).

Alternatively, use Docker Desktop with WSL2 backend (Section 2).

**WSL2 memory note**: By default, WSL2 is limited to 50% of system RAM. For 3D
simulations, increase this in `%USERPROFILE%\.wslconfig`:

```ini
[wsl2]
memory=12GB
```

---

## 4. Pre-Flight Validation Script

Run this script after installation to verify that all components are functional.
Save it as `preflight_check.py` and execute with `python preflight_check.py`.

```python
#!/usr/bin/env python3
"""Pre-flight validation for the CFD bioreactor simulation environment.

Checks all required and optional dependencies, reports pass/fail for each,
and provides actionable fix instructions on failure.

Exit code 0: all critical checks pass.
Exit code 1: one or more critical checks failed.
"""

import sys
import importlib

PASS = "PASS"
FAIL = "FAIL"
WARN = "WARN"

results = []


def check(name, status, message=""):
    """Record a check result."""
    results.append((name, status, message))
    symbol = {"PASS": "+", "FAIL": "!", "WARN": "~"}[status]
    print(f"  [{symbol}] {name}: {status}" + (f" -- {message}" if message else ""))


def main():
    print("=" * 64)
    print("CFD Bioreactor Environment Pre-Flight Check")
    print("=" * 64)
    print()

    # ---- Check 1: Python version ----
    print("Checking Python version...")
    v = sys.version_info
    if v >= (3, 10):
        check("Python version", PASS, f"{v.major}.{v.minor}.{v.micro}")
    else:
        check("Python version", FAIL,
              f"{v.major}.{v.minor}.{v.micro} (need >= 3.10). "
              "Install Python 3.10+ or use the conda environment.")
    print()

    # ---- Check 2: dolfinx ----
    print("Checking FEniCSx (dolfinx)...")
    try:
        import dolfinx
        ver = dolfinx.__version__
        major_minor = tuple(int(x) for x in ver.split(".")[:2])
        if major_minor >= (0, 10):
            check("dolfinx", PASS, f"version {ver}")
        else:
            check("dolfinx", FAIL,
                  f"version {ver} (need >= 0.10). "
                  "Run: conda install -c conda-forge fenics-dolfinx=0.10.*")
    except ImportError:
        check("dolfinx", FAIL,
              "Not installed. Run: conda install -c conda-forge fenics-dolfinx=0.10.* "
              "OR use Docker: docker pull dolfinx/dolfinx:v0.10.0.0")
    print()

    # ---- Check 3: gmsh + OpenCASCADE kernel ----
    print("Checking gmsh and OpenCASCADE kernel...")
    try:
        import gmsh
        gmsh.initialize()
        api_ver = gmsh.GMSH_API_VERSION
        check("gmsh import", PASS, f"API version {api_ver}")

        # Test OCC kernel: create a box, synchronize, and verify entity exists
        try:
            gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
            gmsh.model.occ.synchronize()
            entities = gmsh.model.getEntities(dim=3)
            if len(entities) >= 1:
                check("gmsh OCC kernel", PASS,
                      "addBox + synchronize + entity verification succeeded")
            else:
                check("gmsh OCC kernel", FAIL,
                      "addBox succeeded but no 3D entities found after synchronize. "
                      "OpenCASCADE kernel may be broken. "
                      "Reinstall: conda install -c conda-forge gmsh python-gmsh")
        except Exception as e:
            check("gmsh OCC kernel", FAIL,
                  f"OpenCASCADE test failed: {e}. "
                  "STEP/IGES import will not work. "
                  "Reinstall: conda install -c conda-forge gmsh python-gmsh")
        finally:
            gmsh.finalize()
    except ImportError:
        check("gmsh import", FAIL,
              "Not installed. Run: conda install -c conda-forge gmsh python-gmsh")
        check("gmsh OCC kernel", FAIL, "Skipped (gmsh not available)")
    print()

    # ---- Check 4: PyVista ----
    print("Checking PyVista...")
    try:
        import pyvista as pv
        check("pyvista", PASS, f"version {pv.__version__}")

        # Detect rendering backend
        try:
            plotter = pv.Plotter(off_screen=True)
            plotter.close()
            check("pyvista rendering", PASS, "off-screen rendering available")
        except Exception as e:
            check("pyvista rendering", WARN,
                  f"Off-screen rendering failed: {e}. "
                  "Visualization will require export to VTK and external viewer. "
                  "Try: pip install pyvista[jupyter] OR set DISPLAY env variable.")
    except ImportError:
        check("pyvista", WARN,
              "Not installed. Visualization will use VTK export only. "
              "Install: conda install -c conda-forge pyvista")
    print()

    # ---- Check 5: PETSc + MUMPS solver ----
    print("Checking PETSc and MUMPS solver...")
    try:
        from petsc4py import PETSc
        petsc_ver = ".".join(str(x) for x in PETSc.Sys.getVersion())
        check("petsc4py", PASS, f"version {petsc_ver}")

        # Check MUMPS availability
        mumps_available = False
        try:
            # Method 1: Use hasExternalPackage if available
            if hasattr(PETSc.Sys, "hasExternalPackage"):
                mumps_available = PETSc.Sys.hasExternalPackage("mumps")
            else:
                # Method 2: Try creating a MUMPS solver instance
                ksp = PETSc.KSP().create()
                ksp.setType("preonly")
                pc = ksp.getPC()
                pc.setType("lu")
                pc.setFactorSolverType("mumps")
                mumps_available = True
                ksp.destroy()
        except Exception:
            mumps_available = False

        if mumps_available:
            check("MUMPS solver", PASS, "Direct solver available")
        else:
            check("MUMPS solver", WARN,
                  "MUMPS not detected. Direct solver unavailable. "
                  "2D problems will use LU (slower). 3D problems will use "
                  "iterative solvers. Reinstall PETSc with MUMPS: "
                  "conda install -c conda-forge petsc4py=*=*mumps*")
    except ImportError:
        check("petsc4py", FAIL,
              "Not installed. Run: conda install -c conda-forge petsc4py")
        check("MUMPS solver", FAIL, "Skipped (petsc4py not available)")
    print()

    # ---- Check 6: mpi4py ----
    print("Checking MPI...")
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        check("mpi4py", PASS,
              f"rank {comm.Get_rank()} of {comm.Get_size()}")
    except ImportError:
        check("mpi4py", FAIL,
              "Not installed. Run: conda install -c conda-forge mpi4py mpich")
    print()

    # ---- Check 7: System RAM ----
    print("Checking system resources...")
    try:
        import psutil
        ram_gb = psutil.virtual_memory().total / (1024 ** 3)
        if ram_gb >= 16:
            check("System RAM", PASS,
                  f"{ram_gb:.1f} GB (sufficient for 3D simulations)")
        elif ram_gb >= 8:
            check("System RAM", WARN,
                  f"{ram_gb:.1f} GB (sufficient for 2D; 3D may require "
                  "coarse meshes or iterative solvers)")
        else:
            check("System RAM", WARN,
                  f"{ram_gb:.1f} GB (limited to small 2D problems; "
                  "consider a machine with >= 16 GB for 3D)")
    except ImportError:
        check("System RAM", WARN,
              "psutil not installed; cannot check RAM. "
              "Install: pip install psutil")
    print()

    # ---- Summary ----
    print("=" * 64)
    print("SUMMARY")
    print("=" * 64)

    critical_failures = [r for r in results if r[1] == FAIL]
    warnings = [r for r in results if r[1] == WARN]
    passes = [r for r in results if r[1] == PASS]

    print(f"  Passed:   {len(passes)}")
    print(f"  Warnings: {len(warnings)}")
    print(f"  Failed:   {len(critical_failures)}")
    print()

    if critical_failures:
        print("CRITICAL FAILURES (must fix before running simulations):")
        for name, _, msg in critical_failures:
            print(f"  - {name}: {msg}")
        print()
        print("STATUS: FAIL -- resolve the issues above and re-run this script.")
        return 1

    if warnings:
        print("WARNINGS (simulation will run with reduced capabilities):")
        for name, _, msg in warnings:
            print(f"  - {name}: {msg}")
        print()

    print("STATUS: PASS -- environment is ready for CFD bioreactor simulation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

---

## 5. Version Compatibility Matrix

FEniCSx is a multi-package ecosystem. All components must be from the same release
series to ensure ABI and API compatibility.

| Component       | Required Version | Notes                                    |
|-----------------|------------------|------------------------------------------|
| dolfinx         | 0.10.*           | Core finite element library              |
| basix           | 0.10.*           | Finite element definitions               |
| UFL             | 2024.2.*         | Form language (versioned separately)     |
| FFCx            | 0.10.*           | Form compiler                            |
| PETSc           | 3.21+            | Linear algebra backend                   |
| petsc4py        | 3.21+            | Python bindings for PETSc                |
| gmsh            | >= 4.11          | Mesh generation (OpenCASCADE required)   |
| python-gmsh     | >= 4.11          | Python API for gmsh                      |
| PyVista         | >= 0.42          | 3D visualization                         |
| numpy           | >= 1.24          | Array operations                         |
| scipy           | >= 1.10          | Sparse matrix utilities                  |
| mpi4py          | >= 3.1           | MPI parallel communication               |
| Python          | >= 3.10          | Required by dolfinx 0.10                 |

**Critical rule**: All dolfinx/basix/FFCx packages must share the same minor
version (0.10). Mixing 0.9 basix with 0.10 dolfinx will produce import errors
or silent numerical errors.

The conda `environment.yml` enforces these constraints automatically. If installing
manually, verify with:

```bash
python -c "
import dolfinx, basix
print(f'dolfinx: {dolfinx.__version__}')
print(f'basix:   {basix.__version__}')
assert dolfinx.__version__.startswith('0.10'), 'Version mismatch!'
assert basix.__version__.startswith('0.10'), 'Version mismatch!'
print('Version compatibility: OK')
"
```

---

## 6. Graceful Degradation Modes

Not every component is required for every task. The skill adapts its capabilities
based on what is available.

### 6.1 dolfinx Missing -- Mesh-Only Mode

If FEniCSx is not installed but gmsh is available, the skill can still:
- Import STEP/IGES geometry
- Generate and refine meshes
- Assign physical groups (boundary markers)
- Export meshes to `.msh`, `.vtk`, or `.xdmf` format
- Perform mesh quality analysis

The skill cannot solve any physics equations without dolfinx.

**What to tell the user**: "FEniCSx is not installed. I can generate meshes and
export them for use in another solver. To run simulations, install FEniCSx
following the environment setup guide."

### 6.2 PyVista Missing -- Export-to-VTK Mode

If PyVista is not installed, the skill:
- Runs all simulations normally
- Exports results to `.pvd` / `.xdmf` files instead of interactive plots
- Provides ParaView instructions for visualization

**What to tell the user**: "PyVista is not available. Results have been exported
to [filename].pvd. Open this file in ParaView for visualization."

### 6.3 gmsh OCC Kernel Missing -- Parametric Geometry Only

If gmsh is installed but the OpenCASCADE kernel is non-functional:
- STEP/IGES import will fail
- The skill can still create parametric geometry using the built-in gmsh kernel
  (rectangles, circles, extrusions) -- sufficient for 2D channel flow and simple
  3D geometries
- Mesh generation from parametric geometry works normally

**What to tell the user**: "The gmsh OpenCASCADE kernel is not available. STEP
file import is disabled. I can create parametric geometries (rectangles, cylinders,
channels) directly. To enable STEP import, reinstall gmsh:
`conda install -c conda-forge gmsh python-gmsh`."

### 6.4 Nothing Installed -- Setup Instructions Mode

If no scientific computing packages are available:
- Display the full environment setup instructions
- Offer to generate the `environment.yml` file
- Offer the Docker alternative
- **STOP** -- do not attempt to generate simulation scripts that cannot run

**What to tell the user**: "No simulation packages are installed. Before we can
run any CFD simulations, you need to set up the environment. Here are your
options: [display conda and Docker instructions from this document]."

### Degradation Summary Table

| Missing Component | Capability Lost          | Remaining Capabilities          |
|-------------------|--------------------------|---------------------------------|
| dolfinx           | All physics solving      | Mesh generation, geometry, export |
| PyVista           | Interactive visualization | Simulation + VTK/XDMF export   |
| gmsh OCC          | STEP/IGES import         | Parametric geometry + simulation |
| gmsh (entirely)   | All mesh generation      | Nothing useful -- install gmsh  |
| petsc4py          | All solving              | Mesh generation only            |
| MUMPS             | Direct solver            | Iterative solvers (CG, GMRES)  |
| Everything        | Everything               | Display setup instructions, STOP |
