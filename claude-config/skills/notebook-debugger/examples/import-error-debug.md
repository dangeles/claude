# Import Error Debug Example: Environment Conflict Resolution

## Problem Report

**User complaint**: "Cell 1 fails with `ModuleNotFoundError: No module named 'scanpy'`. I installed it yesterday!"

**Notebook**: `scrna_analysis.ipynb`
**Cell that fails**: Cell 1 (imports)
**Environment**: Jupyter Lab 3.6, Python 3.11, macOS

---

## Phase 1: Diagnose

### Initial Investigation

**Test 1: Check error message**
```python
# Cell 1:
import scanpy as sc
import numpy as np
import pandas as pd

# Error:
# ModuleNotFoundError: No module named 'scanpy'
```

**Test 2: Check Python environment**
```python
# Cell 0 (add before imports):
import sys
print(f"Python: {sys.version}")
print(f"Executable: {sys.executable}")

# Output:
# Python: 3.11.2 (main, Feb  8 2023, 14:49:24)
# Executable: /Users/name/anaconda3/bin/python
```

**Observation**: Using base anaconda environment, not project-specific environment

**Test 3: Verify package installation**
```bash
# In terminal:
pip list | grep scanpy
# (no output - scanpy not installed in this environment)

# Check where it was installed:
micromamba env list
# * base                     /Users/name/anaconda3
#   scrna_project            /Users/name/anaconda3/envs/scrna_project

micromamba activate scrna_project
pip list | grep scanpy
# scanpy    1.9.3
```

**Root Cause Identified**: Scanpy installed in `scrna_project` environment, but Jupyter kernel is using `base` environment

---

## Phase 2: Isolate

### Hypothesis: Wrong Kernel Selected

**Test: Check available kernels**
```bash
# In terminal:
jupyter kernelspec list

# Output:
# Available kernels:
#   python3    /Users/name/anaconda3/share/jupyter/kernels/python3
```

**Observation**: Only base environment kernel is available. Project environment kernel not registered.

### Verify Environment Has Correct Packages

```bash
# Activate project environment:
micromamba activate scrna_project

# Check all packages:
pip list

# Output includes:
# scanpy                1.9.3
# anndata               0.8.0
# leidenalg             0.9.1
# python-igraph         0.10.4
```

**Conclusion**: Packages are installed correctly, but Jupyter can't see them because kernel not registered

---

## Phase 3: Fix

### Solution: Register Project Environment as Jupyter Kernel

```bash
# Activate project environment:
micromamba activate scrna_project

# Install ipykernel in project environment:
micromamba install ipykernel

# Register environment as Jupyter kernel:
python -m ipykernel install --user --name=scrna_project --display-name="Python (scrna_project)"

# Output:
# Installed kernelspec scrna_project in /Users/name/Library/Jupyter/kernels/scrna_project
```

**Verify kernel registered**:
```bash
jupyter kernelspec list

# Output:
# Available kernels:
#   python3          /Users/name/anaconda3/share/jupyter/kernels/python3
#   scrna_project    /Users/name/Library/Jupyter/kernels/scrna_project
```

### Alternative Fix: Install in Current Environment

If you don't want separate environments:

```bash
# Activate current environment (base):
micromamba activate base

# Install scanpy:
pip install scanpy

# Restart kernel in notebook
```

**Trade-offs**:
- ✅ Quick, no kernel management needed
- ❌ Pollutes base environment
- ❌ Version conflicts with other projects
- ❌ Not reproducible for others

---

## Phase 4: Verify

### Verification Test 1: Switch Kernel in Notebook

**In Jupyter Lab**:
1. Click "Kernel" menu → "Change kernel..." → Select "Python (scrna_project)"
2. Restart kernel
3. Run Cell 1

```python
# Cell 1:
import scanpy as sc
import numpy as np
import pandas as pd

print(f"Scanpy version: {sc.__version__}")
print(f"Python executable: {sys.executable}")

# Output:
# Scanpy version: 1.9.3 ✓
# Python executable: /Users/name/anaconda3/envs/scrna_project/bin/python ✓
```

### Verification Test 2: All Dependencies Work

```python
# Cell 2: Test key dependencies
import scanpy as sc
import leidenalg  # Clustering algorithm
import igraph  # Graph library for Leiden

print("All dependencies loaded successfully ✓")

# Test basic functionality:
import anndata as ad
adata = ad.AnnData(np.random.rand(100, 50))
print(f"Created test AnnData: {adata.shape} ✓")
```

### Verification Test 3: Environment Reproducible

**Document environment** (critical for reproducibility):

```python
# Cell 3: Document environment
import sys
import scanpy as sc

print("## Environment Info")
print(f"Python: {sys.version}")
print(f"Scanpy: {sc.__version__}")

# Generate requirements file:
!pip freeze > requirements.txt
print("\n✓ requirements.txt created")
```

**Check requirements.txt**:
```
anndata==0.8.0
leidenalg==0.9.1
numpy==1.24.3
pandas==2.0.1
scanpy==1.9.3
...
```

---

## Phase 5: Document

### Updated Cell 0 (Environment Check)

```python
# Cell 0: Environment Setup and Verification
"""
Environment Requirements
------------------------
This notebook requires the 'scrna_project' micromamba environment.

Setup:
    micromamba env create -f environment.yml
    micromamba activate scrna_project
    python -m ipykernel install --user --name=scrna_project

In Jupyter: Kernel → Change kernel → Python (scrna_project)
"""

import sys
print(f"Python: {sys.version}")
print(f"Executable: {sys.executable}")

# Verify we're in correct environment:
assert 'scrna_project' in sys.executable, \
    f"Wrong environment! Using {sys.executable}. Switch kernel to 'Python (scrna_project)'"

print("✓ Correct environment loaded")
```

### Added Markdown Cell (Setup Instructions)

```markdown
## Setup Instructions

### First-time setup:

1. **Create environment**:
   ```bash
   micromamba env create -f environment.yml
   ```

2. **Register Jupyter kernel**:
   ```bash
   micromamba activate scrna_project
   python -m ipykernel install --user --name=scrna_project --display-name="Python (scrna_project)"
   ```

3. **Launch Jupyter**:
   ```bash
   jupyter lab
   ```

4. **Select kernel**: Kernel → Change kernel → Python (scrna_project)

### Troubleshooting:

- **"ModuleNotFoundError"**: Check kernel (top right corner shows which kernel is active)
- **Wrong kernel listed**: Re-run kernel registration command
- **Kernel not starting**: `micromamba activate scrna_project && jupyter lab`
```

### Created environment.yml

```yaml
# environment.yml
name: scrna_project
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - pip
  - pip:
    - scanpy==1.9.3
    - anndata==0.8.0
    - leidenalg==0.9.1
    - python-igraph==0.10.4
    - numpy==1.24.3
    - pandas==2.0.1
    - matplotlib==3.7.1
    - seaborn==0.12.2
    - jupyter
    - ipykernel
```

**Usage**:
```bash
# Others can recreate exact environment:
micromamba env create -f environment.yml
micromamba activate scrna_project
python -m ipykernel install --user --name=scrna_project
```

---

## Summary

| Phase | Action | Outcome |
|-------|--------|---------|
| Diagnose | Checked Python executable, verified package location | Using base env, scanpy in project env |
| Isolate | Listed Jupyter kernels | Project environment kernel not registered |
| Fix | Registered project environment as Jupyter kernel | Kernel available in Jupyter |
| Verify | Switched kernel, tested imports, documented environment | All imports work, reproducible ✓ |
| Document | Added environment check, setup instructions, environment.yml | Others can reproduce ✓ |

---

## Key Lessons

1. **Jupyter kernels != micromamba environments**: Installing a package in an environment doesn't automatically make it available to Jupyter
2. **Always register kernels**: After creating a project environment, register it with `python -m ipykernel install`
3. **Environment assertion**: Add assertion in first cell to catch wrong environment immediately
4. **Document setup**: Include environment.yml so others can reproduce
5. **Check executable path**: `sys.executable` shows which Python Jupyter is actually using

---

## Prevention for Future Notebooks

**Environment management checklist**:
- [ ] Create dedicated micromamba environment for project
- [ ] Install ipykernel in the environment: `micromamba install ipykernel`
- [ ] Register kernel: `python -m ipykernel install --user --name=myproject`
- [ ] Add environment assertion to first cell of notebook
- [ ] Export environment: `micromamba env export > environment.yml`
- [ ] Document setup instructions in notebook markdown

**Automation script** (add to project):
```bash
#!/bin/bash
# setup_notebook_env.sh

ENV_NAME=$1

if [ -z "$ENV_NAME" ]; then
    echo "Usage: ./setup_notebook_env.sh <environment_name>"
    exit 1
fi

echo "Creating environment: $ENV_NAME"
micromamba env create -f environment.yml -n $ENV_NAME

echo "Activating environment"
micromamba activate $ENV_NAME

echo "Registering Jupyter kernel"
python -m ipykernel install --user --name=$ENV_NAME --display-name="Python ($ENV_NAME)"

echo "✓ Setup complete. Select kernel 'Python ($ENV_NAME)' in Jupyter"
```

---

## Related Issues

**Issue: Package installed but import still fails**
- Check: `pip list | grep package-name` (is it installed?)
- Check: `sys.executable` (are you in the right environment?)
- Check: Kernel name in top-right corner of Jupyter

**Issue: Multiple versions of same package**
- Check: `pip list | grep package-name` (which version?)
- Fix: `pip uninstall package-name && pip install package-name==1.2.3`

**Issue: Kernel keeps dying on import**
- Check: Package version compatibility
- Try: `micromamba install package-name` instead of `pip install` (sometimes resolves dependencies better)
- Check: `jupyter lab` terminal output for error messages
