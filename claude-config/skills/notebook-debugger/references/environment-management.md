# Environment Management for Jupyter Notebooks

Best practices for managing Python environments with Jupyter.

---

## Core Concepts

**Environment**: Isolated Python installation with specific package versions

**Kernel**: Jupyter's connection to a Python interpreter

**Key principle**: One project = one environment = one kernel

---

## Conda vs Pip

### When to use Conda

✅ Scientific computing (numpy, scipy, pandas)
✅ Complex dependencies (requires compiled libraries)
✅ Cross-platform consistency needed
✅ Managing multiple Python versions

### When to use Pip

✅ Pure Python packages
✅ Packages not on conda
✅ Lightweight environments
✅ Integration with requirements.txt

### Best practice: Hybrid approach

```bash
# Create environment with conda:
conda create -n myproject python=3.11

# Install scientific packages with conda:
conda install numpy pandas matplotlib

# Install specialized packages with pip:
pip install scanpy pydeseq2
```

---

## Creating Environments

### With Conda

```bash
# Basic environment:
conda create -n myproject python=3.11

# With packages:
conda create -n myproject python=3.11 numpy pandas jupyter

# From file:
conda env create -f environment.yml

# Clone existing:
conda create --name newenv --clone oldenv
```

### With Virtualenv (pip)

```bash
# Create:
python3 -m venv myproject_env

# Activate:
source myproject_env/bin/activate  # macOS/Linux
myproject_env\Scripts\activate  # Windows

# Install packages:
pip install jupyter numpy pandas

# From file:
pip install -r requirements.txt
```

---

## Registering Kernels

### Register Conda Environment

```bash
# Activate environment:
conda activate myproject

# Install ipykernel:
conda install ipykernel
# Or: pip install ipykernel

# Register kernel:
python -m ipykernel install --user --name=myproject --display-name="Python (myproject)"
```

### Register Virtualenv

```bash
# Activate virtualenv:
source myproject_env/bin/activate

# Install ipykernel:
pip install ipykernel

# Register kernel:
python -m ipykernel install --user --name=myproject --display-name="Python (myproject)"
```

### Managing Kernels

```bash
# List installed kernels:
jupyter kernelspec list

# Remove kernel:
jupyter kernelspec uninstall myproject

# Kernel location (macOS/Linux):
~/.local/share/jupyter/kernels/

# Kernel location (Windows):
%APPDATA%\jupyter\kernels\
```

---

## Exporting Environments

### Conda: environment.yml

```bash
# Export full environment:
conda env export > environment.yml

# Export without builds (more portable):
conda env export --no-builds > environment.yml

# Export only explicit packages:
conda env export --from-history > environment.yml
```

**Example environment.yml**:
```yaml
name: myproject
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - numpy=1.24.3
  - pandas=2.0.1
  - jupyter
  - pip
  - pip:
    - scanpy==1.9.3
    - pydeseq2==0.4.0
```

### Pip: requirements.txt

```bash
# Export all packages:
pip freeze > requirements.txt

# Export only top-level (cleaner):
pip list --format=freeze > requirements.txt

# With specific versions:
pip freeze | grep -E "numpy|pandas|scanpy" > requirements.txt
```

**Example requirements.txt**:
```
numpy==1.24.3
pandas==2.0.1
matplotlib==3.7.1
scanpy==1.9.3
jupyter
```

---

## Reproducibility Best Practices

### Principle: Explicit is better than implicit

**Good**:
```yaml
# environment.yml
dependencies:
  - python=3.11.2
  - numpy=1.24.3
  - pandas=2.0.1
```

**Bad**:
```yaml
# environment.yml
dependencies:
  - python>=3.10
  - numpy
  - pandas
```

### Pin Critical Packages

```
# Critical for reproducibility:
numpy==1.24.3  # Exact version

# Less critical (syntax tools):
black>=23.0.0  # Minimum version OK
```

### Include Python Version

```yaml
# environment.yml
name: myproject
dependencies:
  - python=3.11.2  # Specific Python version
```

```
# requirements.txt
# Python 3.11.2
numpy==1.24.3
```

---

## Common Workflows

### Workflow 1: Start New Project

```bash
# 1. Create environment:
conda create -n myproject python=3.11
conda activate myproject

# 2. Install packages:
conda install jupyter numpy pandas matplotlib
pip install scanpy

# 3. Register kernel:
python -m ipykernel install --user --name=myproject

# 4. Export environment:
conda env export --no-builds > environment.yml

# 5. Launch Jupyter:
jupyter lab

# 6. In notebook: Select kernel "Python (myproject)"
```

### Workflow 2: Reproduce Existing Project

```bash
# With conda:
conda env create -f environment.yml
conda activate myproject
python -m ipykernel install --user --name=myproject
jupyter lab

# With pip:
python3 -m venv myproject_env
source myproject_env/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name=myproject
jupyter lab
```

### Workflow 3: Update Environment

```bash
# Add new package:
conda activate myproject
conda install new-package  # Or: pip install new-package

# Update environment file:
conda env export --no-builds > environment.yml

# Commit changes:
git add environment.yml
git commit -m "Add new-package dependency"
```

### Workflow 4: Sync Environment

```bash
# Update from environment.yml:
conda env update -f environment.yml --prune

# Or recreate:
conda deactivate
conda env remove -n myproject
conda env create -f environment.yml
```

---

## Troubleshooting

### "ModuleNotFoundError" but package is installed

**Cause**: Wrong kernel selected

**Solution**:
```python
# Check which Python Jupyter is using:
import sys
print(sys.executable)

# Should show: /path/to/myproject/bin/python
# If not, switch kernel in Jupyter
```

### Kernel won't start

**Cause**: Missing ipykernel

**Solution**:
```bash
conda activate myproject
conda install ipykernel
python -m ipykernel install --user --name=myproject
```

### Package version conflict

**Cause**: Incompatible package versions

**Solution with conda**:
```bash
# Let conda resolve:
conda install package-a package-b

# If fails, check constraints:
conda search package-a
```

**Solution with pip**:
```bash
# Check dependencies:
pip show package-a

# Downgrade if needed:
pip install package-a==1.0.0
```

### "Kernel appears to have died"

**Cause**: Package incompatibility or missing dependency

**Diagnostic**:
```bash
# Check kernel logs:
jupyter --paths  # Find runtime dir
cat /path/to/runtime/kernel-*.json
```

**Solution**:
```bash
# Reinstall kernel:
jupyter kernelspec uninstall myproject
conda activate myproject
python -m ipykernel install --user --name=myproject
```

---

## Advanced Topics

### Multiple Python Versions

```bash
# Python 3.10 environment:
conda create -n py310 python=3.10 jupyter
conda activate py310
python -m ipykernel install --user --name=py310

# Python 3.11 environment:
conda create -n py311 python=3.11 jupyter
conda activate py311
python -m ipykernel install --user --name=py311

# Both kernels available in Jupyter
```

### Minimal Environments

```bash
# Only essential packages:
conda create -n minimal python=3.11 jupyter ipykernel
conda activate minimal
pip install numpy pandas  # Only what you need

# Benefits: Faster, fewer conflicts, clearer dependencies
```

### Shared Environments

```bash
# System-wide kernel (requires sudo):
python -m ipykernel install --name=shared

# Team environment on shared server:
conda env create -f environment.yml --prefix /shared/envs/myproject
python -m ipykernel install --prefix=/shared/envs/myproject --name=team_project
```

---

## Docker for Ultimate Reproducibility

### Why Docker?

- Includes OS-level dependencies
- 100% reproducible across machines
- Includes system libraries (gcc, fortran, etc.)

### Example Dockerfile

```dockerfile
FROM python:3.11-slim

# Install system dependencies:
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy environment file:
COPY requirements.txt .

# Install Python packages:
RUN pip install --no-cache-dir -r requirements.txt

# Install Jupyter:
RUN pip install jupyter

# Create workspace:
WORKDIR /workspace

# Launch Jupyter:
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--allow-root", "--no-browser"]
```

### Build and run:

```bash
# Build image:
docker build -t myproject .

# Run container:
docker run -p 8888:8888 -v $(pwd):/workspace myproject
```

---

## Automation Scripts

### setup_environment.sh

```bash
#!/bin/bash
# Setup script for new project

ENV_NAME=$1

if [ -z "$ENV_NAME" ]; then
    echo "Usage: ./setup_environment.sh <environment_name>"
    exit 1
fi

echo "Creating conda environment: $ENV_NAME"
conda env create -f environment.yml -n $ENV_NAME

echo "Activating environment"
conda activate $ENV_NAME

echo "Registering Jupyter kernel"
python -m ipykernel install --user --name=$ENV_NAME --display-name="Python ($ENV_NAME)"

echo "✓ Environment ready. Launch Jupyter with: jupyter lab"
```

### update_requirements.sh

```bash
#!/bin/bash
# Update requirements after installing new packages

if [ -f environment.yml ]; then
    echo "Updating environment.yml..."
    conda env export --no-builds > environment.yml
fi

if [ -f requirements.txt ]; then
    echo "Updating requirements.txt..."
    pip freeze > requirements.txt
fi

echo "✓ Requirements updated"
```

---

## Environment Testing

### Verify reproducibility:

```bash
# 1. Export environment:
conda env export --no-builds > environment.yml

# 2. Create test environment:
conda env create -f environment.yml -n test_env

# 3. Run notebook in test environment:
conda activate test_env
python -m ipykernel install --user --name=test_env
jupyter lab

# 4. In notebook: Switch to test_env kernel, run all cells

# 5. Clean up:
conda deactivate
conda env remove -n test_env
jupyter kernelspec uninstall test_env
```

---

## Best Practices Checklist

- [ ] One environment per project
- [ ] Pin critical package versions
- [ ] Include Python version in environment file
- [ ] Register environment as Jupyter kernel
- [ ] Export environment file after adding packages
- [ ] Commit environment.yml or requirements.txt to git
- [ ] Test reproducibility on clean environment
- [ ] Document setup steps in README or notebook
- [ ] Use conda for scientific packages, pip for others
- [ ] Remove unused environments regularly

---

## Quick Reference

```bash
# Conda basics:
conda create -n NAME python=VERSION
conda activate NAME
conda install PACKAGE
conda env export > environment.yml
conda env create -f environment.yml

# Kernel management:
python -m ipykernel install --user --name=NAME
jupyter kernelspec list
jupyter kernelspec uninstall NAME

# Pip basics:
python3 -m venv NAME
source NAME/bin/activate
pip install PACKAGE
pip freeze > requirements.txt
pip install -r requirements.txt
```

---

## Resources

- [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/)
- [Jupyter Kernels](https://jupyter-client.readthedocs.io/en/stable/kernels.html)
- [Python Packaging Guide](https://packaging.python.org/en/latest/)
- [IPykernel Documentation](https://ipykernel.readthedocs.io/en/stable/)
