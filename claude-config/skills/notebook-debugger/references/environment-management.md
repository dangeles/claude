# Environment Management for Jupyter Notebooks

Best practices for managing Python environments with Jupyter.

---

## Core Concepts

**Environment**: Isolated Python installation with specific package versions

**Kernel**: Jupyter's connection to a Python interpreter

**Key principle**: One project = one environment = one kernel

---

## micromamba vs Pip

### When to use micromamba

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
# Create environment with micromamba:
micromamba create -n myproject python=3.11

# Install scientific packages with micromamba:
micromamba install numpy pandas matplotlib

# Install specialized packages with pip:
pip install scanpy pydeseq2
```

---

## Creating Environments

### With micromamba

```bash
# Basic environment:
micromamba create -n myproject python=3.11

# With packages:
micromamba create -n myproject python=3.11 numpy pandas jupyter

# From file:
micromamba env create -f environment.yml

# Clone existing:
micromamba create --name newenv --clone oldenv
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

### Register micromamba Environment

```bash
# Activate environment:
micromamba activate myproject

# Install ipykernel:
micromamba install ipykernel

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

### micromamba: environment.yml

```bash
# Export micromamba packages:
micromamba env export > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# Export without builds (more portable):
micromamba env export --no-builds > environment.yml

# Export only explicit packages:
micromamba env export --from-history > environment.yml
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
micromamba create -n myproject python=3.11
micromamba activate myproject

# 2. Install packages:
micromamba install jupyter numpy pandas matplotlib
pip install scanpy

# 3. Register kernel:
python -m ipykernel install --user --name=myproject

# 4. Export environment:
# Export micromamba packages:
micromamba env export --no-builds > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# 5. Launch Jupyter:
jupyter lab

# 6. In notebook: Select kernel "Python (myproject)"
```

### Workflow 2: Reproduce Existing Project

```bash
# With micromamba:
micromamba env create -f environment.yml
micromamba activate myproject
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
micromamba activate myproject
micromamba install new-package  # Or: pip install new-package

# Update environment file:
# Export micromamba packages:
micromamba env export --no-builds > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# Commit changes:
git add environment.yml
git commit -m "Add new-package dependency"
```

### Workflow 4: Sync Environment

```bash
# For incremental updates (updates already-installed packages only):
micromamba env update -f environment.yml

# For full environment recreation:
micromamba env create --yes -f environment.yml

# Or remove and recreate:
micromamba deactivate
micromamba env remove -n myproject
micromamba env create -f environment.yml
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
micromamba activate myproject
micromamba install ipykernel
python -m ipykernel install --user --name=myproject
```

### Package version conflict

**Cause**: Incompatible package versions

**Solution with micromamba**:
```bash
# Let micromamba resolve:
micromamba install package-a package-b

# If fails, check constraints:
micromamba search package-a
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
micromamba activate myproject
python -m ipykernel install --user --name=myproject
```

---

## Advanced Topics

### Multiple Python Versions

```bash
# Python 3.10 environment:
micromamba create -n py310 python=3.10 jupyter
micromamba activate py310
python -m ipykernel install --user --name=py310

# Python 3.11 environment:
micromamba create -n py311 python=3.11 jupyter
micromamba activate py311
python -m ipykernel install --user --name=py311

# Both kernels available in Jupyter
```

### Minimal Environments

```bash
# Only essential packages:
micromamba create -n minimal python=3.11 jupyter ipykernel
micromamba activate minimal
micromamba install numpy pandas  # Only what you need

# Benefits: Faster, fewer conflicts, clearer dependencies
```

### Shared Environments

```bash
# System-wide kernel (requires sudo):
python -m ipykernel install --name=shared

# Team environment on shared server:
micromamba env create -f environment.yml --prefix /shared/envs/myproject
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

echo "Creating micromamba environment: $ENV_NAME"
micromamba env create -f environment.yml -n $ENV_NAME

echo "Activating environment"
micromamba activate $ENV_NAME

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
    micromamba env export --no-builds > environment.yml
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
# Export micromamba packages:
micromamba env export --no-builds > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# 2. Create test environment:
micromamba env create -f environment.yml -n test_env

# 3. Run notebook in test environment:
micromamba activate test_env
python -m ipykernel install --user --name=test_env
jupyter lab

# 4. In notebook: Switch to test_env kernel, run all cells

# 5. Clean up:
micromamba deactivate
micromamba env remove -n test_env
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
- [ ] Use micromamba for scientific packages, pip for others
- [ ] Remove unused environments regularly

---

## Quick Reference

```bash
# micromamba basics:
micromamba create -n NAME python=VERSION
micromamba activate NAME
micromamba install PACKAGE
micromamba env export > environment.yml
micromamba env create -f environment.yml

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

- [micromamba Documentation](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
- [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/) (for reference; micromamba uses compatible commands)
- [Jupyter Kernels](https://jupyter-client.readthedocs.io/en/stable/kernels.html)
- [Python Packaging Guide](https://packaging.python.org/en/latest/)
- [IPykernel Documentation](https://ipykernel.readthedocs.io/en/stable/)
