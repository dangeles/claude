# Jupyter Troubleshooting Guide

Quick reference for common Jupyter notebook issues and solutions.

---

## Kernel Issues

### Kernel Won't Start

**Symptoms**: Kernel status stuck on "Starting" or "Connecting"

**Common causes**:
1. Wrong Python environment
2. Missing ipykernel package
3. Port conflict
4. Jupyter server issue

**Solutions**:
```bash
# Check if ipykernel is installed:
pip list | grep ipykernel

# Install if missing:
pip install ipykernel

# Register kernel:
python -m ipykernel install --user --name=myenv

# Restart Jupyter:
jupyter lab stop
jupyter lab
```

---

### Kernel Keeps Dying

**Symptoms**: "Kernel died unexpectedly" message, kernel restarts

**Common causes**:
1. Memory error (most common)
2. Segmentation fault (C extension bug)
3. Infinite loop/recursion
4. Incompatible package versions

**Diagnostic steps**:
```python
# Check memory usage:
import sys
print(f"Object size: {sys.getsizeof(obj) / 1e9:.2f} GB")

# Check for infinite loops:
# Add to suspect cell:
import sys
sys.setrecursionlimit(100)  # Will catch infinite recursion faster

# Check package versions:
import package_name
print(package_name.__version__)
```

**Solutions**:
- **Memory error**: Reduce data size, use chunking, or Dask
- **Segmentation fault**: Update packages, check for C extension bugs
- **Infinite loop**: Add breakpoint() before suspect code
- **Version conflict**: Check requirements, update/downgrade packages

---

### Kernel Not Responding

**Symptoms**: Cell shows `[*]` forever, kernel busy

**Common causes**:
1. Long-running operation
2. Blocked on I/O (network, disk)
3. Deadlock in multithreading
4. Waiting for user input (input() call)

**Solutions**:
```python
# Add timeout to operations:
import signal
signal.alarm(60)  # Timeout after 60 seconds

# Check if waiting for input:
# Look for input() or raw_input() calls in code

# Interrupt kernel:
# Jupyter menu: Kernel → Interrupt
# Or press 'i' twice quickly (JupyterLab)

# If interrupt doesn't work, restart:
# Kernel → Restart
```

---

## Import Errors

### ModuleNotFoundError

**Symptom**: `ModuleNotFoundError: No module named 'X'`

**Diagnostic checklist**:
```python
# 1. Check Python executable:
import sys
print(sys.executable)
# Should show correct environment path

# 2. Check if package is installed:
!pip list | grep package-name

# 3. Check Python path:
print(sys.path)
# Should include package location

# 4. Try importing in terminal:
# Exit Jupyter, activate environment, try:
python -c "import package_name"
```

**Solutions**:
```bash
# Install in correct environment:
pip install package-name

# Or register kernel for project environment:
conda activate myproject
python -m ipykernel install --user --name=myproject

# Switch kernel in Jupyter:
# Kernel → Change kernel → myproject
```

---

### ImportError (package found but import fails)

**Symptom**: `ImportError: cannot import name 'X' from 'Y'`

**Common causes**:
1. Wrong package version
2. Circular import
3. Missing dependencies
4. Name conflict with local file

**Solutions**:
```python
# Check version:
import package_name
print(package_name.__version__)

# Check where package is loaded from:
print(package_name.__file__)

# Upgrade/downgrade:
!pip install package-name==1.2.3

# Check for name conflicts:
import os
print(os.listdir('.'))  # Look for files with same name as package
```

---

## Memory Issues

### MemoryError

**Symptom**: `MemoryError: Unable to allocate X GB`

**Common causes**:
1. Loading entire dataset into memory
2. Creating large arrays
3. Accumulating data in loops
4. Memory leak

**Solutions**:
```python
# Use chunking:
import pandas as pd
for chunk in pd.read_csv('large.csv', chunksize=10000):
    process(chunk)

# Use Dask for out-of-memory computing:
import dask.dataframe as dd
df = dd.read_csv('large.csv')

# Explicit cleanup:
del large_object
import gc
gc.collect()

# Monitor memory:
import tracemalloc
tracemalloc.start()
# ... your code ...
snapshot = tracemalloc.take_snapshot()
for stat in snapshot.statistics('lineno')[:10]:
    print(stat)
```

---

### Notebook File Too Large

**Symptom**: Notebook slow to open, >100MB file size

**Common causes**:
1. Large outputs (printed DataFrames, images)
2. Many print statements in loops
3. Embedded images/plots

**Solutions**:
```bash
# Clear all outputs:
# Edit → Clear All Outputs

# Strip outputs from file:
jupyter nbconvert --clear-output --inplace notebook.ipynb

# Or use nbstripout (automatic):
pip install nbstripout
nbstripout notebook.ipynb

# Add to git hooks:
nbstripout --install
```

**Prevent**:
```python
# Limit output size:
pd.set_option('display.max_rows', 20)

# Don't print in loops:
# Bad:
for i in range(1000):
    print(df)  # Huge notebook!

# Good:
for i in range(1000):
    if i % 100 == 0:
        print(f"Progress: {i}/1000")
```

---

## Cell Execution Issues

### Cell Never Completes

**Symptom**: Cell shows `[*]` forever

**Diagnostic steps**:
1. Check if CPU usage is high (still running) or zero (stuck)
2. Check if disk/network I/O is active
3. Interrupt kernel and check traceback

**Solutions**:
```python
# Add progress tracking:
from tqdm import tqdm
for item in tqdm(large_list):
    process(item)

# Add timeouts:
import signal

def timeout_handler(signum, frame):
    raise TimeoutError("Operation timed out")

signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(60)  # 60 second timeout
try:
    long_operation()
finally:
    signal.alarm(0)  # Cancel timeout
```

---

### Execution Order Problems

**Symptom**: Works manually, fails on "Restart & Run All"

**Cause**: Cells run out of order during development

**Diagnostic**:
```python
# Check execution numbers:
# Look at [n] next to each cell
# Should be sequential: [1], [2], [3], ...
# If not: [8], [12], [3] → cells run out of order

# Check for undefined variables:
print([var for var in dir() if not var.startswith('_')])
```

**Solutions**:
1. Reorganize cells in logical order
2. Add assertions: `assert 'df' in dir(), "Run Cell 3 first"`
3. Test with "Restart & Run All" after changes
4. Use markdown headers to organize sections

---

## Output Issues

### Output Not Displaying

**Symptom**: Cell runs but no output appears

**Common causes**:
1. Forgot to print/display
2. Output captured by IPython
3. Backend issue (matplotlib)

**Solutions**:
```python
# Explicit display:
from IPython.display import display
display(df)

# For matplotlib:
%matplotlib inline
import matplotlib.pyplot as plt
plt.plot([1, 2, 3])
plt.show()  # Explicit show()

# Check if output is captured:
# Look for IPython.utils.capture in code
```

---

### Plots Not Showing

**Symptom**: Matplotlib code runs but no plot appears

**Solutions**:
```python
# Enable inline backend:
%matplotlib inline

# Or interactive backend:
%matplotlib widget

# Explicit show():
import matplotlib.pyplot as plt
plt.plot([1, 2, 3])
plt.show()  # Required in scripts, optional in notebooks

# Check backend:
import matplotlib
print(matplotlib.get_backend())
# Should be: 'module://ipykernel.pylab.backend_inline'
```

---

## Performance Issues

### Slow Cell Execution

**Symptom**: Cell takes minutes instead of seconds

**Diagnostic**:
```python
# Time cell execution:
%%time
your_code_here()

# Profile line-by-line:
%load_ext line_profiler
%lprun -f function_name function_name()

# Check for loops over DataFrames:
# Bad (slow):
for i in range(len(df)):
    df.loc[i, 'new'] = df.loc[i, 'old'] * 2

# Good (fast):
df['new'] = df['old'] * 2
```

---

### Jupyter Lab Slow to Load

**Solutions**:
```bash
# Clear workspace:
rm -rf ~/.jupyter/lab/workspaces/*

# Disable extensions:
jupyter labextension disable extension-name

# Build clean:
jupyter lab clean
jupyter lab build

# Check for large notebooks:
du -sh *.ipynb | sort -h
```

---

## Connection Issues

### Can't Connect to Kernel

**Symptom**: "Connecting to kernel..." forever

**Solutions**:
```bash
# Check if Jupyter server is running:
jupyter lab list

# Kill zombie kernels:
jupyter kernelspec list
jupyter kernel list  # Shows running kernels
# Kill stuck kernel by PID

# Restart Jupyter:
jupyter lab stop
jupyter lab

# Check port conflicts:
lsof -i :8888  # Check if port 8888 is in use
```

---

### Browser Can't Connect

**Symptom**: "Unable to connect" in browser

**Solutions**:
```bash
# Check Jupyter is running:
jupyter lab list

# Check firewall:
# Allow connections on port 8888

# Use explicit IP:
jupyter lab --ip=127.0.0.1

# Check token:
jupyter lab list  # Shows token
# Copy full URL with token to browser
```

---

## Environment Issues

### Wrong Python Version

**Symptom**: Code fails with version-specific error

**Diagnostic**:
```python
import sys
print(sys.version)
print(sys.executable)
```

**Solutions**:
```bash
# Create environment with specific Python:
conda create -n myenv python=3.11
conda activate myenv
python -m ipykernel install --user --name=myenv

# Switch kernel in notebook:
# Kernel → Change kernel → myenv
```

---

### Package Version Conflicts

**Symptom**: "Package A requires B>=1.0, but you have B==0.9"

**Solutions**:
```bash
# Check dependencies:
pip show package-name

# Update packages:
pip install --upgrade package-name

# Install specific versions:
pip install package-a==1.0 package-b==2.0

# Use conda to resolve:
conda install package-name  # Better dependency resolution
```

---

## Quick Diagnostic Checklist

When notebook fails, check in this order:

1. **Kernel running?** Look at kernel indicator (top right)
2. **Correct kernel?** Check kernel name matches project
3. **Imports work?** Try importing in first cell
4. **Execution order?** Check cell numbers are sequential
5. **Memory OK?** Check Activity Monitor/Task Manager
6. **Outputs cleared?** Try Edit → Clear All Outputs
7. **Fresh kernel?** Try Kernel → Restart & Run All

---

## Emergency Recovery

### Notebook Won't Open

**Symptom**: Browser shows error, notebook won't load

**Solutions**:
```bash
# Check file integrity:
python -m json.tool notebook.ipynb > /dev/null
# If error: JSON is corrupted

# View raw file:
cat notebook.ipynb

# Recover from checkpoint:
ls .ipynb_checkpoints/
cp .ipynb_checkpoints/notebook-checkpoint.ipynb notebook_recovered.ipynb

# Strip outputs and try again:
jupyter nbconvert --clear-output --to notebook --output=recovered.ipynb notebook.ipynb
```

---

### Lost All Work

**Symptom**: Notebook file deleted or corrupted

**Recovery options**:
1. Check `.ipynb_checkpoints/` directory (autosaves)
2. Check git history if versioned
3. Check system trash/recycle bin
4. Check Time Machine or system backups

**Prevention**:
```bash
# Enable autosave:
# Jupyter automatically saves to .ipynb_checkpoints/

# Use version control:
git init
git add notebook.ipynb
git commit -m "Initial commit"

# Backup before major changes:
cp notebook.ipynb notebook_backup_$(date +%Y%m%d).ipynb
```

---

## Useful Jupyter Magics

```python
# Time execution:
%time single_statement()
%%time
# entire cell

# Debug:
%debug  # After error, drop into debugger

# Enable automatic debugging:
%pdb on

# Load external code:
%load script.py

# Run external script:
%run script.py

# List variables:
%who  # Names only
%whos  # With details

# Clear variables:
%reset  # Clear all variables

# System commands:
!ls
!pip install package
```

---

## Prevention Best Practices

1. **Use version control**: `git add notebook.ipynb && git commit`
2. **Test reproducibility**: Regularly run "Restart & Run All"
3. **Document environment**: Keep `requirements.txt` or `environment.yml`
4. **Clear outputs**: Before committing, clear large outputs
5. **Organize cells**: Use markdown headers, logical sections
6. **Add assertions**: Check dependencies exist before use
7. **Monitor memory**: Use chunking for large datasets
8. **Set random seeds**: For reproducible results
