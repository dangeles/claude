# Debugging Tools Reference

Quick guide to debugging tools for Python, bioinformatics, and system-level issues.

---

## Built-in Python Debugging

### print() Debugging

**When to use**: Quick inspection, simple bugs

```python
# Basic:
print(f"Value of x: {x}")

# Debug syntax (Python 3.8+):
print(f"{x = }")  # Prints "x = 42"

# Multiple variables:
print(f"{x = }, {y = }, {x + y = }")

# Type inspection:
print(f"{type(x) = }, {isinstance(x, int) = }")
```

**Pros**: Fast, no setup, works everywhere
**Cons**: Clutters code, manual cleanup, limited info

---

### pdb (Python Debugger)

**When to use**: Interactive debugging, stepping through code

```python
# Add breakpoint:
import pdb; pdb.set_trace()  # Python <3.7
breakpoint()  # Python 3.7+

# Common commands once in pdb:
# n (next) - Execute next line
# s (step) - Step into function
# c (continue) - Continue until next breakpoint
# l (list) - Show code context
# p variable - Print variable
# pp variable - Pretty print variable
# w (where) - Show stack trace
# u (up) - Move up stack
# d (down) - Move down stack
# q (quit) - Exit debugger
```

**Example session**:
```python
def calculate(x, y):
    breakpoint()  # Debugger stops here
    result = x / y
    return result

calculate(10, 0)  # Will trigger pdb before division by zero
# In pdb:
# (Pdb) p y
# 0
# (Pdb) # Ah! Division by zero will happen next line
```

**Pros**: Interactive, full control, inspect state
**Cons**: Slower workflow, terminal-based interface

---

### logging Module

**When to use**: Production code, long-running processes, filtering by severity

```python
import logging

# Configure logging:
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename='debug.log'
)

# Use in code:
logging.debug("Detailed info for debugging")
logging.info("General information")
logging.warning("Something unexpected")
logging.error("Error occurred")
logging.critical("Critical failure")

# In production, set level=INFO to hide DEBUG messages
```

**Pros**: Configurable, persistent, different severity levels
**Cons**: Requires setup, harder to parse than print

---

## Memory Profiling

### memory_profiler

**When to use**: Tracking memory usage line-by-line

```bash
micromamba install memory_profiler
```

```python
from memory_profiler import profile

@profile
def my_function():
    a = [1] * (10 ** 6)  # Allocate 1M integers
    b = [2] * (2 * 10 ** 7)  # Allocate 20M integers
    del b
    return a

my_function()
```

```bash
python -m memory_profiler script.py
```

**Output**:
```
Line #    Mem usage    Increment   Line Contents
================================================
     3   38.816 MiB   38.816 MiB   @profile
     4                             def my_function():
     5   46.492 MiB    7.676 MiB       a = [1] * (10 ** 6)
     6  199.117 MiB  152.625 MiB       b = [2] * (2 * 10 ** 7)
     7   46.629 MiB -152.488 MiB       del b
     8   46.629 MiB    0.000 MiB       return a
```

**Pros**: Line-by-line memory tracking, spots leaks
**Cons**: Slow (~10x overhead)

---

### tracemalloc (Built-in)

**When to use**: Finding memory leaks, tracking allocations

```python
import tracemalloc

tracemalloc.start()

# ... your code that might leak memory ...

snapshot = tracemalloc.take_snapshot()
top_stats = snapshot.statistics('lineno')

print("Top 10 memory allocations:")
for stat in top_stats[:10]:
    print(stat)

# Output shows:
# /path/to/file.py:42: 152.6 MiB
# /path/to/file.py:58: 45.3 MiB
```

**Pros**: Built-in, shows allocation locations
**Cons**: Performance overhead during profiling

---

## Performance Profiling

### cProfile (Built-in)

**When to use**: Finding slow functions

```bash
python -m cProfile -o output.prof script.py
```

```python
# Or in code:
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()

# ... code to profile ...

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats('cumulative')
stats.print_stats(10)  # Top 10 slowest
```

**Pros**: Built-in, shows time per function
**Cons**: Function-level only (not line-level)

---

### line_profiler

**When to use**: Line-by-line timing (like memory_profiler but for time)

```bash
micromamba install line_profiler
```

```python
@profile
def slow_function():
    total = 0
    for i in range(10000):
        total += i ** 2  # Is this slow?
    return total
```

```bash
kernprof -l -v script.py
```

**Pros**: Shows exactly which lines are slow
**Cons**: Requires decorator, adds overhead

---

## Bioinformatics-Specific Tools

### scanpy/anndata Inspection

```python
import scanpy as sc

# Check AnnData structure:
print(adata)  # Shows obs, var, uns, layers

# Check data types:
print(adata.X)  # Sparse or dense?
print(type(adata.X))  # scipy.sparse.csr_matrix?

# Check dimensions:
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# Check for NaN/inf:
import numpy as np
print(f"NaN in X: {np.isnan(adata.X.data).any()}")  # For sparse
print(f"Inf in X: {np.isinf(adata.X.data).any()}")

# Check gene names:
print(adata.var_names[:10])  # First 10 genes

# Check obs/var columns:
print(adata.obs.columns.tolist())
print(adata.var.columns.tolist())
```

### pandas Inspection

```python
import pandas as pd

# Check DataFrame structure:
df.info()  # Dtypes, non-null counts, memory usage

# Check for issues:
print(f"Duplicates: {df.duplicated().sum()}")
print(f"NaN values:\n{df.isna().sum()}")

# Check data types (common source of bugs):
print(df.dtypes)
# Look for 'object' dtype where you expect numeric

# Check unique values (spot encoding issues):
for col in df.columns:
    print(f"{col}: {df[col].nunique()} unique values")
    if df[col].nunique() < 10:
        print(f"  Values: {df[col].unique()}")

# Memory usage:
print(f"Total memory: {df.memory_usage(deep=True).sum() / 1e6:.2f} MB")
```

---

## System-Level Debugging

### strace (Linux - System call tracing)

**When to use**: File access issues, network issues, mysterious hangs

```bash
# Trace all system calls:
strace python script.py

# Trace only file operations:
strace -e trace=file python script.py

# Trace only network operations:
strace -e trace=network python script.py

# Save to file:
strace -o trace.log python script.py
```

**Example output**:
```
open("/path/to/file.csv", O_RDONLY) = 3  # Successfully opened, file descriptor 3
read(3, "gene,expression\nGeneA,5.2\n", 4096) = 45  # Read 45 bytes
```

**Useful for**:
- "Why can't it find this file?" → See what paths are actually being opened
- "Why is this hanging?" → See what system call it's stuck on

---

### lsof (List Open Files)

**When to use**: "Too many open files" error, file locks

```bash
# List all open files by process:
lsof -p <pid>

# List what's using a specific file:
lsof /path/to/file.csv

# List all open files by Python:
lsof -c python
```

**Common bug**: Opening files in loop without closing
```python
# Bug:
for file in files:
    f = open(file)  # Never closed! Runs out of file descriptors
    process(f)

# Fix:
for file in files:
    with open(file) as f:  # Automatically closed
        process(f)
```

---

## Network Debugging

### tcpdump / Wireshark

**When to use**: API calls failing, network issues

```bash
# Capture HTTP traffic:
sudo tcpdump -i any -s 0 -A 'tcp port 80'

# Capture to file for Wireshark:
sudo tcpdump -i any -s 0 -w capture.pcap 'tcp port 443'
```

**Useful for**:
- "Is the request actually being sent?"
- "What response is the server sending?"

---

## Debugging Distributed/Parallel Code

### logging in multiprocessing

```python
import logging
import multiprocessing

def worker(x):
    logger = logging.getLogger(f"worker-{multiprocessing.current_process().name}")
    logger.info(f"Processing {x}")
    return x ** 2

# Configure logging:
logging.basicConfig(level=logging.INFO)

# Use pool:
with multiprocessing.Pool() as pool:
    results = pool.map(worker, range(10))
```

**Tip**: Include process ID/name in logs to track which worker did what.

---

## Debugging Jupyter Notebooks

### %debug Magic

```python
# After error in cell:
%debug

# Drops into pdb at the point of error
# Can inspect variables, stack trace
```

### %pdb Magic

```python
# Auto-enter debugger on exception:
%pdb on

# Now any error automatically starts pdb
```

### %%time Magic

```python
%%time
# Code in this cell
slow_operation()

# Prints execution time
```

---

## Remote Debugging

### pdb over SSH

```python
# On remote server:
import pdb
pdb.set_trace()  # Will wait for you to connect

# Connect from local:
ssh -t user@server python script.py
# Now in pdb on remote server
```

---

## Visual Debugging Tools

### objgraph (Object Reference Visualization)

**When to use**: Finding reference cycles, understanding object relationships

```bash
micromamba install objgraph
```

```python
import objgraph

# Show most common types:
objgraph.show_most_common_types()

# Show what references an object:
objgraph.show_refs([my_object], filename='refs.png')

# Find reference cycles:
objgraph.show_chain(
    objgraph.find_backref_chain(my_object, inspect.isclass),
    filename='chain.png'
)
```

**Useful for**: Memory leak investigations

---

## Summary: When to Use Which Tool

| Problem | Tool | Why |
|---------|------|-----|
| Quick inspection | print() | Fast, simple |
| Step through code | pdb / breakpoint() | Interactive control |
| Production logging | logging module | Persistent, filterable |
| Memory leak | tracemalloc, memory_profiler | Shows allocations |
| Slow code | cProfile, line_profiler | Shows timing |
| File access issues | strace, lsof | System-level view |
| API/network issues | tcpdump, Wireshark | See actual traffic |
| Jupyter notebooks | %debug, %pdb | Notebook-specific |
| Reference cycles | objgraph | Visual object graphs |
| Pandas data issues | df.info(), df.isna() | Data quality checks |
