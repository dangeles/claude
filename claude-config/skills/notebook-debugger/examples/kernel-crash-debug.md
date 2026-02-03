# Kernel Crash Debug Example: Memory Error in Data Processing

## Problem Report

**User complaint**: "Notebook kernel crashes when running cell 15. No error message, just dies."

**Notebook**: `rna_seq_analysis.ipynb`
**Cell that crashes**: Cell 15 (data aggregation)
**Environment**: Jupyter Lab 3.6, Python 3.11, 16GB RAM machine

---

## Phase 1: Diagnose

### Initial Investigation

**Test 1: Can we reproduce?**
```
Action: Kernel → Restart & Run All
Result: Kernel crashes at cell 15 ✓ Reproducible
```

**Test 2: What's in cell 15?**
```python
# Cell 15:
# Aggregate expression by cell type
expression_by_celltype = {}
for celltype in adata.obs['celltype'].unique():
    mask = adata.obs['celltype'] == celltype
    subset = adata[mask, :]
    expression_by_celltype[celltype] = subset.X.mean(axis=0)
```

**Test 3: Check data size**
```python
# Cell 14 (before crash):
print(f"adata shape: {adata.shape}")
print(f"Memory usage: {adata.X.data.nbytes / 1e9:.2f} GB")

# Output:
# adata shape: (50000, 20000)  # 50k cells, 20k genes
# Memory usage: 4.2 GB (sparse matrix)
```

**Observation**: 50,000 cells, 20 cell types → ~2,500 cells per type. Creating dense arrays for each subset.

---

## Phase 2: Isolate

### Hypothesis: Memory Error During Subset Operations

**Test: Run cell 15 with smaller data**
```python
# Test with first 1000 cells only:
adata_small = adata[:1000, :]
expression_by_celltype = {}
for celltype in adata_small.obs['celltype'].unique():
    mask = adata_small.obs['celltype'] == celltype
    subset = adata_small[mask, :]
    expression_by_celltype[celltype] = subset.X.mean(axis=0)

# Result: ✓ Works! No crash with 1000 cells
```

**Conclusion**: Memory-related, triggered by full dataset size.

### Memory Profiling

**Add memory tracking**:
```python
import sys

expression_by_celltype = {}
for i, celltype in enumerate(adata.obs['celltype'].unique()):
    mask = adata.obs['celltype'] == celltype
    subset = adata[mask, :]

    # Track memory:
    print(f"Cell type {i+1}/20: {celltype}")
    print(f"  Subset size: {sys.getsizeof(subset.X) / 1e9:.2f} GB")

    expression_by_celltype[celltype] = subset.X.mean(axis=0)

# Output (before crash):
# Cell type 1/20: T_cell
#   Subset size: 0.4 GB
# Cell type 2/20: B_cell
#   Subset size: 0.3 GB
# ...
# Cell type 8/20: Macrophage
#   Subset size: 0.5 GB
# [KERNEL CRASHED]
```

**Observation**: Accumulated 8 subsets × ~0.4GB = ~3.2GB, plus original 4.2GB = 7.4GB. Close to system limit.

### Root Cause Identified

**Problem**: Creating dense subsets from sparse matrix, accumulating in memory

```python
# Each iteration:
subset = adata[mask, :]  # subset.X is sparse (good)
subset.X.mean(axis=0)  # .mean() converts to dense (bad!)

# Dense array size per cell type:
# (2500 cells) × (20000 genes) × (8 bytes/float64) = 400 MB
# 20 cell types × 400 MB = 8 GB total
# Original data (4.2 GB) + dense subsets (8 GB) = 12.2 GB → exceeds 16GB RAM
```

---

## Phase 3: Fix

### Solution: Keep Operations Sparse

```python
# Original (crashes):
expression_by_celltype = {}
for celltype in adata.obs['celltype'].unique():
    mask = adata.obs['celltype'] == celltype
    subset = adata[mask, :]
    expression_by_celltype[celltype] = subset.X.mean(axis=0)  # Dense!

# Fixed (stays sparse):
import numpy as np
import scipy.sparse

expression_by_celltype = {}
for celltype in adata.obs['celltype'].unique():
    mask = adata.obs['celltype'] == celltype
    subset = adata[mask, :]

    # Keep sparse:
    if scipy.sparse.issparse(subset.X):
        mean_expr = np.array(subset.X.mean(axis=0)).flatten()
    else:
        mean_expr = subset.X.mean(axis=0)

    expression_by_celltype[celltype] = mean_expr
```

**Memory improvement**:
- Before: 12.2 GB peak (crashes)
- After: 4.5 GB peak (works!)

### Alternative Solution: Process One at a Time

```python
# If keeping sparse doesn't help, process and clear:
expression_by_celltype = {}
for celltype in adata.obs['celltype'].unique():
    mask = adata.obs['celltype'] == celltype
    subset = adata[mask, :]
    expression_by_celltype[celltype] = subset.X.mean(axis=0)

    del subset, mask  # Explicit cleanup
    import gc; gc.collect()  # Force garbage collection
```

---

## Phase 4: Verify

### Verification Test 1: Full Dataset

```python
# Restart kernel, run all cells with fix
# Cell 15 now succeeds ✓

# Check results:
print(f"Cell types processed: {len(expression_by_celltype)}")
print(f"Expression shape: {expression_by_celltype['T_cell'].shape}")

# Output:
# Cell types processed: 20 ✓
# Expression shape: (20000,) ✓
```

### Verification Test 2: Memory Stable

```python
import tracemalloc
tracemalloc.start()

# Run cell 15
expression_by_celltype = {}
for celltype in adata.obs['celltype'].unique():
    # ... fixed code ...

snapshot = tracemalloc.take_snapshot()
top_stats = snapshot.statistics('lineno')
print(f"Peak memory: {snapshot.statistics('lineno')[0].size / 1e9:.2f} GB")

# Output: Peak memory: 4.8 GB ✓ (well below 16GB limit)
```

### Verification Test 3: Results Correct

```python
# Compare with reference implementation:
reference = {}
for celltype in ['T_cell', 'B_cell']:  # Test subset
    mask = adata.obs['celltype'] == celltype
    reference[celltype] = adata[mask, :].X.mean(axis=0)

# Check our fixed version matches:
for celltype in reference:
    diff = np.abs(expression_by_celltype[celltype] - reference[celltype]).max()
    print(f"{celltype}: max diff = {diff}")
    assert diff < 1e-10, "Results don't match!"

# Output:
# T_cell: max diff = 0.0 ✓
# B_cell: max diff = 0.0 ✓
```

---

## Phase 5: Document

### Updated Cell 15 (with documentation)

```python
# Cell 15: Aggregate expression by cell type
# NOTE: Keep operations sparse to avoid memory error
# With 50k cells × 20k genes, dense conversion would require 12+ GB

import numpy as np
import scipy.sparse

expression_by_celltype = {}
for celltype in adata.obs['celltype'].unique():
    mask = adata.obs['celltype'] == celltype
    subset = adata[mask, :]

    # Keep sparse (critical for large datasets):
    if scipy.sparse.issparse(subset.X):
        mean_expr = np.array(subset.X.mean(axis=0)).flatten()
    else:
        mean_expr = subset.X.mean(axis=0)

    expression_by_celltype[celltype] = mean_expr

print(f"Processed {len(expression_by_celltype)} cell types")
```

### Added Markdown Cell (Known Issues)

```markdown
## Memory Requirements

**This notebook requires ~8GB RAM minimum**

**Cell 15** (aggregation step) can crash kernel if insufficient memory:
- Expected memory usage: 4-5 GB peak
- If kernel crashes: Reduce dataset size using `adata = adata[:10000, :]` before cell 15

**Sparse matrix operations**:
- Data is stored as sparse matrix to save memory (4.2 GB vs 8 GB dense)
- Operations that convert to dense (`.mean()`, `.std()`) must be handled carefully
- Current implementation keeps data sparse throughout
```

### Environment Documentation

```python
# Added to cell 1:
"""
Environment Requirements
------------------------
numpy==1.24.3
pandas==2.0.1
scanpy==1.9.3
scipy==1.10.1

Install: pip install -r requirements.txt

Memory: Minimum 8GB RAM recommended
Runtime: ~5 minutes on 16GB machine
"""

import sys
print(f"Python: {sys.version}")
print(f"Scanpy: {sc.__version__}")
```

---

## Summary

| Phase | Action | Outcome |
|-------|--------|---------|
| Diagnose | Reproduced crash, identified cell 15 | Consistent crash at aggregation step |
| Isolate | Tested with smaller data, added memory profiling | Memory accumulation during loop |
| Root Cause | Sparse → dense conversion in `.mean()` | 12GB peak exceeds 16GB RAM |
| Fix | Keep operations sparse | Peak reduced to 4.8GB |
| Verify | Full dataset test, memory monitoring | Notebook runs successfully ✓ |
| Document | Added comments, known issues, requirements | Reproducible for others ✓ |

---

## Key Lessons

1. **Sparse matrices are memory-efficient but fragile**: One operation can convert to dense
2. **Monitor memory in loops**: Accumulation problems not obvious from single iteration
3. **Test with full dataset**: Small test data may hide memory issues
4. **Document memory requirements**: Future users need to know system requirements
5. **Explicit garbage collection helps**: `del` + `gc.collect()` in tight loops

---

## Prevention for Future Notebooks

**Memory-aware development checklist**:
- [ ] Profile memory usage during development: `import tracemalloc`
- [ ] Test with full dataset before declaring "done"
- [ ] Keep sparse operations sparse: Check `.issparse()` before operations
- [ ] Document memory requirements in notebook header
- [ ] Add assertions for sparse data: `assert scipy.sparse.issparse(adata.X)`
