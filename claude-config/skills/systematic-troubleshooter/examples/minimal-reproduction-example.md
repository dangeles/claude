# Minimal Reproducible Example: How to Create One

## Why Minimal Examples Matter

**Problem**: Debugging in a 10,000-line codebase with 30 dependencies takes hours and often leads to wrong conclusions.

**Solution**: Strip away everything unrelated to the bug. Often, the act of simplifying reveals the root cause immediately.

---

## Example 1: Import Error in Production Pipeline

### Original (Complex, Hard to Debug)

**Context**: 500-line bioinformatics pipeline fails with import error

```python
# pipeline.py (500 lines)
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from src.utils.preprocessing import normalize_counts  # FAILS HERE
from src.utils.plotting import plot_umap
from src.models.clustering import run_leiden
# ... 490 more lines ...
```

**Error**:
```
ModuleNotFoundError: No module named 'src.utils.preprocessing'
```

**Why it's hard**:
- 500 lines to read
- 30+ imports to check
- Complex directory structure to understand
- Might be environment, path, or actual missing file

### Minimal Reproduction (Isolated the Issue)

**Steps to minimize**:
1. Create new file with ONLY the failing import
2. Remove all other code
3. Test in same environment

```python
# test_import.py (2 lines)
from src.utils.preprocessing import normalize_counts
print("Import successful!")
```

**Run**:
```bash
python test_import.py
# ModuleNotFoundError: No module named 'src.utils.preprocessing'
```

**Result**: Same error with 2 lines instead of 500. Now we can focus on the import mechanism.

**Root cause discovered**: File exists at `src/utils/preprocessing.py` but `src/` is not a package (missing `__init__.py`).

**Fix**:
```bash
touch src/__init__.py
touch src/utils/__init__.py
```

---

## Example 2: DataFrame Operation Produces Wrong Results

### Original (Complex)

```python
def analyze_expression(adata, gene_list, normalization='scran', filter_cells=True,
                       min_genes=200, min_cells=3, target_sum=1e4, log_transform=True,
                       highly_variable=True, n_top_genes=2000, batch_key=None):
    """Complex function with 10 parameters and 50 lines of logic"""
    # Filter cells
    if filter_cells:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalize
    if normalization == 'scran':
        # ... 15 lines ...
    elif normalization == 'seurat':
        # ... 10 lines ...

    # More processing...
    # ... 20 more lines ...

    # Final step that gives wrong answer
    result = adata[:, gene_list].X.mean(axis=0)  # BUG SOMEWHERE
    return result

# When called: Returns wrong mean values (off by 10x)
result = analyze_expression(adata, ['CD8A', 'CD4'], normalization='scran',
                           filter_cells=True, highly_variable=True)
```

**Problem**: Function returns values 10x too high. Where's the bug?

### Minimal Reproduction

**Strategy**: Remove everything except the failing operation

```python
import numpy as np
import scanpy as sc

# Create minimal test data (3 cells, 2 genes)
adata = sc.AnnData(np.array([[1, 2],
                              [3, 4],
                              [5, 6]]))
adata.var_names = ['CD8A', 'CD4']

# Just the failing operation
result = adata[:, ['CD8A', 'CD4']].X.mean(axis=0)
print("Result:", result)  # Expect [3, 4]
print("Actual:", result)  # Got [ 3.  4.], correct!

# So the mean() operation works correctly. Bug must be in earlier processing.
# Let's test normalization:

sc.pp.normalize_total(adata, target_sum=1e4)
result = adata[:, ['CD8A', 'CD4']].X.mean(axis=0)
print("After normalization:", result)  # Expect [3e4, 4e4] if input was [3, 4]
print("Actual:", result)  # Got [3e4, 4e4] - so normalization is 10x multiplier!

# AHA! Original data wasn't [1,2,3,4,5,6], it was [0.1, 0.2, ...] (CPM)
# Normalizing CPM data (already normalized) applies normalization again → 10x error
```

**Root cause**: Data was already normalized (CPM), but function normalized again.

**Fix**: Check if data is already normalized before normalizing:
```python
if not adata.uns.get('normalized'):
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.uns['normalized'] = True
```

---

## Example 3: Intermittent Test Failure

### Original (Unpredictable)

```python
def test_clustering():
    """Test that runs leiden clustering and checks result."""
    adata = load_test_data()  # 10,000 cells
    preprocess_pipeline(adata)  # 5 steps
    sc.tl.leiden(adata, resolution=0.5)

    # Fails ~30% of the time
    assert len(adata.obs['leiden'].unique()) == 8  # Expected 8 clusters
```

**Problem**: Test fails intermittently. Sometimes 8 clusters, sometimes 7, sometimes 9.

### Minimal Reproduction

**Strategy**: Isolate the random element

```python
import scanpy as sc
import numpy as np

# Create minimal test case with known structure
np.random.seed(42)  # Set seed
adata = sc.datasets.pbmc68k_reduced()  # Standard test data

sc.tl.leiden(adata, resolution=0.5)  # WITHOUT seed
print("Run 1:", len(adata.obs['leiden'].unique()))  # 8 clusters

sc.tl.leiden(adata, resolution=0.5)  # WITHOUT seed
print("Run 2:", len(adata.obs['leiden'].unique()))  # 9 clusters (different!)

# AHA! leiden() is non-deterministic without random_state
sc.tl.leiden(adata, resolution=0.5, random_state=0)
print("Run 3:", len(adata.obs['leiden'].unique()))  # 8 clusters

sc.tl.leiden(adata, resolution=0.5, random_state=0)
print("Run 4:", len(adata.obs['leiden'].unique()))  # 8 clusters (same!)
```

**Root cause**: Missing `random_state` parameter in leiden clustering.

**Fix**:
```python
sc.tl.leiden(adata, resolution=0.5, random_state=0)  # Deterministic
```

---

## Template for Creating Minimal Reproductions

```python
"""
Minimal Reproducible Example for [Bug Description]

Original issue: [1-sentence description]
Minimal test: [What this file does]
"""

# 1. Minimal imports (only what's needed)
import package_with_bug

# 2. Minimal data (smallest dataset that triggers bug)
data = create_minimal_test_data()

# 3. Minimal code (just the failing operation)
result = operation_that_fails(data)

# 4. Assertion showing the bug
print(f"Expected: {expected_result}")
print(f"Actual: {result}")
assert result == expected_result, "BUG: Results don't match!"
```

---

## Minimization Checklist

When creating a minimal reproduction:

- [ ] **Remove unrelated imports**: If removing an import doesn't change the error, delete it
- [ ] **Use toy data**: Replace real data with smallest example that triggers bug (3 rows vs 1 million)
- [ ] **Remove unrelated operations**: If a line doesn't affect the bug, delete it
- [ ] **Inline functions**: Replace function calls with their minimal equivalent
- [ ] **Remove error handling**: Unless error handling is the bug, remove try/except
- [ ] **Hard-code variables**: Replace config files and variables with literal values
- [ ] **Single file**: Combine everything into one file if possible
- [ ] **Remove dependencies**: Use built-ins instead of external libraries when possible

**Goal**: 5-20 lines that trigger the exact same bug as the 1,000-line original.

---

## Benefits of Minimal Reproductions

1. **Faster debugging**: 10 lines vs 1,000 lines
2. **Easier to share**: Paste into Stack Overflow or GitHub issue
3. **Forces understanding**: Can't minimize if you don't understand the bug
4. **Often reveals solution**: Act of simplifying often makes root cause obvious
5. **Better tests**: Minimal reproduction becomes a regression test
