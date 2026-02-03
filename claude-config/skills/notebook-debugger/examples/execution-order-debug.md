# Execution Order Debug Example: Cell Dependency Issue

## Problem Report

**User complaint**: "Notebook works when I run cells manually, but fails on 'Restart & Run All' with `NameError: name 'df' is not defined`"

**Notebook**: `data_pipeline.ipynb`
**Cell that fails**: Cell 5 (when run with Restart & Run All)
**Environment**: Jupyter Lab 3.6, Python 3.11

---

## Phase 1: Diagnose

### Initial Investigation

**Test 1: Reproduce error with "Restart & Run All"**
```
Action: Kernel → Restart & Run All
Result: Cell 5 fails with NameError: name 'df' is not defined
```

**Test 2: Check execution numbers**

Looking at cell execution numbers (the `[n]` next to each cell):

```python
[12] # Cell 1: Import libraries
[13] # Cell 2: Set parameters
[8]  # Cell 3: Load data
[14] # Cell 4: Process data
[15] # Cell 5: Analyze data
```

**Observation**: Cell 3 has execution number `[8]`, earlier than cells 1-2. This means cells were run out of order during development.

**Test 3: Examine Cell 5 code**
```python
# Cell 5: Analyze data
result = df.groupby('category').agg({
    'value': ['mean', 'std', 'count']
})
print(result)

# Error on Restart & Run All:
# NameError: name 'df' is not defined
```

**Test 4: Check where 'df' is defined**

Looking through all cells:

```python
# Cell 1: Import libraries
import pandas as pd
import numpy as np

# Cell 2: Set parameters
DATA_PATH = 'data/raw_data.csv'
MIN_VALUE = 0

# Cell 3: Load data (NOT RUN - marked with comment)
# df = pd.read_csv(DATA_PATH)  # TEMP: Using test data instead

# Cell 4: Process data
# Uses df_filtered (defined elsewhere?)

# Cell 7: Create test data (FOUND IT!)
df = pd.DataFrame({
    'category': ['A', 'B', 'C'] * 10,
    'value': np.random.randn(30)
})
```

**Root Cause Identified**: Cell 7 creates `df`, but Cell 5 tries to use `df`. During development, user ran Cell 7 before Cell 5. On "Restart & Run All", cells run in displayed order (1→2→3→4→5), so Cell 7 never runs before Cell 5.

---

## Phase 2: Isolate

### Hypothesis: Hidden State from Out-of-Order Execution

**Test: Map dependencies**

Creating dependency graph:

```
Cell 1: import pandas, numpy → (no dependencies)
Cell 2: DATA_PATH, MIN_VALUE → (no dependencies)
Cell 3: [COMMENTED OUT] would create df
Cell 4: uses df_filtered → depends on ???
Cell 5: uses df → depends on Cell 7 (out of order!)
Cell 7: creates df → (no dependencies)
```

**Problem 1**: Cell 5 depends on Cell 7, but Cell 7 comes after Cell 5
**Problem 2**: Cell 4 uses `df_filtered` but it's never defined
**Problem 3**: Cell 3 is commented out but should load real data

### Verify Hidden Variables

**Add diagnostic cell**:
```python
# New Cell 6: Check namespace
print("Variables in namespace:")
print([var for var in dir() if not var.startswith('_')])

# Output (after running cells out of order):
# ['DATA_PATH', 'MIN_VALUE', 'df', 'df_filtered', 'np', 'pd', 'result']

# Output (after Restart & Run All):
# ['DATA_PATH', 'MIN_VALUE', 'np', 'pd']
# (df and df_filtered missing!)
```

**Conclusion**: `df` and `df_filtered` exist from previous runs, but aren't created by cells 1-5 in order

---

## Phase 3: Fix

### Solution 1: Reorganize Cells in Logical Order

**Fixed cell order**:

```python
# Cell 1: Import libraries
import pandas as pd
import numpy as np

# Cell 2: Set parameters
DATA_PATH = 'data/raw_data.csv'
MIN_VALUE = 0

# Cell 3: Load data (MOVED FROM CELL 7, UNCOMMENTED)
# Option A: Load real data
df = pd.read_csv(DATA_PATH)

# Option B: Create test data (if file doesn't exist)
# Uncomment for testing:
# df = pd.DataFrame({
#     'category': ['A', 'B', 'C'] * 10,
#     'value': np.random.randn(30)
# })

print(f"Loaded data: {df.shape}")

# Cell 4: Filter data (FIXED - renamed from df_filtered)
df_filtered = df[df['value'] >= MIN_VALUE]
print(f"After filtering: {df_filtered.shape}")

# Cell 5: Process data (UPDATED - uses df_filtered correctly)
df_processed = df_filtered.copy()
df_processed['value_squared'] = df_processed['value'] ** 2

# Cell 6: Analyze data (RENUMBERED - was Cell 5)
result = df_processed.groupby('category').agg({
    'value': ['mean', 'std', 'count']
})
print(result)
```

**Deleted**: Old Cell 7 (test data - merged into Cell 3 as comment)

### Solution 2: Add Dependency Checks

**Defensive programming** to catch execution order issues:

```python
# Cell 5: Process data
assert 'df_filtered' in dir(), \
    "Error: df_filtered not defined. Did you run Cell 4?"

df_processed = df_filtered.copy()
df_processed['value_squared'] = df_processed['value'] ** 2
```

**Benefits**:
- Fails fast with clear error message
- Tells user exactly what's wrong
- Better than cryptic `NameError`

### Solution 3: Add Section Headers

**Markdown cells** to organize workflow:

```markdown
# Data Pipeline

## 1. Setup
Cells 1-2: Import libraries and set parameters

## 2. Load Data
Cell 3: Load raw data from CSV

## 3. Preprocessing
Cells 4-5: Filter and process data

## 4. Analysis
Cell 6: Aggregate by category
```

---

## Phase 4: Verify

### Verification Test 1: Restart & Run All

```
Action: Kernel → Restart & Clear Output
Action: Kernel → Restart & Run All
Result: All cells execute successfully ✓
```

**Check execution numbers**:
```
[1] # Cell 1: Import libraries
[2] # Cell 2: Set parameters
[3] # Cell 3: Load data
[4] # Cell 4: Filter data
[5] # Cell 5: Process data
[6] # Cell 6: Analyze data
```

All sequential ✓

### Verification Test 2: Clean Environment Test

**Create new virtualenv and test**:
```bash
# Terminal:
conda create -n test_pipeline python=3.11
conda activate test_pipeline
pip install pandas numpy jupyter
jupyter lab

# Open notebook, run Restart & Run All
# Result: Works ✓
```

### Verification Test 3: No Hidden State

**Check namespace after each cell**:
```python
# Add to end of each cell during verification:
print(f"Variables after Cell X: {[v for v in dir() if not v.startswith('_')]}")
```

**Expected output** (after Cell 6):
```
Variables after Cell 1: ['np', 'pd']
Variables after Cell 2: ['DATA_PATH', 'MIN_VALUE', 'np', 'pd']
Variables after Cell 3: ['DATA_PATH', 'MIN_VALUE', 'df', 'np', 'pd']
Variables after Cell 4: ['DATA_PATH', 'MIN_VALUE', 'df', 'df_filtered', 'np', 'pd']
Variables after Cell 5: ['DATA_PATH', 'MIN_VALUE', 'df', 'df_filtered', 'df_processed', 'np', 'pd']
Variables after Cell 6: ['DATA_PATH', 'MIN_VALUE', 'df', 'df_filtered', 'df_processed', 'np', 'pd', 'result']
```

All variables introduced in expected order ✓

---

## Phase 5: Document

### Added Markdown Cell (Workflow Overview)

```markdown
## Notebook Structure

**This notebook must be run in order (Restart & Run All).**

### Workflow:
1. **Setup** (Cells 1-2): Import libraries, set parameters
2. **Load** (Cell 3): Read data from CSV or create test data
3. **Preprocess** (Cells 4-5): Filter and transform data
4. **Analyze** (Cell 6): Aggregate and summarize results

### Running the notebook:
- **Recommended**: Kernel → Restart & Run All
- **During development**: Run cells in displayed order
- **Testing changes**: Restart kernel before re-running

### Troubleshooting:
- `NameError: name 'X' is not defined`: Restart kernel and run all cells in order
- File not found: Update DATA_PATH in Cell 2
```

### Added Cell Comments

```python
# Cell 3: Load data
# NOTE: This cell must run BEFORE cells 4-6
# Defines: df

df = pd.read_csv(DATA_PATH)
print(f"✓ Loaded {len(df)} rows")
```

```python
# Cell 4: Filter data
# NOTE: Requires df from Cell 3
# Defines: df_filtered

assert 'df' in dir(), "Run Cell 3 first to load data"
df_filtered = df[df['value'] >= MIN_VALUE]
print(f"✓ Filtered to {len(df_filtered)} rows")
```

### Added Jupyter Notebook Metadata

**Enable execution order enforcement** (in notebook JSON):

```json
{
  "metadata": {
    "kernelspec": {
      "display_name": "Python 3",
      "language": "python",
      "name": "python3"
    },
    "execution": {
      "allow_errors": false,
      "timeout": 60
    }
  }
}
```

---

## Summary

| Phase | Action | Outcome |
|-------|--------|---------|
| Diagnose | Check execution numbers, find where df defined | Cell 5 uses df from Cell 7 (out of order) |
| Isolate | Map dependencies between cells | Cell 7 must run before Cell 5 |
| Fix | Reorganize cells, add assertions, add sections | Cells now in logical order |
| Verify | Restart & Run All, clean environment test | All cells execute sequentially ✓ |
| Document | Add workflow overview, cell comments, assertions | Clear execution requirements ✓ |

---

## Key Lessons

1. **Execution numbers reveal history**: `[8]` appearing before `[12]` means cells run out of order
2. **"Restart & Run All" is the test**: Only reliable way to verify reproducibility
3. **Hidden state is dangerous**: Variables in memory but not in visible cells
4. **Assertions catch errors early**: Better than `NameError` late in pipeline
5. **Document execution order**: Make it clear cells must run in sequence

---

## Prevention for Future Notebooks

**Development discipline checklist**:
- [ ] Run cells in displayed order during development
- [ ] Test with "Restart & Run All" after every major change
- [ ] Add assertions for critical variables: `assert 'df' in dir()`
- [ ] Use markdown headers to organize sections
- [ ] Clear all outputs before committing: Edit → Clear All Outputs
- [ ] Document any cells that can be run out of order (e.g., visualization cells)

**Code pattern for dependency checks**:
```python
# Template for cells with dependencies:
def check_dependency(var_name, source_cell):
    """Check if required variable exists"""
    if var_name not in dir():
        raise RuntimeError(
            f"Variable '{var_name}' not defined. "
            f"Run {source_cell} first."
        )

# Usage:
check_dependency('df', 'Cell 3')
check_dependency('df_filtered', 'Cell 4')
```

**Automation**:
```python
# Add to first cell to enable strict mode:
import sys

class StrictNotebook:
    """Detect out-of-order execution"""
    def __init__(self):
        self.expected_cell = 1

    def mark_cell(self, cell_number):
        if cell_number != self.expected_cell:
            print(f"⚠️  Warning: Running Cell {cell_number}, expected Cell {self.expected_cell}")
            print(f"   Consider restarting kernel and running in order")
        self.expected_cell = cell_number + 1

strict = StrictNotebook()

# In each cell:
strict.mark_cell(1)  # Cell 1
strict.mark_cell(2)  # Cell 2
# etc.
```

---

## Related Issues

**Issue: "Works for me but not colleague"**
- Likely cause: You ran cells out of order, they didn't
- Fix: Restart & Run All, verify it works, then share

**Issue: "Notebook was working yesterday, now broken"**
- Check: Did you restart kernel?
- Check: Are execution numbers sequential?
- Fix: Restart & Run All

**Issue: In-place operations breaking re-runs**
```python
# Bad: Modifies df in place
df.dropna(inplace=True)
# Re-running this cell removes MORE rows each time

# Good: Returns new DataFrame
df_clean = df.dropna()
# Re-running this always produces same result
```

**Issue: Random results on each run**
```python
# Problem:
df['random'] = np.random.randn(len(df))
# Different values each run

# Solution: Set random seed
np.random.seed(42)
df['random'] = np.random.randn(len(df))
# Reproducible
```
