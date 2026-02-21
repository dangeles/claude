# Notebook Best Practices for Reproducibility

Guidelines for creating notebooks that run reliably and can be reproduced by others.

---

## The Reproducibility Standard

**Definition**: Notebook is reproducible if "Restart & Run All" succeeds on a fresh environment

**Test**:
1. `Kernel → Restart & Clear Output`
2. `Kernel → Restart & Run All`
3. All cells execute successfully ✓
4. Outputs match expected results ✓

---

## Environment Documentation

### Essential: Document Your Environment

**Always include** (choose one):

**Option A: environment.yml (micromamba)**
```yaml
name: myproject
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11.2
  - numpy=1.24.3
  - pandas=2.0.1
  - matplotlib=3.7.1
  - jupyter
```

**Option B: requirements.txt (pip)**
```
# Python 3.11.2
numpy==1.24.3
pandas==2.0.1
matplotlib==3.7.1
jupyter
```

**Generate automatically**:
```bash
# micromamba:
# Export micromamba packages:
micromamba env export --no-builds > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# Pip:
pip freeze > requirements.txt
```

### Include Setup Instructions

**Add to first markdown cell**:

```markdown
# Project Title

## Environment Setup

```bash
# Create environment:
micromamba env create -f environment.yml
micromamba activate myproject

# Register Jupyter kernel:
python -m ipykernel install --user --name=myproject

# Launch notebook:
jupyter lab
```

## In Jupyter:
- Select kernel: `Python (myproject)` (top right corner)
- Run: Kernel → Restart & Run All

## Expected Runtime
- Full notebook: ~10 minutes
- Memory required: ~4GB
```

---

## Cell Organization

### Principle: Linear Narrative Flow

Cells should read like a story: setup → load → process → analyze → visualize → export

### Template Structure

```python
# Cell 1: Setup and Environment Check
"""
Environment Requirements
------------------------
See environment.yml for full dependencies.
"""
import sys
print(f"Python: {sys.version}")
assert 'myproject' in sys.executable, "Wrong environment!"

# Cell 2: Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

# Cell 3: Parameters (all configurable values here)
DATA_PATH = 'data/input.csv'
OUTPUT_PATH = 'results/output.csv'
RANDOM_SEED = 42
MIN_THRESHOLD = 0.5

# Cell 4: Load Data
df = pd.read_csv(DATA_PATH)
print(f"Loaded {len(df)} rows, {len(df.columns)} columns")

# Cell 5: Validate Data
assert not df.empty, "Dataset is empty"
assert 'value' in df.columns, "Missing 'value' column"
print("✓ Data validation passed")

# Cell 6-10: Processing steps
# ... (each with clear purpose)

# Cell 11: Export Results
df_results.to_csv(OUTPUT_PATH, index=False)
print(f"✓ Results saved to {OUTPUT_PATH}")
```

### Use Markdown Headers

```markdown
# Project Title

## 1. Setup
Cells 1-3: Environment check, imports, parameters

## 2. Data Loading
Cell 4: Load raw data from CSV

## 3. Preprocessing
Cells 5-7: Clean and validate data

## 4. Analysis
Cells 8-10: Statistical analysis and modeling

## 5. Visualization
Cells 11-13: Generate figures

## 6. Export
Cell 14: Save results
```

---

## Execution Order Best Practices

### Rule: Cells Must Run Sequentially

**Bad** (hidden dependencies):
```python
# Cell 5:
result = df.groupby('category').mean()  # Uses df

# Cell 10:
df = pd.read_csv('data.csv')  # Defines df
```

**Good** (explicit order):
```python
# Cell 3:
df = pd.read_csv('data.csv')  # Defines df first

# Cell 5:
result = df.groupby('category').mean()  # Uses df later
```

### Add Dependency Checks

```python
# Cell that depends on previous cells:
assert 'df' in dir(), "Run Cell 3 first to load data"
assert 'model' in dir(), "Run Cell 8 first to train model"
```

### Test After Every Major Change

```
After adding/modifying cells:
1. Kernel → Restart & Clear Output
2. Kernel → Restart & Run All
3. Verify all cells succeed
```

---

## Variable Naming

### Avoid In-Place Operations

**Problem**: Re-running cell gives different result

```python
# Bad (in-place):
df.dropna(inplace=True)  # Modifies df
# Re-running removes MORE rows (df already cleaned)

# Good (new variable):
df_clean = df.dropna()  # Returns new DataFrame
# Re-running always produces same result
```

### Use Descriptive Names

```python
# Bad:
df1 = pd.read_csv('data.csv')
df2 = df1[df1['value'] > 0]
df3 = df2.groupby('category').mean()

# Good:
df_raw = pd.read_csv('data.csv')
df_filtered = df_raw[df_raw['value'] > 0]
df_summary = df_filtered.groupby('category').mean()
```

### Track Data Transformations

```python
# Show pipeline:
print(f"Raw data: {len(df_raw)} rows")
df_filtered = df_raw[df_raw['value'] > 0]
print(f"After filtering: {len(df_filtered)} rows ({len(df_filtered)/len(df_raw)*100:.1f}%)")
df_clean = df_filtered.dropna()
print(f"After cleaning: {len(df_clean)} rows ({len(df_clean)/len(df_filtered)*100:.1f}%)")
```

---

## Randomness and Reproducibility

### Always Set Random Seeds

```python
# At top of notebook (Cell 2 or 3):
import numpy as np
import random

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

# For scikit-learn:
from sklearn.model_selection import train_test_split
train, test = train_test_split(df, random_state=RANDOM_SEED)

# For TensorFlow:
import tensorflow as tf
tf.random.set_seed(RANDOM_SEED)

# For PyTorch:
import torch
torch.manual_seed(RANDOM_SEED)
```

### Document Non-Deterministic Operations

```markdown
## Note: Non-Deterministic Results

Cell 8 uses GPU acceleration which may produce slightly different results across runs due to floating-point precision.

Expected variation: ±0.01 in accuracy metric.
```

---

## Memory Management

### Load Data Efficiently

```python
# Bad (loads entire 10GB file):
df = pd.read_csv('huge_file.csv')

# Good (chunked loading):
chunks = []
for chunk in pd.read_csv('huge_file.csv', chunksize=10000):
    # Process each chunk
    processed = chunk[chunk['value'] > 0]
    chunks.append(processed)
df = pd.concat(chunks, ignore_index=True)

# Better (only needed columns):
df = pd.read_csv('huge_file.csv', usecols=['id', 'value', 'category'])

# Best (use Dask for out-of-memory):
import dask.dataframe as dd
df = dd.read_csv('huge_file.csv')
result = df[df['value'] > 0].compute()
```

### Clean Up Large Objects

```python
# After processing large intermediate result:
large_intermediate = process_data(df)
final_result = summarize(large_intermediate)

del large_intermediate  # Free memory
import gc
gc.collect()
```

### Monitor Memory Usage

```python
# Add to cells that might use lots of memory:
import sys
print(f"DataFrame size: {sys.getsizeof(df) / 1e6:.2f} MB")
print(f"Memory usage:\n{df.memory_usage(deep=True)}")
```

---

## File Paths

### Use Relative Paths

```python
# Bad (hardcoded, won't work for others):
df = pd.read_csv('/Users/alice/myproject/data/input.csv')

# Good (relative to notebook):
from pathlib import Path
data_dir = Path('data')
df = pd.read_csv(data_dir / 'input.csv')
```

### Define Paths as Parameters

```python
# Cell 3: Parameters
from pathlib import Path

# Paths:
DATA_DIR = Path('data')
OUTPUT_DIR = Path('results')

# Input files:
INPUT_FILE = DATA_DIR / 'raw_data.csv'
METADATA_FILE = DATA_DIR / 'metadata.csv'

# Output files:
RESULTS_FILE = OUTPUT_DIR / 'results.csv'
FIGURES_DIR = OUTPUT_DIR / 'figures'

# Create output directories:
OUTPUT_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)
```

---

## Output Management

### Limit Output Size

```python
# Configure pandas display:
pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 10)

# Don't print in loops:
# Bad:
for i in range(1000):
    print(df.iloc[i])  # Huge notebook!

# Good:
for i in range(1000):
    if i % 100 == 0:
        print(f"Progress: {i}/1000")
```

### Save Figures, Don't Just Display

```python
# Good practice:
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Save high-resolution:
fig.savefig('figures/plot.png', dpi=300, bbox_inches='tight')
fig.savefig('figures/plot.pdf', bbox_inches='tight')  # Vector format

plt.show()
```

### Clear Outputs Before Committing

```bash
# Manual:
# Edit → Clear All Outputs

# Automated with nbstripout:
micromamba install nbstripout
nbstripout notebook.ipynb

# Add to git hooks:
nbstripout --install
```

---

## Error Handling

### Add Validation Checks

```python
# Validate input data:
assert not df.empty, "Dataset is empty"
assert 'required_column' in df.columns, "Missing required column"
assert df['value'].notna().all(), "Found NaN values in 'value' column"

# Validate processing:
assert len(df_processed) > 0, "No data remaining after processing"
assert df_processed['value'].min() >= 0, "Found negative values"
```

### Informative Error Messages

```python
# Bad:
assert len(df) > 100

# Good:
assert len(df) > 100, f"Insufficient data: {len(df)} rows (need >100)"
```

### Try-Except for External Resources

```python
# Robust file loading:
from pathlib import Path

try:
    df = pd.read_csv(DATA_PATH)
except FileNotFoundError:
    print(f"Error: File not found: {DATA_PATH}")
    print(f"Current directory: {Path.cwd()}")
    print(f"Expected location: {Path(DATA_PATH).absolute()}")
    raise
```

---

## Documentation

### Cell Comments

```python
# Good cell comments explain WHY, not WHAT:

# Bad:
# Load data from CSV
df = pd.read_csv('data.csv')

# Good:
# Load patient data
# Note: This dataset uses 1-based indexing for patient IDs
df = pd.read_csv('data.csv')
```

### Markdown Explanations

```markdown
## Data Preprocessing

We apply the following filters:
1. Remove rows with missing values (typically <5% of data)
2. Filter for quality score >0.8 (excludes low-quality measurements)
3. Log-transform expression values to normalize distribution

**Expected result**: ~90-95% of original rows retained
```

### Session Info

```python
# Final cell: Document exact environment
import sys
import numpy as np
import pandas as pd
import matplotlib

print("## Session Info")
print(f"Python: {sys.version}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")
print(f"Matplotlib: {matplotlib.__version__}")

# Alternatively:
!pip list
```

---

## Testing Reproducibility

### Pre-Flight Checklist

Before sharing notebook:
- [ ] Environment documented (environment.yml or requirements.txt)
- [ ] Setup instructions in first markdown cell
- [ ] All file paths are relative
- [ ] Random seeds set
- [ ] "Restart & Run All" succeeds
- [ ] Outputs cleared (Edit → Clear All Outputs)
- [ ] Session info included

### Fresh Environment Test

```bash
# Create test environment:
micromamba env create -f environment.yml -n test_env
micromamba activate test_env
python -m ipykernel install --user --name=test_env

# Run notebook:
jupyter lab
# Select test_env kernel, run Restart & Run All

# Should succeed without modification

# Clean up:
micromamba deactivate
micromamba env remove -n test_env
jupyter kernelspec uninstall test_env
```

---

## Common Pitfalls to Avoid

### 1. Global State Pollution

```python
# Bad (modifies global state):
import matplotlib.pyplot as plt
plt.figure()
plt.plot(x1, y1)
# ... many cells later ...
plt.plot(x2, y2)  # Accidentally adds to same figure!

# Good (explicit figure management):
fig, ax = plt.subplots()
ax.plot(x1, y1)
plt.show()

fig, ax = plt.subplots()  # New figure
ax.plot(x2, y2)
plt.show()
```

### 2. Import Inside Loops

```python
# Bad (slow):
for file in files:
    import pandas as pd  # Re-imports every iteration!
    df = pd.read_csv(file)

# Good (import once):
import pandas as pd
for file in files:
    df = pd.read_csv(file)
```

### 3. Assuming Execution Order

```python
# Bad (assumes Cell 5 ran before Cell 3):
# Cell 3:
result = process(intermediate_data)

# Cell 5:
intermediate_data = load_and_prepare()

# Good (dependencies clear):
# Cell 3:
intermediate_data = load_and_prepare()

# Cell 5:
result = process(intermediate_data)
```

### 4. Hardcoded Timestamps

```python
# Bad (timestamp in filename):
df.to_csv('results_2024-01-15.csv')
# Next run creates different file!

# Good (overwrite or parameterize):
df.to_csv('results.csv')
# Or:
from datetime import datetime
timestamp = datetime.now().strftime('%Y%m%d')
df.to_csv(f'results_{timestamp}.csv')
# And document this behavior
```

---

## Version Control

### What to Commit

✅ **Commit**:
- Notebook file (.ipynb)
- environment.yml or requirements.txt
- README with setup instructions
- Small data files (<10MB)

❌ **Don't commit**:
- Large data files (>10MB)
- Outputs (clear before commit)
- .ipynb_checkpoints/
- __pycache__/

### .gitignore for Notebooks

```gitignore
# Jupyter
.ipynb_checkpoints/
*/.ipynb_checkpoints/*

# Python
__pycache__/
*.py[cod]

# Data
data/raw/*
data/processed/*
!data/raw/.gitkeep
!data/processed/.gitkeep

# Results
results/*
!results/.gitkeep

# Environment
.env
```

### Commit Messages

```bash
# Good commit messages for notebooks:
git commit -m "Add exploratory data analysis notebook"
git commit -m "Fix cell execution order in preprocessing"
git commit -m "Update environment.yml: add scanpy dependency"
git commit -m "Clear outputs before commit"
```

---

## Advanced: Automation

### Pre-commit Hook

```bash
# .git/hooks/pre-commit
#!/bin/bash

# Strip outputs before commit:
jupyter nbconvert --clear-output --inplace notebooks/*.ipynb

# Test reproducibility:
jupyter nbconvert --execute --to notebook --inplace notebooks/analysis.ipynb

exit 0
```

### CI/CD Testing

```yaml
# .github/workflows/test-notebooks.yml
name: Test Notebooks

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2

      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.11'

      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          micromamba install nbconvert pytest

      - name: Test notebooks
        run: |
          jupyter nbconvert --execute --to notebook notebooks/*.ipynb
```

---

## Quick Reference

**Reproducibility Checklist**:
- [ ] Environment file exists (environment.yml or requirements.txt)
- [ ] Setup instructions documented
- [ ] Relative file paths only
- [ ] Random seeds set
- [ ] Cells run sequentially (no out-of-order dependencies)
- [ ] "Restart & Run All" succeeds
- [ ] Outputs cleared before sharing
- [ ] Session info included

**Before Sharing**:
```bash
# 1. Clear outputs:
jupyter nbconvert --clear-output --inplace notebook.ipynb

# 2. Export environment:
# Export micromamba packages:
micromamba env export --no-builds > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# 3. Test in fresh environment:
micromamba env create -f environment.yml -n test
micromamba activate test
jupyter lab
# Run "Restart & Run All"
```

---

## Resources

- [Jupyter Best Practices](https://jupyter-notebook.readthedocs.io/en/stable/notebook.html)
- [Ten Simple Rules for Reproducible Research](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003285)
- [nbstripout Documentation](https://github.com/kynan/nbstripout)
- [papermill: Parameterize Notebooks](https://papermill.readthedocs.io/)
