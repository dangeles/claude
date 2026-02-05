# Notebook Templates

Templates and guidelines for generating Jupyter notebooks with pseudocode.

## Detail Levels

### Simple (High-Level Intent)

Use for exploratory analyses where implementation details are flexible.

```python
# Cell 1: Setup
# Import necessary libraries for single-cell analysis

# Cell 2: Load Data
# Load the AnnData object containing single-cell RNA-seq data

# Cell 3: Quality Control
# Filter cells and genes based on standard QC metrics

# Cell 4: Normalization
# Normalize counts and apply log transformation

# Cell 5: Visualization
# Create UMAP embedding and visualize by condition
```

**Characteristics:**
- One-line comments per logical step
- No specific function calls
- No parameters
- Suitable for: brainstorming, early planning

### Standard (API-Level)

Use for defined analyses where methods are chosen but parameters need tuning.

```python
# Cell 1: Setup and Imports
# TODO: Import the following libraries
# - scanpy (sc): Single-cell analysis
# - pandas (pd): Data manipulation
# - numpy (np): Numerical operations
# - matplotlib.pyplot (plt): Visualization

# Cell 2: Load Data
# TODO: Load AnnData from file
# Function: sc.read_h5ad(path)
# Expected input: "{data_dir}/raw_counts.h5ad"
# Expected shape: (n_cells, n_genes) approximately (50000, 20000)

# Cell 3: Quality Control
# TODO: Apply QC filters
# Step 1: Calculate QC metrics
#   sc.pp.calculate_qc_metrics(adata, inplace=True)
#
# Step 2: Filter cells
#   sc.pp.filter_cells(adata, min_genes=200)
#   sc.pp.filter_cells(adata, max_genes=5000)
#
# Step 3: Filter genes
#   sc.pp.filter_genes(adata, min_cells=3)

# Cell 4: Normalization
# TODO: Normalize and transform
# Step 1: Total count normalization
#   sc.pp.normalize_total(adata, target_sum=1e4)
#
# Step 2: Log transformation
#   sc.pp.log1p(adata)

# Cell 5: Dimensionality Reduction
# TODO: Compute PCA and UMAP
# Step 1: Identify highly variable genes
#   sc.pp.highly_variable_genes(adata, n_top_genes=2000)
#
# Step 2: PCA
#   sc.tl.pca(adata, n_comps=50)
#
# Step 3: UMAP
#   sc.pp.neighbors(adata, n_neighbors=15)
#   sc.tl.umap(adata)
```

**Characteristics:**
- Multiple lines per step with numbered sub-steps
- Specific function names and key parameters
- Expected inputs/outputs documented
- Suitable for: implementation-ready plans

### Complex (Implementation Skeleton)

Use for critical analyses requiring error handling and validation.

```python
# Cell 1: Setup and Configuration
# TODO: Import libraries and set configuration
#
# Required imports:
#   import scanpy as sc
#   import pandas as pd
#   import numpy as np
#   from scipy import stats
#   import warnings
#
# Configuration:
#   sc.settings.verbosity = 3
#   sc.settings.figdir = "{output_dir}/figures/"
#   warnings.filterwarnings('ignore', category=FutureWarning)
#
# Parameters (adjust as needed):
#   MIN_GENES_PER_CELL = 200
#   MAX_GENES_PER_CELL = 5000
#   MIN_CELLS_PER_GENE = 3
#   N_TOP_GENES = 2000
#   N_PCS = 50
#   N_NEIGHBORS = 15

# Cell 2: Load and Validate Data
# TODO: Load data with validation
#
# def load_and_validate(path: str) -> sc.AnnData:
#     """Load AnnData with validation checks."""
#
#     # Check file exists
#     if not os.path.exists(path):
#         raise FileNotFoundError(f"Data file not found: {path}")
#
#     # Load data
#     adata = sc.read_h5ad(path)
#
#     # Validation checks
#     assert adata.n_obs > 0, "No cells in dataset"
#     assert adata.n_vars > 0, "No genes in dataset"
#
#     # Log data summary
#     print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
#     print(f"Conditions: {adata.obs['condition'].unique()}")
#
#     return adata
#
# adata = load_and_validate("{data_dir}/raw_counts.h5ad")

# Cell 3: Quality Control with Diagnostics
# TODO: QC with diagnostic plots
#
# # Calculate QC metrics
# sc.pp.calculate_qc_metrics(adata, inplace=True)
#
# # Diagnostic plots BEFORE filtering
# fig, axes = plt.subplots(1, 3, figsize=(12, 4))
#
# # Plot 1: Genes per cell distribution
# axes[0].hist(adata.obs['n_genes_by_counts'], bins=50)
# axes[0].axvline(MIN_GENES_PER_CELL, color='r', linestyle='--')
# axes[0].axvline(MAX_GENES_PER_CELL, color='r', linestyle='--')
# axes[0].set_title('Genes per cell')
#
# # Plot 2: Counts per cell distribution
# axes[1].hist(adata.obs['total_counts'], bins=50)
# axes[1].set_title('Counts per cell')
#
# # Plot 3: Cells per gene distribution
# axes[2].hist(adata.var['n_cells_by_counts'], bins=50)
# axes[2].axvline(MIN_CELLS_PER_GENE, color='r', linestyle='--')
# axes[2].set_title('Cells per gene')
#
# plt.tight_layout()
# plt.savefig(f"{sc.settings.figdir}/qc_before_filtering.png")
#
# # Record cell counts before filtering
# n_cells_before = adata.n_obs
# n_genes_before = adata.n_vars
#
# # Apply filters
# sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
# sc.pp.filter_cells(adata, max_genes=MAX_GENES_PER_CELL)
# sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
#
# # Report filtering results
# print(f"Cells: {n_cells_before} -> {adata.n_obs} ({n_cells_before - adata.n_obs} removed)")
# print(f"Genes: {n_genes_before} -> {adata.n_vars} ({n_genes_before - adata.n_vars} removed)")
#
# # Validate sufficient data remains
# assert adata.n_obs >= 100, f"Too few cells remaining: {adata.n_obs}"
# assert adata.n_vars >= 1000, f"Too few genes remaining: {adata.n_vars}"

# Cell 4: Differential Expression with Multiple Testing Correction
# TODO: DE analysis with proper statistics
#
# def run_de_analysis(adata, groupby: str, groups: list) -> pd.DataFrame:
#     """Run differential expression with BH correction."""
#
#     # Run rank_genes_groups
#     sc.tl.rank_genes_groups(
#         adata,
#         groupby=groupby,
#         groups=groups,
#         method='wilcoxon',  # Non-parametric, doesn't assume normality
#         pts=True  # Include percent expressed
#     )
#
#     # Extract results
#     result = sc.get.rank_genes_groups_df(adata, group=groups[0])
#
#     # Apply Benjamini-Hochberg correction
#     from statsmodels.stats.multitest import multipletests
#     _, result['pvals_adj_bh'], _, _ = multipletests(
#         result['pvals'],
#         method='fdr_bh',
#         alpha=0.05
#     )
#
#     # Filter significant
#     significant = result[
#         (result['pvals_adj_bh'] < 0.05) &
#         (abs(result['logfoldchanges']) > 0.5)
#     ]
#
#     print(f"Significant genes: {len(significant)} of {len(result)}")
#
#     return significant
#
# de_results = run_de_analysis(adata, groupby='condition', groups=['treatment', 'control'])
# de_results.to_csv(f"{output_dir}/de_results.csv")
```

**Characteristics:**
- Full function definitions with docstrings
- Type hints and validation
- Error handling with assertions
- Diagnostic plots and logging
- Multiple testing correction implemented
- Suitable for: production-ready implementation

## Notebook Structure

### Standard Notebook Layout

```
1. Title and Overview (Markdown)
   - Analysis name
   - Goal (one sentence)
   - Statistical approach summary

2. Setup and Imports (Code)
   - Library imports
   - Configuration
   - Parameter definitions

3. Data Loading (Code)
   - Load data
   - Initial validation
   - Data summary

4. Exploratory / QC (Code)
   - Quality metrics
   - Diagnostic plots
   - Filtering decisions

5-N. Analysis Steps (Code/Markdown alternating)
   - Each major step has:
     - Markdown: What and why
     - Code: Implementation
     - Output: Validation/visualization

N+1. Results Summary (Markdown)
   - Key findings
   - Output files generated
   - Next steps

N+2. Session Info (Code)
   - Package versions
   - Runtime environment
```

### Notebook Metadata

All generated notebooks include metadata:

```json
{
  "metadata": {
    "kernelspec": {
      "display_name": "Python 3",
      "language": "python",
      "name": "python3"
    },
    "language_info": {
      "name": "python",
      "version": "3.9"
    },
    "scientific_analysis_architect": {
      "version": "1.0.0",
      "session_id": "session-20260204-143022-12345",
      "chapter": 1,
      "notebook_number": 2,
      "analysis_type": "differential_expression",
      "detail_level": "standard",
      "generated_at": "2026-02-04T15:15:00Z"
    }
  },
  "nbformat": 4,
  "nbformat_minor": 5
}
```

## Variable Naming

### Conventions

| Type | Convention | Example |
|------|------------|---------|
| Data objects | Descriptive snake_case | `adata_filtered`, `de_results` |
| Parameters | UPPER_SNAKE_CASE | `MIN_GENES`, `N_NEIGHBORS` |
| Functions | verb_noun snake_case | `load_data()`, `run_de_analysis()` |
| Paths | descriptive with suffix | `output_dir`, `data_path` |
| Figures | fig_descriptive | `fig_qc`, `fig_umap` |

### Consistency Rules

1. **Reuse standard names across notebooks**:
   - `adata` for main AnnData object
   - `results` for analysis output DataFrames
   - `fig, ax` for matplotlib figures

2. **Document data flow**:
   ```python
   # Input: adata (raw counts, shape: n_cells x n_genes)
   # Output: adata_normalized (log-normalized, same shape)
   ```

3. **Use meaningful intermediate names**:
   ```python
   # Good
   highly_variable_genes = sc.pp.highly_variable_genes(adata, n_top_genes=2000)

   # Avoid
   hvg = sc.pp.highly_variable_genes(adata, n_top_genes=2000)
   ```

## Cell Organization

### Markdown Cells

Use markdown to explain:
- **What**: Brief description of the step
- **Why**: Statistical/biological rationale
- **Notes**: Assumptions, caveats, references

Example:
```markdown
## Differential Expression Analysis

We use the Wilcoxon rank-sum test for differential expression because:
1. Non-parametric: doesn't assume normal distribution
2. Robust to outliers common in single-cell data
3. Widely used and well-validated for scRNA-seq

Multiple testing correction: Benjamini-Hochberg (FDR < 0.05)
```

### Code Cells

Structure each code cell:
1. **Comment header** (what this cell does)
2. **Code body** (implementation or pseudocode)
3. **Validation/output** (print statements, assertions)

Example:
```python
# DIFFERENTIAL EXPRESSION: Run DE analysis between conditions
# Statistical test: Wilcoxon rank-sum
# Multiple testing: Benjamini-Hochberg correction

sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon')
de_results = sc.get.rank_genes_groups_df(adata, group='treatment')

# Apply BH correction
from statsmodels.stats.multitest import multipletests
_, de_results['pvals_adj'], _, _ = multipletests(de_results['pvals'], method='fdr_bh')

# Report significant genes
n_sig = (de_results['pvals_adj'] < 0.05).sum()
print(f"Significant genes (FDR < 0.05): {n_sig}")
```

## Validation

All notebooks must pass nbformat validation:

```python
import nbformat

def validate_notebook(path: str) -> bool:
    """Validate notebook structure."""

    with open(path) as f:
        nb = nbformat.read(f, as_version=4)

    # Structural validation
    nbformat.validate(nb)

    # Content validation
    assert len(nb.cells) > 0, "Notebook has no cells"
    assert nb.cells[0].cell_type == "markdown", "First cell should be markdown (title)"

    # Check for required metadata
    assert "kernelspec" in nb.metadata, "Missing kernelspec"
    assert "scientific_analysis_architect" in nb.metadata, "Missing provenance metadata"

    return True
```

## Template Selection Logic

```python
def select_detail_level(analysis: dict) -> str:
    """Select appropriate detail level for analysis."""

    # Complex: Critical analyses with statistical rigor requirements
    if analysis.get("statistical_complexity") == "high":
        return "complex"
    if analysis.get("critical_for_conclusions"):
        return "complex"
    if analysis.get("multiple_testing_required"):
        return "complex"

    # Simple: Exploratory or flexible implementation
    if analysis.get("exploratory"):
        return "simple"
    if analysis.get("implementation_flexible"):
        return "simple"

    # Default: Standard
    return "standard"
```

## Output File Structure

```
{output_dir}/
├── chapter1_data-atlas/
│   ├── notebook1_1_quality-control.ipynb
│   ├── notebook1_2_normalization.ipynb
│   └── notebook1_3_clustering.ipynb
├── chapter2_hypothesis-testing/
│   ├── notebook2_1_differential-expression.ipynb
│   └── notebook2_2_pathway-enrichment.ipynb
└── chapter3_mechanism-exploration/
    ├── notebook3_1_trajectory-analysis.ipynb
    └── notebook3_2_gene-regulatory-networks.ipynb
```

### Naming Convention

`notebook{chapter}_{number}_{slug}.ipynb`

- `{chapter}`: Chapter number (1-7)
- `{number}`: Notebook number within chapter (1-N)
- `{slug}`: Kebab-case analysis name
