# Reproducibility Templates

Verbatim code and markdown blocks for the Reproducibility Standards in `SKILL.md`. Copy these into notebooks; SKILL.md carries the rationale and the "when/why" prose.

## Environment Documentation (Mandatory)

Start every notebook with environment documentation:

```python
# %%
# Computational Environment
import sys
import numpy as np
import pandas as pd
import scanpy as sc  # or relevant packages

print("=" * 60)
print("COMPUTATIONAL ENVIRONMENT")
print("=" * 60)
print(f"Python: {sys.version}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")
print(f"Scanpy: {sc.__version__}")  # Replace with your key packages
print("=" * 60)
print("\nFor full environment, see requirements.txt")
```

Create environment files before starting analysis:

```bash
# For micromamba users (recommended for bioinformatics):
# Export micromamba packages:
micromamba env export > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# For pip users:
pip freeze > requirements.txt

# Document which file to use in notebook
```

In notebook markdown cell:

```markdown
## Computational Environment

- **Kernel**: Python 3.11 (bio-analysis-env)
- **Environment file**: `environment.yml` (recreate with `micromamba env create -f environment.yml`)
- **Key packages**: scanpy==1.10.0, numpy==1.26.3, pandas==2.2.0, scipy==1.12.0
- **Execution date**: 2026-01-29
```

## Random Seed Setting (Mandatory for Stochastic Processes)

Set seeds in setup cell:

```python
# %%
# Random seeds for reproducibility
import numpy as np
import random

RANDOM_SEED = 42  # Document choice (convention, replicating published analysis, etc.)

# Core Python/NumPy
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

# Scanpy (single-cell analysis)
import scanpy as sc
sc.settings.seed = RANDOM_SEED

# PyTorch (if using deep learning)
import torch
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

# TensorFlow (if using)
import tensorflow as tf
tf.random.set_seed(RANDOM_SEED)

print(f"Random seed set to {RANDOM_SEED} for reproducibility")
```

Document in notebook:

```markdown
## Stochastic Operations
This analysis uses:
- UMAP (random initialization, seed=42)
- Leiden clustering (random walk, seed=42)
- 1000-iteration permutation test (seed=42)

All seeds set to 42 for reproducibility.
```

## Session Info Output (Mandatory)

End every notebook with comprehensive session info:

```python
# %%
# Session Information for Reproducibility
import session_info

session_info.show(
    dependencies=True,
    html=False
)

# Alternative for single-cell workflows:
# import scanpy as sc
# sc.logging.print_versions()

# Alternative for base Python:
# import sys
# import pkg_resources
# print(f"Python: {sys.version}")
# for pkg in ['numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn']:
#     print(f"{pkg}: {pkg_resources.get_distribution(pkg).version}")
```

## File Path Best Practices

Use relative paths and variables:

```python
# %%
from pathlib import Path

# Define all paths at top of notebook
DATA_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
RESULTS_DIR = Path("results/analysis_2026-01-29")
FIGURES_DIR = RESULTS_DIR / "figures"

# Create output directories
for directory in [PROCESSED_DIR, RESULTS_DIR, FIGURES_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Use variables throughout
counts_file = DATA_DIR / "counts_matrix.h5ad"
metadata_file = DATA_DIR / "sample_metadata.csv"
output_file = PROCESSED_DIR / "normalized_counts.h5ad"
figure_file = FIGURES_DIR / "umap_clusters.pdf"

print(f"Data directory: {DATA_DIR.resolve()}")
print(f"Results directory: {RESULTS_DIR.resolve()}")
```

Never use hardcoded absolute paths:

```python
# ❌ BAD (non-reproducible):
adata = sc.read_h5ad("/Users/yourname/project/data/counts.h5ad")
plt.savefig("/Users/yourname/Desktop/figure.pdf")

# ✅ GOOD (reproducible):
adata = sc.read_h5ad(DATA_DIR / "counts.h5ad")
plt.savefig(FIGURES_DIR / "umap_clusters.pdf")
```

## Data Provenance Documentation

Document data sources in notebook:

```markdown
## Data Sources

### Input Data
- **File**: `data/raw/GSE123456_counts.h5ad`
- **Source**: GEO accession GSE123456
- **Download date**: 2026-01-15
- **Download command**: `wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123456`
- **Original publication**: Smith et al. (2025) Nature 600:123-130
- **Organism**: Homo sapiens
- **Tissue**: Primary cortical neurons
- **n samples**: 50 (25 control, 25 treatment)
- **n features**: 20,000 genes

### Reference Data
- **Genome build**: GRCh38 (hg38)
- **Gene annotations**: GENCODE v42
- **Downloaded**: 2026-01-10 from https://www.gencodegenes.org/
```

## Integration with notebook-writer Skill

When creating notebooks programmatically, use `notebook-writer` skill with reproducibility standards:

```python
from pathlib import Path

# Use notebook-writer to create template
cells = [
    {'type': 'markdown', 'content': '## Computational Environment\n...'},
    {'type': 'code', 'content': 'import sys\nprint(f"Python: {sys.version}")'},
    {'type': 'markdown', 'content': '## Data Loading\n...'},
    # ... analysis cells ...
    {'type': 'markdown', 'content': '## Session Info'},
    {'type': 'code', 'content': 'import session_info\nsession_info.show()'}
]

# Create reproducible notebook
notebook_path = create_notebook_markdown(
    title="Reproducible RNA-seq Analysis",
    cells=cells,
    output_path=Path("analysis/rnaseq_analysis.md")
)
```

## Bioinformatics-Specific Reproducibility Considerations

Organism and reference versions:

```python
# Document in code cell
ORGANISM = "Homo sapiens"
GENOME_BUILD = "GRCh38"  # or "mm39" for mouse, "dm6" for fly, etc.
ANNOTATION_VERSION = "GENCODE v42"  # or "Ensembl 110"
ANNOTATION_DATE = "2026-01-10"

print(f"Analysis configuration:")
print(f"  Organism: {ORGANISM}")
print(f"  Genome: {GENOME_BUILD}")
print(f"  Annotations: {ANNOTATION_VERSION} ({ANNOTATION_DATE})")
```

Bioinformatics tools (if used):

```markdown
## External Tools
- **STAR aligner**: v2.7.11a (for read mapping)
- **MACS2**: v2.2.9.1 (for peak calling)
- **bedtools**: v2.31.0 (for interval operations)

All tools available in micromamba environment (see environment.yml).
```

Data processing parameters:

```python
# Document all filtering/QC thresholds
QC_PARAMS = {
    'min_genes_per_cell': 200,
    'min_cells_per_gene': 3,
    'max_pct_mt': 15,  # percent mitochondrial reads
    'min_counts': 1000,
    'highly_variable_genes': 2000,
    'n_pcs': 50,  # principal components
    'umap_neighbors': 15,
    'leiden_resolution': 0.8
}

print("Quality control parameters:")
for param, value in QC_PARAMS.items():
    print(f"  {param}: {value}")
```
