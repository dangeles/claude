---
name: notebook-writer
last_updated: 2026-06-09
description: Use when creating reproducible Jupyter notebooks — structuring analysis narrative, documenting cells, capturing environment/session info, and producing git-friendly (Jupytext) notebooks. NOT for debugging broken notebooks or kernel crashes (use notebook-debugger) or implementing the underlying analysis logic (use bioinformatician).
success_criteria:
  - Notebook structured with clear sections and narrative
  - Code cells properly documented with explanatory markdown
  - Reproducibility ensured (can run end-to-end)
  - Environment requirements documented
  - Outputs properly displayed and interpreted
  - Git-friendly format (Jupytext markdown)
---

# Notebook Writer Skill

You are a specialist in creating well-structured Jupyter notebooks for scientific analyses and documentation.

## When to Use This Skill

Use this skill when:
- Creating parameter sweeps or sensitivity analyses
- Documenting calculations with reproducible code
- Generating analysis reports that combine code, results, and interpretation
- Packaging agent work (Calculator, Researcher) into shareable notebooks

## Notebook Format: Jupytext Markdown

We use Jupytext-compatible Markdown for notebooks to enable git-friendly version control.

### Cell Markers

- **Markdown cells**: Regular Markdown text (no special marker)
- **Code cells**: Start with `# %%` on its own line

### Example Structure

```markdown
---
jupyter:
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Analysis Title

Brief description of what this notebook does.

## Section 1: Data Loading

# %%
import pandas as pd
import numpy as np

# %%
data = pd.read_csv('data.csv')
data.head()

## Section 2: Analysis

Explanation of the analysis approach.

# %%
# Perform calculation
result = np.mean(data['value'])
print(f"Mean: {result:.2f}")
```

## Python Utility API

Many projects provide `src/utils/notebook_builder.py` with helper functions for programmatic notebook creation: `create_notebook_markdown`, `create_parameter_sweep_notebook`, `create_analysis_report_notebook`, and `validate_notebook`. Signatures, parameters, and worked examples are in `references/notebook-builder-api.md`.

## Workflow

The Jupytext round trip, in order:

1. **Create notebook** using the utility or manual Markdown
2. **Edit `.md` file** directly (agents write Markdown well)
3. **Convert to `.ipynb`**: `python3 -m jupytext --to ipynb notebook.md`
4. **Run in Jupyter**: `jupyter notebook notebook.ipynb`
5. **Sync changes back**: `python3 -m jupytext --sync notebook.ipynb` (bidirectional)

## Jupytext Configuration

Projects should include `.jupytext.toml` in repository root:

```toml
# Jupytext configuration
# Enables git-friendly notebook version control

# Pair markdown and ipynb files
# Use myst format which supports # %% cell markers
formats = "md:myst,ipynb"
```

This tells Jupytext to:
- Recognize `.md` files as notebooks
- Use MyST Markdown format (supports `# %%` markers)
- Auto-sync with `.ipynb` when either is modified

## Git Tracking Strategy

**Recommended `.gitignore` configuration:**
```gitignore
# Track .md notebooks (Jupytext source), ignore generated .ipynb files
*.ipynb
.ipynb_checkpoints/
```

**What's tracked:**
- `.md` notebook files (human-readable source)
- Not `.ipynb` files (generated, binary JSON)
- Not `.ipynb_checkpoints/` (Jupyter temp files)

**Rationale:** `.md` files produce readable git diffs. `.ipynb` files are JSON with embedded outputs and can be regenerated from `.md`.

## Common Operations

### Create from scratch (manual)
1. Write Markdown file with `# %%` markers
2. Add YAML frontmatter (kernel info)
3. Convert: `python3 -m jupytext --to ipynb file.md`

### Convert existing notebook to Markdown
```bash
python3 -m jupytext --to md:myst notebook.ipynb
```

### Edit existing notebook
**Option 1: Edit .md file directly** (recommended for agents)
```bash
# Edit notebook.md in text editor
# Then convert:
python3 -m jupytext --to ipynb notebook.md
```

**Option 2: Edit in Jupyter, sync back**
```bash
jupyter notebook notebook.ipynb
# Make changes in Jupyter
# Sync back to .md:
python3 -m jupytext --sync notebook.ipynb
```

### Validate structure
```bash
python3 -c "
from pathlib import Path
from src.utils.notebook_builder import validate_notebook
validate_notebook(Path('notebook.ipynb'))
print('✓ Valid')
"
```

### Convert multiple notebooks
```bash
# Convert all .md notebooks in a directory
python3 -m jupytext --to ipynb analysis/*.md

# Or sync all paired notebooks
python3 -m jupytext --sync analysis/*.ipynb
```

## Narrative and structure

A notebook is a document that happens to execute. Give it a title that states its purpose, group imports at the top, and let markdown carry the reasoning: what a calculation is for before the cell, what the numbers mean after it. Meaningful variable names, units in comments and axis labels, and figures saved to files all make the notebook readable months later. How finely to section the narrative is a judgment call — follow the shape of the analysis rather than a fixed template.

## Reproducibility Standards

Scientific notebooks should let another researcher recreate the computational environment, rerun the analysis and get identical results, and understand the data sources and processing steps.

### Environment Documentation

Include an environment documentation cell:

```python
# %%
# Environment Information
# Run: pip freeze > requirements.txt
# Or: micromamba env export > environment.yml

import sys
import numpy as np
import pandas as pd
import scanpy as sc  # Example for single-cell analysis

print(f"Python: {sys.version}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")
print(f"Scanpy: {sc.__version__}")

# Include this output in your notebook for documentation
```

**Create environment files:**

```bash
# For pip users:
pip freeze > requirements.txt

# For micromamba users:
# Export micromamba packages:
micromamba env export > environment.yml

# Export pip-installed packages separately (micromamba export does not include pip packages):
pip freeze > pip-requirements.txt

# Include these files in your repository
```

**Document kernel selection:**
```markdown
## Computational Environment

- **Kernel**: Python 3.11 (project-env)
- **Dependencies**: See `requirements.txt` for full package list
- **Critical packages**: scanpy==1.10.0, numpy==1.26.3, pandas==2.2.0
```

### Random Seed Setting

For any stochastic process, set random seeds:

```python
# %%
# Set random seeds for reproducibility
import numpy as np
import random

RANDOM_SEED = 42  # Document why this value was chosen (convention, previous analysis, etc.)

np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)

# For machine learning:
import torch
torch.manual_seed(RANDOM_SEED)

# For scanpy:
import scanpy as sc
sc.settings.seed = RANDOM_SEED

print(f"Random seed set to {RANDOM_SEED}")
```

**Stochastic processes requiring seeds:**
- UMAP, t-SNE (dimensionality reduction)
- Random forest, neural networks (machine learning)
- Monte Carlo simulations
- Random sampling or bootstrapping
- Graph algorithms with random initialization

### Session Info Output

End the notebook with a session info cell:

```python
# %%
# Session Information (for reproducibility)
import session_info

session_info.show(
    dependencies=True,
    html=False
)

# Alternative for single-cell analysis:
# import scanpy as sc
# sc.logging.print_versions()
```

This captures:
- Python version
- Operating system
- Package versions (all dependencies)
- Execution timestamp

### File Path Best Practices

Use relative paths and variables:

```python
# %%
from pathlib import Path

# Define paths at the top of the notebook
DATA_DIR = Path("data/raw")
RESULTS_DIR = Path("results/analysis_2025-01-29")

# Ensure output directories exist
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Use variables throughout
input_file = DATA_DIR / "counts.csv"
output_file = RESULTS_DIR / "normalized_counts.csv"
```

Avoid hardcoded absolute paths:
```python
# BAD:
data = pd.read_csv("/Users/yourname/project/data.csv")

# GOOD:
data = pd.read_csv(DATA_DIR / "data.csv")
```

### End-to-end execution

Before sharing or archiving, run the notebook on a fresh kernel (Restart & Run All) and confirm it completes; if it is a re-run of an existing analysis, compare results against the expected output. That execution is the check that matters — it exercises the environment file, the seeds, and the path handling at once. Document data sources (where to download, version, date) alongside the environment files.

## AI assistance inside Jupyter

JupyterLab's `%%ai` magics, chat UI, and the JetBrains AI Assistant can help with boilerplate and exploration. See `references/jupyter-ai-integration.md`.

## Project-Specific Usage

Many projects have a `docs/NOTEBOOK-WORKFLOW.md` or similar document with project-specific examples and patterns. Check your project's documentation for:
- Domain-specific notebook templates
- Agent integration patterns (which skills create notebooks)
- Directory structure conventions (where to save notebooks)
- Project-specific best practices

## Troubleshooting

### Issue: Jupytext can't find format

**Error:** `Format 'percent' is not associated to extension '.md'`

**Fix:** Use `md:myst` format in `.jupytext.toml` (not `md:percent`). MyST Markdown supports `# %%` markers.

```toml
formats = "md:myst,ipynb"  # Correct
```

### Issue: Sync not working

**Symptom:** Changes to `.ipynb` don't appear in `.md`

**Solution:**
1. Check `.jupytext.toml` exists and has correct format
2. Run sync explicitly: `python3 -m jupytext --sync notebook.ipynb`
3. Check both files exist (create `.md` first if needed)

### Issue: Validation fails

**Error:** `nbformat.ValidationError`

**Causes:**
- Missing cell IDs (required in nbformat v4.5+)
- Invalid JSON structure
- Missing required fields

**Solution:** Use `notebook_builder.py` utility functions which handle validation automatically.

### Issue: Git shows .ipynb files

**Symptom:** `.ipynb` files appearing in `git status`

**Fix:** Ensure `.gitignore` contains `*.ipynb`. Check with:
```bash
git check-ignore -v notebook.ipynb
```

### Issue: Cells not recognized

Code cells need the `# %%` marker, and the YAML frontmatter (kernelspec) must be exact — see the example structure above. For `.ipynb` files created directly, `notebook_builder.validate_notebook()` runs the nbformat check.

## Dependencies

**Required packages:**
```bash
micromamba install jupytext nbformat
```

**Check installation:**
```bash
pip3 list | grep -E "(jupytext|nbformat)"
```

**Tested versions:**
- jupytext: 1.19+
- nbformat: 5.10+
- Python: 3.9+

## Integration with Other Skills

- **notebook-debugger**: verify end-to-end execution, fix kernel/environment failures
- **bioinformatician**: apply reproducibility standards to computational biology analyses
- **copilot**: review notebook code for correctness and reproducibility compliance
- **Quantitative analysis skills**: package calculations as reproducible notebooks with parameter sweeps
- **Research skills**: document literature-derived parameters with citations in data notebooks
- **Planning skills**: generate protocol notebooks with expected results and analysis templates

---

Remember: Notebooks are for interactive exploration and reproducible documentation. For production code, use Python modules in `src/`.

**For project-specific examples and patterns**, see your project's documentation (often `docs/NOTEBOOK-WORKFLOW.md` or similar).
