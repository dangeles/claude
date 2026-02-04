# Structure Template: Research Projects

## Overview

This reference defines directory structure templates for research projects (scientific research, academic work, data analysis). Based on the Cookiecutter Data Science template and academic best practices.

## Standard Research Structure

### Full Template

```
project/
├── data/
│   ├── raw/                  # Original, immutable data
│   │   └── README.md         # Data provenance
│   ├── interim/              # Intermediate data
│   ├── processed/            # Final, canonical datasets
│   └── external/             # Third-party data sources
├── notebooks/                # Jupyter notebooks
│   ├── 01-data-exploration.ipynb
│   ├── 02-data-processing.ipynb
│   └── 03-analysis.ipynb
├── docs/
│   ├── literature/           # Paper notes, reviews
│   │   ├── review-{topic}.md
│   │   └── {author}-{year}.md
│   ├── analysis/             # Analysis documents
│   │   └── analysis-{question}.md
│   ├── protocols/            # Lab protocols
│   │   └── protocol-{procedure}.md
│   └── meetings/             # Meeting notes
│       └── {YYYY-MM-DD}-meeting.md
├── models/                   # Trained models
│   └── {model-name}/
│       ├── model.pkl
│       └── config.yaml
├── references/               # Data dictionaries, manuals
│   ├── data-dictionary.md
│   └── variable-definitions.md
├── reports/
│   ├── figures/              # Generated graphics
│   │   └── fig-{NN}-{desc}.png
│   └── drafts/               # Report drafts
│       └── draft-{YYYY-MM-DD}.md
├── src/                      # Source code (if any)
│   └── {package_name}/
│       ├── __init__.py
│       ├── data.py           # Data loading
│       ├── features.py       # Feature engineering
│       ├── models.py         # Model code
│       └── visualization.py  # Plotting
├── tests/                    # Tests for src/
│   └── test_data.py
├── environment.yml           # Conda environment
├── requirements.txt          # Pip requirements
├── Makefile                  # Automation
├── README.md
├── CLAUDE.md
└── .gitignore
```

### Minimal Structure (Small Projects)

```
project/
├── data/
│   ├── raw/
│   └── processed/
├── notebooks/
│   └── analysis.ipynb
├── docs/
│   └── notes.md
├── reports/
│   └── figures/
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Computational Research

### Machine Learning Project

```
project/
├── data/
│   ├── raw/
│   ├── processed/
│   └── features/             # Feature store
├── notebooks/
│   ├── 01-eda.ipynb
│   ├── 02-features.ipynb
│   ├── 03-training.ipynb
│   └── 04-evaluation.ipynb
├── src/
│   ├── __init__.py
│   ├── data/
│   │   ├── __init__.py
│   │   ├── load.py
│   │   └── preprocess.py
│   ├── features/
│   │   ├── __init__.py
│   │   └── build.py
│   ├── models/
│   │   ├── __init__.py
│   │   ├── train.py
│   │   └── predict.py
│   └── visualization/
│       └── __init__.py
├── models/
│   └── {experiment-name}/
├── experiments/              # Experiment tracking
│   └── {experiment-id}/
├── configs/
│   └── config.yaml
├── tests/
├── MLproject                 # MLflow
├── README.md
├── CLAUDE.md
└── .gitignore
```

### Bioinformatics Project

```
project/
├── data/
│   ├── raw/
│   │   ├── sequences/        # FASTA, FASTQ
│   │   └── annotations/      # GFF, GTF
│   ├── processed/
│   │   ├── alignments/       # BAM, SAM
│   │   └── counts/           # Feature counts
│   └── external/
│       └── databases/        # Reference databases
├── workflows/                # Snakemake/Nextflow
│   ├── Snakefile
│   └── rules/
├── notebooks/
├── scripts/                  # Analysis scripts
│   ├── preprocess.py
│   └── analysis.R
├── envs/                     # Conda environments
│   └── environment.yml
├── results/
│   ├── figures/
│   └── tables/
├── docs/
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Lab-Based Research

### Experimental Research

```
project/
├── data/
│   ├── raw/
│   │   ├── {YYYYMMDD}-experiment-{N}/
│   │   └── README.md
│   └── processed/
├── notebooks/
│   └── {YYYY-MM-DD}-{description}.ipynb
├── docs/
│   ├── protocols/
│   │   └── protocol-{procedure}.md
│   ├── lab-notebook/
│   │   └── {YYYY-MM-DD}-entry.md
│   └── literature/
├── reports/
│   ├── figures/
│   └── manuscript/
│       ├── main.tex
│       └── figures/
├── analysis/                 # Analysis scripts
│   └── {experiment}/
├── README.md
├── CLAUDE.md
└── .gitignore
```

### Wet Lab + Computational

```
project/
├── data/
│   ├── raw/
│   │   ├── wetlab/           # Instrument outputs
│   │   └── computational/    # Simulations
│   └── processed/
├── protocols/                # Lab protocols
│   └── protocol-{name}.md
├── experiments/              # Organized by experiment
│   └── {experiment-name}/
│       ├── data/
│       ├── analysis/
│       └── notes.md
├── notebooks/
├── src/                      # Computational tools
├── docs/
├── reports/
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Directory Purposes

### Data Directories

| Directory | Purpose | Rules |
|-----------|---------|-------|
| `data/raw/` | Original data | NEVER modify |
| `data/interim/` | Intermediate processing | Reproducible |
| `data/processed/` | Final datasets | Used for analysis |
| `data/external/` | Third-party data | Document source |

### Documentation Directories

| Directory | Purpose | Contents |
|-----------|---------|----------|
| `docs/literature/` | Paper notes | Reviews, summaries |
| `docs/analysis/` | Analysis docs | Methods, questions |
| `docs/protocols/` | Procedures | Lab protocols |
| `docs/meetings/` | Meeting notes | Dated entries |

### Output Directories

| Directory | Purpose | Contents |
|-----------|---------|----------|
| `reports/figures/` | Generated graphics | Numbered figures |
| `reports/drafts/` | Report drafts | Dated versions |
| `models/` | Trained models | Versioned models |
| `results/` | Analysis results | Tables, summaries |

## Notebook Organization

### Numbering Convention

```
notebooks/
├── 00-{initials}-setup.ipynb           # Environment setup
├── 01-{initials}-data-exploration.ipynb
├── 02-{initials}-data-cleaning.ipynb
├── 03-{initials}-feature-engineering.ipynb
├── 10-{initials}-model-baseline.ipynb
├── 11-{initials}-model-tuning.ipynb
├── 20-{initials}-evaluation.ipynb
└── 30-{initials}-visualization.ipynb
```

**Number Ranges**:
- 00-09: Setup, data loading
- 10-19: Data processing, features
- 20-29: Modeling
- 30-39: Evaluation, visualization
- 40-49: Reporting

## Data Management

### Raw Data README

Every `data/raw/` should have a README:

```markdown
# Raw Data

## Sources

### {dataset-name}
- **Source**: [URL or description]
- **Downloaded**: YYYY-MM-DD
- **Format**: CSV/Parquet/etc.
- **Size**: X MB
- **Description**: Brief description

## Files

| File | Description | Source | Date |
|------|-------------|--------|------|
| data.csv | Main dataset | API | 2024-01-15 |
```

### Data Versioning

| Method | When to Use | Tools |
|--------|-------------|-------|
| Date in filename | Simple projects | Manual |
| Git LFS | Medium files | Git |
| DVC | Large datasets | DVC |
| Cloud storage | Very large | S3, GCS |

## Reproducibility

### Makefile Example

```makefile
.PHONY: data clean

data: data/processed/dataset.csv

data/processed/dataset.csv: data/raw/original.csv
	python src/data/preprocess.py

clean:
	rm -rf data/processed/*

environment:
	conda env create -f environment.yml
```

### Environment Files

```yaml
# environment.yml
name: project-env
channels:
  - conda-forge
dependencies:
  - python=3.10
  - pandas
  - scikit-learn
  - jupyter
  - pip:
    - custom-package
```

## Anti-Patterns

| Pattern | Problem | Solution |
|---------|---------|----------|
| Numbered versions | `data_v2_final.csv` | Use git/DVC |
| Modifying raw data | Irreproducible | Keep raw immutable |
| Code in notebooks only | Hard to reuse | Extract to `src/` |
| No data documentation | Unknown provenance | README in data/ |
| Flat notebook structure | Hard to follow | Number and organize |

## Scaling Guidelines

### Small -> Medium

Add:
- `src/` for reusable code
- Numbered notebooks
- `Makefile` for automation

### Medium -> Large

Add:
- `experiments/` for tracking
- Workflow manager (Snakemake, Nextflow)
- DVC for data versioning
- CI/CD for tests

### Collaboration

Add:
- `CONTRIBUTING.md`
- `CODEOWNERS`
- Protected branches
- Code review process

## References

- [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)
- [Noble 2009: Quick Guide to Organizing Computational Biology Projects](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)
- [Wilson et al. 2017: Good Enough Practices in Scientific Computing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510)
