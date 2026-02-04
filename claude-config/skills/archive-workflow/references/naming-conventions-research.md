# Naming Conventions: Research Projects

## Overview

This reference defines naming conventions for research projects (scientific research, literature reviews, data analysis, academic work). Based on best practices from scientific computing and research data management.

## Document Types

### Literature and Reviews

| Document Type | Convention | Example |
|---------------|------------|---------|
| Literature review | `review-{topic}.md` | `review-oxygen-solubility.md` |
| Paper notes | `{author}-{year}-{topic}.md` | `smith-2024-bioreactors.md` |
| Annotated bibliography | `bibliography-{topic}.md` | `bibliography-membrane-transport.md` |
| Reading list | `reading-{topic}.md` | `reading-machine-learning.md` |

### Analysis Documents

| Document Type | Convention | Example |
|---------------|------------|---------|
| Analysis document | `analysis-{question}.md` | `analysis-growth-rates.md` |
| Methodology notes | `methodology-{technique}.md` | `methodology-pcr.md` |
| Protocol | `protocol-{procedure}.md` | `protocol-cell-culture.md` |
| Comparison | `comparison-{items}.md` | `comparison-algorithms.md` |

### Planning Documents

| Document Type | Convention | Example |
|---------------|------------|---------|
| Experiment plan | `{YYYY-MM-DD}-{title}.md` | `2024-01-15-fermentation-test.md` |
| Project plan | `plan-{project}.md` | `plan-phd-timeline.md` |
| Meeting notes | `{YYYY-MM-DD}-meeting-{topic}.md` | `2024-01-15-meeting-advisor.md` |
| Lab notebook entry | `{YYYY-MM-DD}-notebook-{topic}.md` | `2024-01-15-notebook-experiment-3.md` |

## Data Files

### Raw Data

| Element | Convention | Example |
|---------|------------|---------|
| Raw data files | `{YYYYMMDD}-{source}-{description}.{ext}` | `20240115-sensor-temp-readings.csv` |
| Location | `data/raw/` | Never modified after collection |
| Metadata | `{datafile}-metadata.yaml` | `20240115-sensor-temp-readings-metadata.yaml` |

### Processed Data

| Element | Convention | Example |
|---------|------------|---------|
| Processed files | `{YYYYMMDD}-{description}-processed.{ext}` | `20240115-temp-readings-processed.csv` |
| Location | `data/processed/` | Results of transformation |
| Intermediate | `data/interim/` | In-progress transformations |

### External Data

| Element | Convention | Example |
|---------|------------|---------|
| External sources | `{source}-{YYYYMMDD}-{description}.{ext}` | `ncbi-20240115-gene-sequences.fasta` |
| Location | `data/external/` | Third-party data |
| Attribution | Include README with source | Document provenance |

## Notebooks

### Jupyter Notebooks

| Element | Convention | Example |
|---------|------------|---------|
| Exploration | `{NN}-{initials}-{description}.ipynb` | `01-da-data-exploration.ipynb` |
| Analysis | `{NN}-{initials}-{analysis-type}.ipynb` | `02-da-statistical-analysis.ipynb` |
| Visualization | `{NN}-{initials}-{viz-type}.ipynb` | `03-da-figure-generation.ipynb` |

**Numbering Rules**:
- Start with two-digit number: `01-`, `02-`, etc.
- Numbers indicate intended execution order
- Initials help track authorship in collaborative projects
- Description should be kebab-case

### R Markdown

| Element | Convention | Example |
|---------|------------|---------|
| Analysis | `{NN}-{description}.Rmd` | `01-data-cleaning.Rmd` |
| Report | `report-{topic}.Rmd` | `report-final-results.Rmd` |

## Figures and Outputs

### Figures

| Element | Convention | Example |
|---------|------------|---------|
| Generated figures | `fig-{NN}-{description}.{ext}` | `fig-01-growth-curve.png` |
| Publication figures | `figure-{N}-{description}.{ext}` | `figure-1-main-result.pdf` |
| Location | `reports/figures/` | All generated graphics |

### Tables

| Element | Convention | Example |
|---------|------------|---------|
| Data tables | `table-{NN}-{description}.csv` | `table-01-summary-stats.csv` |
| Publication tables | `table-{N}-{description}.tex` | `table-1-comparison.tex` |

### Reports

| Element | Convention | Example |
|---------|------------|---------|
| Draft reports | `draft-{YYYY-MM-DD}-{title}.md` | `draft-2024-01-15-preliminary-results.md` |
| Final reports | `report-{title}.pdf` | `report-final-analysis.pdf` |
| Location | `reports/` | All report outputs |

## Directory Structure

### Standard Research Directories

```
project/
├── data/
│   ├── raw/              # Immutable original data
│   ├── interim/          # Intermediate transformations
│   ├── processed/        # Final datasets
│   └── external/         # Third-party data
├── notebooks/            # Numbered Jupyter/R notebooks
├── docs/
│   ├── literature/       # Paper notes, reviews
│   └── analysis/         # Analysis documents
├── models/               # Trained models (if ML)
├── references/           # Data dictionaries, manuals
├── reports/
│   ├── figures/          # Generated graphics
│   └── drafts/           # Draft reports
├── src/                  # Source code (if any)
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Date Formats

### When to Use Which Format

| Context | Format | Example |
|---------|--------|---------|
| Filenames (sortable) | `YYYYMMDD` | `20240115` |
| Document titles | `YYYY-MM-DD` | `2024-01-15` |
| Human-readable | `Month DD, YYYY` | `January 15, 2024` |

### Why YYYYMMDD for Files

- Sorts chronologically in file listings
- Unambiguous (no US vs EU confusion)
- Works in filenames (no slashes)

## Versioning

### Preferred: No Versions in Filenames

Use git for versioning instead of:
- ~~`analysis-v2.md`~~
- ~~`data-final-FINAL.csv`~~
- ~~`report-draft-2.docx`~~

### When Versions Are Necessary

For non-git-tracked files or external sharing:
- Use date: `report-20240115.pdf`
- Or semantic: `protocol-1.2.md` (major.minor)

## Special Characters

### Allowed
- Alphanumeric: `a-z`, `A-Z`, `0-9`
- Hyphens: `-` (preferred word separator)
- Underscores: `_` (Python compatibility)
- Periods: `.` (extensions only)

### Not Allowed
- Spaces (use hyphens)
- Special characters: `!@#$%^&*()`
- Slashes: `/\`
- Colons: `:` (Windows incompatible)

## Anti-Patterns

| Bad | Good | Reason |
|-----|------|--------|
| `final_FINAL_v2.docx` | `report-2024-01-15.docx` | Use dates, not "final" |
| `data (copy).csv` | `data-backup-20240115.csv` | No spaces or parentheses |
| `New Folder/` | `analysis-q1/` | Descriptive, no spaces |
| `Paper Notes.md` | `paper-notes-smith-2024.md` | Lowercase, author/date |
| `Figure 1.png` | `fig-01-growth-curve.png` | No spaces, numbered |

## Adaptive Mode

When analyzing existing research projects:

1. **Check for established conventions** (lab/institution standards)
2. **Respect domain patterns** (bioinformatics, physics, etc.)
3. **Preserve provenance** in data file names
4. **Maintain notebook numbering** if already established

## References

- [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)
- [Harvard Research Data Management](https://researchdatamanagement.harvard.edu/)
- [Carnegie Mellon File Naming Guidelines](https://www.library.cmu.edu/datapub/file-naming)
