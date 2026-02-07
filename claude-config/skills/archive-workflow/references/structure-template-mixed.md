# Structure Template: Mixed Projects

## Overview

This reference defines directory structure templates for mixed-type projects that combine code, research, and/or data elements. The mixed template merges relevant directories from the code, research, and data templates, adapting based on the project's primary and secondary types.

## Standard Mixed Structure

### Full Template (Code + Research)

```
project/
├── src/                      # Source code (from code template)
│   └── {package_name}/
│       ├── __init__.py
│       ├── main.py
│       ├── models/
│       ├── services/
│       └── utils/
├── tests/                    # Test code (from code template)
│   ├── unit/
│   └── integration/
├── docs/                     # Documentation (from research template)
│   ├── literature/           # Paper notes, reviews
│   ├── analysis/             # Analysis documents
│   └── guides/               # User guides, API docs
├── data/                     # Data files (from data template)
│   ├── raw/                  # Immutable original data
│   ├── processed/            # Transformed data
│   └── external/             # Third-party data
├── notebooks/                # Jupyter notebooks (from research template)
│   ├── 01-exploration.ipynb
│   └── 02-analysis.ipynb
├── scripts/                  # Utility scripts (from code template)
│   ├── setup.sh
│   └── deploy.sh
├── reports/                  # Generated outputs (from research template)
│   ├── figures/
│   └── drafts/
├── configs/                  # Configuration (from data template)
│   ├── config-dev.yaml
│   └── config-prod.yaml
├── pyproject.toml            # Package config
├── README.md
├── CLAUDE.md
├── CHANGELOG.md
└── .gitignore
```

### Full Template (Research + Data)

```
project/
├── data/                     # Primary data (from data template)
│   ├── raw/
│   │   ├── {source}/
│   │   └── README.md
│   ├── processed/
│   ├── features/             # ML features (if applicable)
│   └── external/
├── notebooks/                # Analysis notebooks (from research template)
│   ├── exploration/
│   └── modeling/
├── docs/                     # Documentation (from research template)
│   ├── literature/
│   ├── analysis/
│   └── protocols/
├── src/                      # Supporting code (from code template)
│   └── {package_name}/
├── pipelines/                # ETL/processing (from data template)
│   ├── extract/
│   ├── transform/
│   └── load/
├── models/                   # Trained models (from research/data)
├── reports/                  # Generated outputs
│   ├── figures/
│   └── tables/
├── tests/
├── configs/
├── README.md
├── CLAUDE.md
└── .gitignore
```

### Minimal Mixed Structure

```
project/
├── src/                      # Source code
├── tests/                    # Tests
├── docs/                     # Documentation
├── data/                     # Data files
│   ├── raw/
│   └── processed/
├── notebooks/                # Notebooks (if applicable)
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Directory Purposes

### Core Directories (Always Present)

| Directory | Source Template | Purpose | Contents |
|-----------|---------------|---------|----------|
| `src/` | Code | Source code | Main application code |
| `tests/` | Code | Test code | Unit, integration tests |
| `docs/` | Research | Documentation | Guides, literature, analysis |
| `data/` | Data | Data files | Raw, processed, external |

### Optional Directories (Based on Secondary Type)

| Directory | When to Include | Purpose |
|-----------|----------------|---------|
| `notebooks/` | Research or data component | Jupyter/R notebooks |
| `scripts/` | Automation needed | Build, deploy, utility scripts |
| `reports/` | Research output needed | Figures, tables, drafts |
| `models/` | ML component | Trained models |
| `pipelines/` | ETL/data processing | Extract, transform, load |
| `configs/` | Multiple environments | Environment configs |
| `experiments/` | Experimental work | Tracking experiments |
| `references/` | Reference materials | Data dictionaries, manuals |

## Adaptation Rules

### Based on `project.secondary_type`

| Primary | Secondary | Key Adaptations |
|---------|-----------|----------------|
| code | research | Add `docs/literature/`, `notebooks/`, `reports/` |
| code | data | Add `data/`, `pipelines/`, `configs/` |
| research | code | Add `src/`, `tests/`, `scripts/` |
| research | data | Add `data/raw/`, `pipelines/`, `models/` |
| data | code | Elevate `src/` and `tests/` to required |
| data | research | Add `docs/literature/`, `notebooks/`, `reports/` |

### When Both Secondary Types Apply

If a project is truly mixed across all three types, include all core directories:
- `src/`, `tests/` (code)
- `docs/`, `notebooks/`, `reports/` (research)
- `data/`, `pipelines/`, `configs/` (data)

## Required Files

| File | Purpose | Always Required |
|------|---------|----------------|
| `README.md` | Project overview | Yes |
| `CLAUDE.md` | AI assistant context | Yes |
| `.gitignore` | Git exclusions | Yes |
| `pyproject.toml` | Package config | If Python code present |
| `package.json` | Package config | If Node.js code present |

## Generated Directories (gitignored)

| Directory | Purpose |
|-----------|---------|
| `dist/`, `build/` | Compiled output |
| `__pycache__/` | Python bytecode |
| `node_modules/` | Node dependencies |
| `.venv/` | Virtual environment |
| `target/` | Rust compiled output |

## Scaling Guidelines

### Small Mixed -> Medium Mixed

Add:
- Proper test organization (`tests/unit/`, `tests/integration/`)
- Documentation structure (`docs/api/`, `docs/guides/`)
- Numbered notebooks
- Scripts directory for automation

### Medium Mixed -> Large Mixed

Add:
- Feature-based organization in `src/`
- Workflow manager for data pipelines
- `CODEOWNERS` for ownership
- `docs/adr/` for architecture decisions
- CI/CD configuration

## Anti-Patterns

| Pattern | Problem | Solution |
|---------|---------|----------|
| Flat structure (everything in root) | Hard to find files | Use type-appropriate subdirectories |
| Code mixed with data | Version control issues | Separate `src/` from `data/` |
| Notebooks as only code | Hard to reuse | Extract reusable code to `src/` |
| No data documentation | Unknown provenance | README in `data/raw/` |
| Deeply nested (>4 levels) | Hard to navigate | Flatten or split by feature |

## References

- structure-template-code.md (for code-specific structures)
- structure-template-research.md (for research-specific structures)
- structure-template-data.md (for data-specific structures)
