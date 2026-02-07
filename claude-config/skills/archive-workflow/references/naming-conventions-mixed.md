# Naming Conventions: Mixed Projects

## Overview

This reference defines naming conventions for mixed-type projects that combine code, research, and/or data elements. When a project spans multiple types, this template provides a unified set of conventions that prioritize the PRIMARY type while accommodating secondary concerns.

## Primary Type Priority

Mixed projects should identify their primary type and use that type's conventions as the default. When conflicts arise between type conventions:

1. **Code-primary mixed**: Use code conventions (snake_case for .py) as default; research conventions for docs
2. **Research-primary mixed**: Use research conventions (kebab-case for .md) as default; code conventions for source
3. **Data-primary mixed**: Use data conventions (snake_case with dates) as default; code conventions for pipeline scripts

## Source Code Files

Inherited from code conventions:

| Element | Convention | Example | Notes |
|---------|------------|---------|-------|
| Python modules | snake_case | `my_module.py` | Lowercase with underscores |
| Python packages | snake_case | `my_package/` | Same as modules |
| Tests | `test_*.py` | `test_utils.py` | pytest convention |
| Scripts | snake_case | `run_analysis.py` | Descriptive verbs |
| JS/TS components | PascalCase | `UserProfile.tsx` | React convention |
| JS/TS utilities | camelCase or kebab-case | `apiClient.ts` | Project-dependent |

## Documentation Files

Inherited from research conventions:

| Document Type | Convention | Example |
|---------------|------------|---------|
| Literature review | `review-{topic}.md` | `review-oxygen-solubility.md` |
| Paper notes | `{author}-{year}-{topic}.md` | `smith-2024-bioreactors.md` |
| Analysis document | `analysis-{question}.md` | `analysis-growth-rates.md` |
| Meeting notes | `{YYYY-MM-DD}-meeting-{topic}.md` | `2024-01-15-meeting-advisor.md` |
| General docs | kebab-case | `api-reference.md` |

## Data Files

Inherited from data conventions:

| Element | Convention | Example |
|---------|------------|---------|
| Raw data | `{source}-{YYYYMMDD}-{description}.{ext}` | `sensor-20240115-readings.csv` |
| Processed data | `{domain}-{description}-{YYYYMMDD}.{ext}` | `sales-summary-20240115.parquet` |
| Metadata | `{datafile}.meta.yaml` | `readings.csv.meta.yaml` |

## Directory Naming

| Context | Convention | Example |
|---------|------------|---------|
| Source code dirs | snake_case (Python) or lowercase (general) | `my_package/`, `src/` |
| Documentation dirs | lowercase | `docs/`, `reports/` |
| Data dirs | lowercase | `data/raw/`, `data/processed/` |
| Test dirs | lowercase | `tests/unit/`, `tests/integration/` |

## Notebooks

Inherited from research conventions:

| Element | Convention | Example |
|---------|------------|---------|
| Exploration | `{NN}-{initials}-{description}.ipynb` | `01-da-data-exploration.ipynb` |
| Analysis | `{NN}-{initials}-{analysis-type}.ipynb` | `02-da-statistical-analysis.ipynb` |

## Configuration Files

| File | Convention |
|------|------------|
| `.gitignore` | Exact, dotfile |
| `.env` | Exact, dotfile |
| `README.md` | Exact, uppercase |
| `CHANGELOG.md` | Exact, uppercase |
| `CLAUDE.md` | Exact, uppercase |
| Environment configs | `config-{env}.yaml` |

## Anti-Patterns

### Critical (Must Fix)

| Pattern | Problem | Solution |
|---------|---------|----------|
| Spaces in names | Breaks shell commands | Use underscores or hyphens |
| Special characters | Encoding issues | Alphanumeric only |
| Starting with number | Breaks Python imports | Prefix with word |
| Mixed case at same level | Confusion | Pick one convention per directory |
| Version numbers in filenames | Maintenance burden | Use git for versioning |

### Context-Specific Rules

| File Location | Convention Priority | Rationale |
|---------------|-------------------|-----------|
| `src/` or code directories | Code conventions | Import compatibility |
| `docs/` or documentation | Research conventions | Readability |
| `data/` directories | Data conventions | Date-based sorting |
| `notebooks/` | Research conventions | Numbered ordering |
| `scripts/` | Code conventions | Shell compatibility |

## Adaptive Mode

When analyzing existing mixed projects:

1. **Identify dominant type** by file count and activity
2. **Check each directory's conventions** independently
3. **If >80% consistent within a directory**: Enforce that pattern
4. **If inconsistent**: Recommend the type-appropriate default for that directory
5. **Document detected patterns** in report

## Conflict Resolution

When code and research conventions conflict:
- **In `src/`**: Code conventions always win (snake_case for .py)
- **In `docs/`**: Research conventions always win (kebab-case for .md)
- **In root directory**: Follow the primary project type
- **If unclear**: Default to the convention that preserves import/build compatibility

## References

- naming-conventions-code.md (for code-specific rules)
- naming-conventions-research.md (for research-specific rules)
- naming-conventions-data.md (for data-specific rules)
