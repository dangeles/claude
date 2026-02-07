# Example: Organizing a Research Project

## Scenario

A user invokes archive-workflow on a bioinformatics research project with:
- Unorganized data files in root
- Paper notes with inconsistent naming
- Jupyter notebooks without numbering
- No clear separation between raw and processed data
- Missing documentation

## Initial State

```
research_project/
├── data.csv                    # Raw data in root
├── processed_data.csv          # Processed data in root
├── Smith2024_paper_notes.md    # Inconsistent naming
├── review oxygen solubility.md # Spaces in filename
├── analysis.ipynb              # No numbering
├── more_analysis.ipynb         # No numbering
├── figures.ipynb               # No numbering
├── experiment_20240110.csv     # Data mixed with code
├── old_analysis.ipynb          # Stale notebook
├── README.md
└── requirements.txt
```

## Workflow Execution

### Phase 1: Project Analysis (library-pm)

```
[library-pm] Detecting project type...

Signals detected:
- .ipynb files: +3 (HIGH, Research/Data Science)
- .csv files: +2 (MEDIUM, Data)
- Paper notes (.md with author-year): +2 (HIGH, Research)
- No src/ or tests/: -2 (Not primarily code)

Classification: RESEARCH
Secondary: Data Science

Quality Gate 1: PASS
```

### Wave 1: Clutter Analysis

```
[Wave 1/4 - Clutter Analysis] Dispatching clutter-analyst via Task tool...

Task tool:
  Description: "Clutter analyst: Scan project for generated files, stale content, and mess"
  archival_context: "skip"
```

**Output: clutter-report.md**

```markdown
# Clutter Report

## Summary
- Total clutter items: 2
- Critical: 0
- Moderate: 1
- Low: 1
- Total Clutter Score: 9

## Stale Content (Moderate)

| Path | Last Modified | Days Stale | Reason | Recommendation |
|------|---------------|------------|--------|----------------|
| old_analysis.ipynb | 2023-08-15 | 153 | Prefix "old_" | Archive or delete |

## Organizational Mess (Low)

| Path | Issue | Recommendation |
|------|-------|----------------|
| review oxygen solubility.md | Spaces in filename | Rename to kebab-case |

## Gitignore Recommendations

Add:
```
.ipynb_checkpoints/
*.pyc
__pycache__/
.DS_Store
```

Quality Gate 2: PASS
```

### Wave 2: Parallel Analysis

```
[Wave 2/4 - Parallel Organization] Dispatching nomenclature-enforcer and structure-organizer simultaneously...

Task tool (Agent 1):
  Description: "Nomenclature enforcer: Audit file and directory naming conventions"
  archival_context: "skip"

Task tool (Agent 2):
  Description: "Structure organizer: Propose directory structure for research project"
  archival_context: "skip"
```

**Output: naming-violations.md**

```markdown
# Naming Violations Report

## Summary
- Files audited: 9
- Violations found: 4
- Critical: 1 (spaces)
- Warning: 3

## File Naming Violations

| Current Name | Suggested Name | Rule Violated | Severity |
|--------------|----------------|---------------|----------|
| review oxygen solubility.md | review-oxygen-solubility.md | No spaces | Critical |
| Smith2024_paper_notes.md | smith-2024-paper-notes.md | Author-year format | Warning |
| analysis.ipynb | 01-da-analysis.ipynb | Notebook numbering | Warning |
| more_analysis.ipynb | 02-da-more-analysis.ipynb | Notebook numbering | Warning |
| figures.ipynb | 03-da-figures.ipynb | Notebook numbering | Warning |

## Detected Patterns
- Data files: Use date prefix (YYYYMMDD)
- Paper notes: Use author-year format (lowercase-hyphen)
- Notebooks: Add number prefix
```

**Output: structure-proposal.md**

```markdown
# Structure Proposal

## Summary
- Project type: Research (Bioinformatics)
- Mode: Adaptive
- Impact: 9 files to reorganize

## Proposed Structure (Cookiecutter Data Science)

```
research_project/
├── data/
│   ├── raw/                  # Immutable original data
│   │   └── 20240110-experiment.csv
│   ├── processed/            # Transformed data
│   │   └── 20240110-experiment-processed.csv
│   └── external/             # Third-party data
├── notebooks/                # Numbered Jupyter notebooks
│   ├── 01-da-data-exploration.ipynb
│   ├── 02-da-analysis.ipynb
│   └── 03-da-figures.ipynb
├── docs/
│   ├── literature/           # Paper notes
│   │   ├── smith-2024-bioreactors.md
│   │   └── review-oxygen-solubility.md
│   └── analysis/             # Analysis documents
├── reports/
│   └── figures/              # Generated figures
├── archive/                  # Old/stale content (optional)
│   └── old_analysis.ipynb
├── README.md
├── CLAUDE.md
└── requirements.txt
```

## Migration Plan

### Phase 1: Create Structure
1. Create data/raw/, data/processed/, data/external/
2. Create notebooks/
3. Create docs/literature/, docs/analysis/
4. Create reports/figures/
5. Create archive/ (optional)

### Phase 2: Move Files
| Current | New | Notes |
|---------|-----|-------|
| data.csv | data/raw/20240110-data.csv | Date prefix added |
| processed_data.csv | data/processed/20240110-data-processed.csv | Standardized name |
| experiment_20240110.csv | data/raw/20240110-experiment.csv | Date moved to prefix |
| analysis.ipynb | notebooks/01-da-analysis.ipynb | Numbered |
| more_analysis.ipynb | notebooks/02-da-more-analysis.ipynb | Numbered |
| figures.ipynb | notebooks/03-da-figures.ipynb | Numbered |
| Smith2024_paper_notes.md | docs/literature/smith-2024-bioreactors.md | Renamed |
| review oxygen solubility.md | docs/literature/review-oxygen-solubility.md | Spaces removed |
| old_analysis.ipynb | archive/old_analysis.ipynb | Preserved |

Quality Gate 3: PASS
```

### Wave 3: Expandability Review

```
[Wave 3/4 - Expandability Review] Dispatching expandability-reviewer...

Task tool:
  Description: "Expandability reviewer: Assess scalability and modularity of proposed structure"
  archival_context: "skip"
```

**Output: expandability-assessment.md**

```markdown
# Expandability Assessment

## Overall Score: 8/10

**Flag**: ADVISORY

## Scalability Assessment

| Factor | Score | Notes |
|--------|-------|-------|
| Data growth | 9/10 | Partitioned data dirs handle growth |
| Notebook scaling | 8/10 | Numbered notebooks maintain order |
| Literature organization | 8/10 | Author-year format scales well |
| Collaboration | 7/10 | Consider initials in notebook names |

## Recommendations

1. **Add data/raw/README.md** - Document data provenance
2. **Consider DVC** - For large data versioning
3. **Add environment.yml** - Reproducible environment

## Conclusion

Proposed structure follows Cookiecutter Data Science pattern. Well-suited for research reproducibility.

Quality Gate 4: PASS
```

### Wave 4: Execution

```
[Wave 4/4 - Synthesis & Execution] Dispatching decision-integrator...

Task tool:
  Description: "Decision integrator: Merge analyst reports, execute approved operations"
  archival_context: "skip"
```

**User approval**: APPROVE WITH EXCLUSIONS
- Exclusion: Keep old_analysis.ipynb in archive/ (don't delete)

**Output: execution-log.md**

```markdown
# Execution Log

## Operations Performed

### Directories Created
- data/raw/ - SUCCESS
- data/processed/ - SUCCESS
- data/external/ - SUCCESS
- notebooks/ - SUCCESS
- docs/literature/ - SUCCESS
- docs/analysis/ - SUCCESS
- reports/figures/ - SUCCESS
- archive/ - SUCCESS

### Files Moved and Renamed
| Original | New | Status |
|----------|-----|--------|
| data.csv | data/raw/20240110-data.csv | SUCCESS |
| processed_data.csv | data/processed/20240110-data-processed.csv | SUCCESS |
| experiment_20240110.csv | data/raw/20240110-experiment.csv | SUCCESS |
| analysis.ipynb | notebooks/01-da-analysis.ipynb | SUCCESS |
| more_analysis.ipynb | notebooks/02-da-more-analysis.ipynb | SUCCESS |
| figures.ipynb | notebooks/03-da-figures.ipynb | SUCCESS |
| Smith2024_paper_notes.md | docs/literature/smith-2024-bioreactors.md | SUCCESS |
| review oxygen solubility.md | docs/literature/review-oxygen-solubility.md | SUCCESS |
| old_analysis.ipynb | archive/old_analysis.ipynb | SUCCESS |

### Files Created
| File | Status |
|------|--------|
| data/raw/README.md | SUCCESS |
| .gitignore | SUCCESS |
| CLAUDE.md | SUCCESS |
| .gitkeep files | SUCCESS |

### Documentation Updated
- README.md: Updated with structure overview
- CLAUDE.md: Created with project context
- data/raw/README.md: Created with provenance template

## Final Statistics
- Files processed: 9
- Files moved: 9
- Directories created: 8
- Documentation files created: 3

## Status: COMPLETE

Quality Gate 5: PASS
```

## Final State

```
research_project/
├── data/
│   ├── raw/
│   │   ├── README.md               # Data provenance
│   │   ├── 20240110-data.csv       # Renamed with date
│   │   └── 20240110-experiment.csv # Renamed with date
│   ├── processed/
│   │   └── 20240110-data-processed.csv
│   └── external/
│       └── .gitkeep
├── notebooks/
│   ├── 01-da-analysis.ipynb        # Numbered
│   ├── 02-da-more-analysis.ipynb   # Numbered
│   └── 03-da-figures.ipynb         # Numbered
├── docs/
│   ├── literature/
│   │   ├── smith-2024-bioreactors.md    # Author-year format
│   │   └── review-oxygen-solubility.md  # Spaces removed
│   └── analysis/
│       └── .gitkeep
├── reports/
│   └── figures/
│       └── .gitkeep
├── archive/
│   └── old_analysis.ipynb          # Preserved
├── README.md                       # Updated
├── CLAUDE.md                       # Created
├── requirements.txt
└── .gitignore                      # Created
```

## Key Outcomes

1. **Data Management**: Clear raw/processed separation with date-prefixed names
2. **Reproducibility**: Numbered notebooks show execution order
3. **Literature Organization**: Consistent author-year naming
4. **Documentation**: README.md updated, CLAUDE.md created, data README added
5. **Archive Strategy**: Old content preserved rather than deleted
6. **Gitignore**: Research-appropriate patterns added

## data/raw/README.md Template

```markdown
# Raw Data

## Provenance

| File | Source | Date | Description |
|------|--------|------|-------------|
| 20240110-data.csv | Lab instrument | 2024-01-10 | Temperature readings |
| 20240110-experiment.csv | Manual entry | 2024-01-10 | Experiment parameters |

## Data Dictionary

See docs/analysis/ for variable definitions.

## Notes

- Raw data is IMMUTABLE - never modify these files
- Process data using notebooks/ and save to data/processed/
```
