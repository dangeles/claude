# Example: Organizing a Mixed Project (Code + Data + Docs)

## Scenario

A user invokes archive-workflow on a machine learning project that combines:
- Python code for ML models
- Jupyter notebooks for exploration
- Large datasets
- Research documentation
- API server code

This is a common "mixed" project type that requires careful organization.

## Initial State

```
ml_project/
├── train.py                    # ML training code
├── model.py                    # Model definitions
├── api.py                      # API server
├── utils.py                    # Utilities
├── test_model.py               # Tests in root
├── exploration.ipynb           # Notebook in root
├── feature_engineering.ipynb   # Notebook in root
├── data.csv                    # 500MB data file (tracked!)
├── model.pkl                   # 200MB model (tracked!)
├── processed_features.parquet  # 150MB features (tracked!)
├── notes.md                    # Unorganized notes
├── paper_draft.md              # Research paper
├── requirements.txt
├── README.md
└── .env                        # Credentials (should not be tracked)
```

## Workflow Execution

### Phase 1: Project Analysis (library-pm)

```
[library-pm] Detecting project type...

Signals detected:
- *.py files with ML imports: +3 (HIGH, Code)
- *.ipynb files: +2 (MEDIUM, Research/Data Science)
- Large data files (*.csv, *.parquet): +2 (HIGH, Data)
- model.pkl: +1 (MEDIUM, ML)
- Paper draft: +1 (MEDIUM, Research)

Classification: MIXED (Code + Data + Research)
Primary: Code (ML)
Secondary: Data Science, Research

Warning: Project has multiple concerns - will use hybrid template.

Quality Gate 1: PASS
```

### Wave 1: Clutter Analysis

**Output: clutter-report.md**

```markdown
# Clutter Report

## Summary
- Total clutter items: 6
- Critical: 4 (large files tracked)
- Moderate: 1
- Low: 1
- Total Clutter Score: 54 (Significant)

## Generated Files (Critical)

| Path | Type | Size | In Git? | Recommendation |
|------|------|------|---------|----------------|
| data.csv | Raw data | 500MB | YES | Remove from git, add to .gitignore |
| model.pkl | Trained model | 200MB | YES | Remove from git, add to .gitignore |
| processed_features.parquet | Features | 150MB | YES | Remove from git, add to .gitignore |
| __pycache__/ | Python cache | 2MB | NO | Add to .gitignore |

## Security Issues (Critical)

| Path | Type | Recommendation |
|------|------|----------------|
| .env | Credentials | Should NOT be tracked - verify gitignore |

## Organizational Mess (Low)

| Path | Issue | Recommendation |
|------|-------|----------------|
| notes.md | Generic name | Rename with topic/date |

## Gitignore Recommendations

```
# Data files (use DVC or similar)
*.csv
*.parquet
*.pkl
data/
models/
!data/.gitkeep
!models/.gitkeep

# Secrets
.env
.env.*
credentials*.json

# Python
__pycache__/
*.py[cod]
.venv/

# Notebooks
.ipynb_checkpoints/
```

## Important Note

This project has 850MB of data/model files tracked in git.
Strongly recommend:
1. Add to .gitignore
2. Use DVC for data versioning
3. Use model registry for model files

Quality Gate 2: PASS
```

### Wave 2: Parallel Analysis

**Output: naming-violations.md**

```markdown
# Naming Violations Report

## Summary
- Files audited: 11
- Violations found: 3
- Warning: 3

## File Naming Violations

| Current | Suggested | Rule | Severity |
|---------|-----------|------|----------|
| notes.md | notes-{topic}.md or {date}-notes.md | Descriptive names | Warning |
| exploration.ipynb | 01-da-exploration.ipynb | Notebook numbering | Warning |
| feature_engineering.ipynb | 02-da-feature-engineering.ipynb | Notebook numbering | Warning |

## Detected Patterns
- Python code: snake_case (consistent, no violations)
- Notebooks: Missing number prefixes
- Docs: Need organization
```

**Output: structure-proposal.md**

```markdown
# Structure Proposal

## Summary
- Project type: Mixed (Code + Data + Research)
- Mode: Adaptive
- Template: Hybrid ML Project

## Proposed Structure

```
ml_project/
├── src/                          # Source code
│   └── ml_project/
│       ├── __init__.py
│       ├── train.py              # Training logic
│       ├── model.py              # Model definitions
│       ├── api.py                # API server
│       └── utils.py              # Utilities
├── tests/                        # Test code
│   └── test_model.py
├── notebooks/                    # Jupyter notebooks
│   ├── 01-da-exploration.ipynb
│   └── 02-da-feature-engineering.ipynb
├── data/                         # Data directory
│   ├── raw/
│   │   └── data.csv             # Raw data (gitignored)
│   ├── processed/
│   │   └── processed_features.parquet  # (gitignored)
│   └── README.md
├── models/                       # Trained models
│   └── model.pkl                # (gitignored)
├── docs/                         # Documentation
│   ├── notes/
│   │   └── notes-project.md
│   └── papers/
│       └── paper_draft.md
├── configs/                      # Configuration
│   └── config.yaml
├── scripts/                      # Utility scripts
│   └── train.sh
├── README.md
├── CLAUDE.md
├── pyproject.toml
├── .env.example                  # Template for .env
├── .gitignore
└── .dvc/                         # DVC config (optional)
```

## Migration Plan

### Phase 1: Create Structure
1. Create src/ml_project/
2. Create tests/
3. Create notebooks/
4. Create data/raw/, data/processed/
5. Create models/
6. Create docs/notes/, docs/papers/
7. Create configs/
8. Create scripts/

### Phase 2: Move Files
| Current | New | Notes |
|---------|-----|-------|
| train.py | src/ml_project/train.py | Core code |
| model.py | src/ml_project/model.py | Core code |
| api.py | src/ml_project/api.py | Core code |
| utils.py | src/ml_project/utils.py | Core code |
| test_model.py | tests/test_model.py | Tests |
| exploration.ipynb | notebooks/01-da-exploration.ipynb | Numbered |
| feature_engineering.ipynb | notebooks/02-da-feature-engineering.ipynb | Numbered |
| data.csv | data/raw/data.csv | Data organization |
| processed_features.parquet | data/processed/features.parquet | Data organization |
| model.pkl | models/model.pkl | Model organization |
| notes.md | docs/notes/notes-project.md | Documentation |
| paper_draft.md | docs/papers/paper_draft.md | Documentation |

### Phase 3: Special Handling
- .env: Keep in root, create .env.example template
- Large files: Update .gitignore, consider DVC

Quality Gate 3: PASS
```

### Wave 3: Expandability Review

**Output: expandability-assessment.md**

```markdown
# Expandability Assessment

## Overall Score: 8/10

**Flag**: ADVISORY with recommendations

## Scalability Assessment

| Factor | Score | Notes |
|--------|-------|-------|
| Code scaling | 9/10 | Proper package structure |
| Data scaling | 7/10 | Needs DVC or similar |
| Model management | 7/10 | Consider MLflow or similar |
| Documentation | 8/10 | Clear separation |

## Modularity Assessment

| Factor | Score | Notes |
|--------|-------|-------|
| Component separation | 9/10 | Clear src/notebooks/data split |
| API isolation | 8/10 | api.py can be own service |
| Config management | 6/10 | Add configs/ directory |

## Recommendations

### High Priority
1. **Add DVC** for data versioning
2. **Add configs/** for ML hyperparameters
3. **Create .env.example** for onboarding

### Medium Priority
4. **Consider MLflow** for experiment tracking
5. **Add Makefile** for common commands
6. **Add docker-compose.yml** for reproducibility

### Low Priority
7. **Add scripts/train.sh** for training automation
8. **Consider splitting api.py** into separate service

## Growth Scenarios

### Scenario: More Data Sources
- data/raw/{source}/ pattern handles this well
- Consider data/external/ for third-party data

### Scenario: Multiple Models
- models/{model_name}/ pattern
- Consider MLflow model registry

### Scenario: Team Growth
- Add CODEOWNERS
- Document in CLAUDE.md

## Conclusion

Hybrid structure is appropriate. Main risks are data management at scale - recommend DVC integration.

Quality Gate 4: PASS
```

### Wave 4: Execution

**User approval**: APPROVE ALL

**Conflict Resolution Applied**:
- Rule 1 (Rename + Move): exploration.ipynb -> notebooks/01-da-exploration.ipynb
- Rule 5 (Clutter Priority): .gitignore updated before moving files

**Output: execution-log.md (summary)**

```markdown
# Execution Log

## Operations Summary

### Directories Created: 11
- src/ml_project/
- tests/
- notebooks/
- data/raw/
- data/processed/
- models/
- docs/notes/
- docs/papers/
- configs/
- scripts/

### Files Moved: 12
All source files moved to src/ml_project/
Tests moved to tests/
Notebooks moved and renamed to notebooks/
Data files moved to data/
Model file moved to models/
Documentation moved to docs/

### Files Created: 6
- src/ml_project/__init__.py
- .gitignore (comprehensive)
- .env.example (template)
- CLAUDE.md
- configs/config.yaml (template)
- data/README.md

### Git Operations
- Removed large files from git tracking
- Files remain on disk, now gitignored
- Preserved git history with git mv

## Final Statistics
- Files processed: 12
- Directories created: 11
- Total disk space now gitignored: 850MB

## Status: COMPLETE

Quality Gate 5: PASS
```

## Final State

```
ml_project/
├── src/
│   └── ml_project/
│       ├── __init__.py
│       ├── train.py
│       ├── model.py
│       ├── api.py
│       └── utils.py
├── tests/
│   └── test_model.py
├── notebooks/
│   ├── 01-da-exploration.ipynb
│   └── 02-da-feature-engineering.ipynb
├── data/
│   ├── raw/
│   │   └── data.csv              # Gitignored
│   ├── processed/
│   │   └── features.parquet      # Gitignored
│   └── README.md
├── models/
│   └── model.pkl                 # Gitignored
├── docs/
│   ├── notes/
│   │   └── notes-project.md
│   └── papers/
│       └── paper_draft.md
├── configs/
│   └── config.yaml
├── scripts/
│   └── .gitkeep
├── README.md                     # Updated
├── CLAUDE.md                     # Created
├── pyproject.toml               # Created (template)
├── .env                         # Kept, gitignored
├── .env.example                 # Created
├── requirements.txt
└── .gitignore                   # Comprehensive
```

## Key Outcomes

1. **Hybrid Structure**: Combines code project (src/tests/) with data science (notebooks/data/) and research (docs/)
2. **Large File Management**: 850MB of data/models now gitignored, DVC recommended
3. **Security**: .env gitignored, .env.example created for onboarding
4. **ML Best Practices**: configs/ for hyperparameters, models/ for artifacts
5. **Clear Separation**: Code vs notebooks vs data vs documentation

## .gitignore (Created)

```gitignore
# Data files
data/raw/*
data/processed/*
!data/*/.gitkeep
!data/*/README.md
*.csv
*.parquet

# Models
models/*
!models/.gitkeep
*.pkl
*.h5
*.onnx

# Secrets
.env
.env.local
credentials*.json

# Python
__pycache__/
*.py[cod]
.venv/
venv/
*.egg-info/
dist/
build/

# Notebooks
.ipynb_checkpoints/

# IDE
.idea/
.vscode/
*.swp

# OS
.DS_Store
```

## .env.example (Created)

```bash
# Database
DATABASE_URL=postgresql://user:password@localhost:5432/db

# API Keys
OPENAI_API_KEY=sk-...

# AWS (if using S3 for data)
AWS_ACCESS_KEY_ID=
AWS_SECRET_ACCESS_KEY=
AWS_DEFAULT_REGION=us-east-1

# MLflow (if using)
MLFLOW_TRACKING_URI=
```

## Next Steps for User

1. **Initialize DVC** (recommended):
   ```bash
   dvc init
   dvc add data/raw/data.csv
   dvc add models/model.pkl
   git add data/raw/data.csv.dvc models/model.pkl.dvc .dvc
   git commit -m "Add DVC tracking for data and models"
   ```

2. **Run tests** to verify moves:
   ```bash
   pytest tests/
   ```

3. **Update imports** in notebooks (may need manual adjustment)
