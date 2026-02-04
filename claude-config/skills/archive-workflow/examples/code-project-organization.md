# Example: Organizing a Python Code Project

## Scenario

A user invokes archive-workflow on a Python codebase with the following issues:
- Files scattered in root directory (main.py, utils.py, helpers.py)
- Tests mixed with source code (test_main.py in root)
- No .gitignore (node_modules from frontend accidentally tracked)
- Inconsistent naming (MyClass.py, other_module.py)
- Stale backup file (backup_main.py.bak)

## Initial State

```
project/
├── main.py              # Source in root
├── utils.py             # Source in root
├── helpers.py           # Source in root
├── MyClass.py           # Wrong naming convention
├── test_main.py         # Test in root
├── backup_main.py.bak   # Stale backup
├── node_modules/        # 342MB, accidentally tracked
├── __pycache__/         # 5MB, not gitignored
├── data.csv             # Data in root
├── README.md
└── requirements.txt
```

## Workflow Execution

### Phase 0: Pre-flight Checks

```
[library-pm] Running pre-flight checks...
- Git status: CLEAN
- Branch: main (not detached)
- Permissions: WRITE OK
- Disk space: 5GB available

Creating session: /tmp/archive-workflow-session-20240115-143022-12345/
```

### Phase 1: Project Analysis (library-pm)

```
[library-pm] Detecting project type...

Signals detected:
- requirements.txt: +3 (HIGH, Python)
- *.py files: +2 (MEDIUM, Python)
- node_modules/: +1 (MEDIUM, JavaScript - secondary)

Classification: CODE (Python) - Primary
Secondary: JavaScript (likely frontend component)

Quality Gate 1: PASS - Project type detected
```

### Wave 1: Clutter Analysis (clutter-analyst)

**Output: clutter-report.md**

```markdown
# Clutter Report

## Summary
- Total clutter items: 4
- Critical: 2 items
- Moderate: 1 item
- Low: 1 item
- Total Clutter Score: 28

## Generated Files (Critical)

| Path | Type | Size | In Git? | Recommendation |
|------|------|------|---------|----------------|
| node_modules/ | npm deps | 342MB | YES | Remove from git, add to .gitignore |
| __pycache__/ | Python cache | 5MB | NO | Add to .gitignore |

## Stale Content (Moderate)

| Path | Last Modified | Days Stale | Reason | Recommendation |
|------|---------------|------------|--------|----------------|
| backup_main.py.bak | 2023-06-15 | 214 | Backup file | Delete |

## Organizational Mess (Low)

| Path | Issue | Recommendation |
|------|-------|----------------|
| data.csv | Data in root | Move to data/ |

## Gitignore Recommendations

Add these patterns:
```
node_modules/
__pycache__/
*.pyc
*.bak
```

Quality Gate 2: PASS - Clutter report complete
```

### Wave 2: Parallel Analysis

**Output: naming-violations.md**

```markdown
# Naming Violations Report

## Summary
- Files audited: 6
- Violations found: 1
- Critical: 0
- Warning: 1

## Detected Patterns
- Dominant file pattern: snake_case (83% of Python files)
- Recommendation: Enforce snake_case

## File Naming Violations

| Current Name | Suggested Name | Rule Violated | Severity | Impact |
|--------------|----------------|---------------|----------|--------|
| MyClass.py | my_class.py | snake_case for Python | Warning | 0 imports |
```

**Output: structure-proposal.md**

```markdown
# Structure Proposal

## Summary
- Project type detected: Code (Python)
- Current depth: 1 level (flat)
- Mode: Adaptive (existing project)
- Impact: 6 files to move

## Current Structure Analysis
- All source files in root (anti-pattern)
- Tests mixed with source
- Data file in root

## Proposed Structure

```
project/
├── src/
│   └── project_name/
│       ├── __init__.py
│       ├── main.py
│       ├── my_class.py
│       ├── utils.py
│       └── helpers.py
├── tests/
│   └── test_main.py
├── data/
│   └── data.csv
├── README.md
├── CLAUDE.md
├── requirements.txt
└── .gitignore
```

## Migration Plan

### Phase 1: Create Structure
1. Create src/project_name/
2. Create tests/
3. Create data/

### Phase 2: Move Files
| Current | New | Impact |
|---------|-----|--------|
| main.py | src/project_name/main.py | 0 imports |
| utils.py | src/project_name/utils.py | 2 imports |
| helpers.py | src/project_name/helpers.py | 1 import |
| MyClass.py | src/project_name/my_class.py | 0 imports (rename + move) |
| test_main.py | tests/test_main.py | 0 refs |
| data.csv | data/data.csv | 1 ref in main.py |

Quality Gate 3: PASS - Both reports complete
```

### Wave 3: Expandability Review (expandability-reviewer)

**Output: expandability-assessment.md**

```markdown
# Expandability Assessment

## Overall Score: 8/10

**Flag**: ADVISORY - Structure is well-suited for growth

## Scalability Assessment

| Factor | Score | Notes |
|--------|-------|-------|
| File growth capacity | 9/10 | src/ structure handles growth well |
| Contributor scaling | 7/10 | Consider CODEOWNERS for teams |
| Module addition | 8/10 | New modules can be added to src/ |
| Data growth | 9/10 | Separate data/ directory |

## Recommendations

### High Priority
1. Add __init__.py to make src/project_name/ a proper package

### Medium Priority
2. Add pyproject.toml for modern packaging

### Low Priority
3. Consider CODEOWNERS if team grows

## Conclusion

The proposed structure is approved. No critical issues.

Quality Gate 4: PASS - Expandability assessment complete
```

### Wave 4: Execution (decision-integrator)

**Pre-Execution: execution-plan.md**

```markdown
# Execution Plan

## Category A: Non-Destructive Operations

### Directories to Create
- src/project_name/
- tests/
- data/

### Files to Move
| Current | New | Reason |
|---------|-----|--------|
| main.py | src/project_name/main.py | Structure |
| utils.py | src/project_name/utils.py | Structure |
| helpers.py | src/project_name/helpers.py | Structure |
| test_main.py | tests/test_main.py | Structure |
| data.csv | data/data.csv | Structure |

### Files to Rename + Move
| Current | New | Reason |
|---------|-----|--------|
| MyClass.py | src/project_name/my_class.py | Naming + Structure |

## Category B: Clutter Cleanup

### Gitignore Patterns to Add
```
node_modules/
__pycache__/
*.pyc
*.bak
.venv/
```

### Files to Remove from Git
- node_modules/ (will remain on disk, just untracked)

## Category C: Deletions (APPROVAL REQUIRED)

| File | Size | Reason | Approval |
|------|------|--------|----------|
| backup_main.py.bak | 2KB | Stale backup (214 days) | [ ] Approve |

## User Action Required
Please respond: APPROVE ALL / APPROVE WITH EXCLUSIONS / REJECT
```

**User Response**: APPROVE ALL

**Post-Execution: execution-log.md**

```markdown
# Execution Log

## Session: archive-workflow-session-20240115-143022-12345

## Execution Started
- Timestamp: 2024-01-15T14:32:15Z
- Rollback point: abc123f

## Operations Performed

### Gitignore Updates
- Created .gitignore
- Added 5 patterns

### Git Operations
- Removed node_modules/ from git tracking
- Files remain on disk but untracked

### Directories Created
| Directory | Status |
|-----------|--------|
| src/project_name/ | SUCCESS |
| tests/ | SUCCESS |
| data/ | SUCCESS |

### Files Moved
| Original | New | Status |
|----------|-----|--------|
| main.py | src/project_name/main.py | SUCCESS |
| utils.py | src/project_name/utils.py | SUCCESS |
| helpers.py | src/project_name/helpers.py | SUCCESS |
| test_main.py | tests/test_main.py | SUCCESS |
| data.csv | data/data.csv | SUCCESS |

### Files Renamed + Moved
| Original | New | Status |
|----------|-----|--------|
| MyClass.py | src/project_name/my_class.py | SUCCESS |

### Files Deleted
| File | Status | Confirmed By |
|------|--------|--------------|
| backup_main.py.bak | SUCCESS | User approval |

### Files Created
| File | Status |
|------|--------|
| src/project_name/__init__.py | SUCCESS |
| .gitignore | SUCCESS |
| CLAUDE.md | SUCCESS |

### Documentation Updated
- README.md: Updated structure section
- CLAUDE.md: Created with project context

## Final Statistics
- Files processed: 8
- Files moved: 6
- Files renamed: 1
- Files deleted: 1
- Directories created: 3
- Gitignore patterns added: 5

## Status: COMPLETE

Quality Gate 5: PASS - All operations successful
```

## Final State

```
project/
├── src/
│   └── project_name/
│       ├── __init__.py       # Created
│       ├── main.py           # Moved from root
│       ├── my_class.py       # Renamed + moved from MyClass.py
│       ├── utils.py          # Moved from root
│       └── helpers.py        # Moved from root
├── tests/
│   └── test_main.py          # Moved from root
├── data/
│   └── data.csv              # Moved from root
├── node_modules/             # Still exists, now gitignored
├── __pycache__/              # Still exists, now gitignored
├── README.md                 # Updated
├── CLAUDE.md                 # Created
├── requirements.txt
└── .gitignore                # Created with patterns
```

## Key Outcomes

1. **Structure**: Proper Python package structure with src/, tests/, data/
2. **Naming**: All files now follow snake_case convention
3. **Gitignore**: Generated files excluded, secrets template added
4. **Documentation**: CLAUDE.md created, README.md updated
5. **Cleanup**: Stale backup removed, clutter gitignored
6. **Git History**: File moves use git mv, history preserved

## Rollback (If Needed)

```bash
# If issues found after execution:
git reset --hard abc123f

# Verify restoration:
git status
ls -la  # Should show original flat structure
```
