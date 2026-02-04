# Clutter Detection Rules

## Overview

This reference defines rules for detecting clutter in projects. Clutter is organized into three categories by severity: Critical (generated files), Moderate (stale content), and Low (organizational mess).

## Category 1: Generated Files (Critical Priority)

### Action: Add to .gitignore (never track)

### Detection Patterns

| Pattern | Type | Typical Size | Confidence |
|---------|------|--------------|------------|
| `node_modules/` | npm dependencies | 100MB-1GB | HIGH |
| `__pycache__/` | Python bytecode | 1-50MB | HIGH |
| `.venv/`, `venv/`, `ENV/` | Python virtualenv | 100MB-500MB | HIGH |
| `build/`, `dist/` | Build outputs | 10MB-500MB | HIGH |
| `*.pyc`, `*.pyo` | Python compiled | <1MB each | HIGH |
| `.next/`, `.nuxt/` | Framework builds | 50MB-200MB | HIGH |
| `target/` | Rust/Java builds | 100MB-1GB | HIGH |
| `.pytest_cache/` | Test cache | 1-10MB | HIGH |
| `.mypy_cache/` | Type checker cache | 1-10MB | HIGH |
| `coverage/`, `htmlcov/` | Coverage reports | 5-50MB | MEDIUM |
| `.tox/`, `.nox/` | Test environments | 50-200MB | HIGH |
| `.eggs/`, `*.egg-info/` | Python packaging | 1-10MB | HIGH |
| `*.egg` | Python eggs | 1-50MB | HIGH |
| `.sass-cache/` | Sass compilation | 1-10MB | MEDIUM |
| `.parcel-cache/` | Parcel bundler | 10-100MB | HIGH |
| `*.so`, `*.dylib`, `*.dll` | Compiled extensions | 1-50MB | HIGH |
| `.terraform/` | Terraform providers | 100MB-500MB | HIGH |

### Special Cases

| Pattern | Detection Rule | Action |
|---------|----------------|--------|
| `node_modules/` in git | Size > 100MB | Add to gitignore, remove from git |
| `venv/` tracked | Has `pyvenv.cfg` | Add to gitignore, remove from git |
| `build/` with artifacts | Contains `.o`, `.a` files | Add to gitignore |

### Confidence Modifiers

- **+HIGH**: Pattern matches exactly known generated path
- **+MEDIUM**: Pattern matches common generated file extension
- **-LOW**: File is small (<1MB) and may be intentional

## Category 2: Stale Content (Moderate Priority)

### Action: Recommend archive or delete (user confirmation required)

### Detection Heuristics

| Signal | Threshold | Confidence | Detection Method |
|--------|-----------|------------|------------------|
| Last commit | >90 days | MEDIUM | `git log -1 --format=%cd -- {file}` |
| Merged branch | Still present | HIGH | `git branch --merged main` |
| Naming pattern | `wip-*`, `draft-*`, `old-*` | MEDIUM | Filename regex |
| Naming pattern | `*-backup`, `*.bak`, `*.old` | MEDIUM | Filename regex |
| Naming pattern | `*-copy`, `* (1)`, `* (2)` | MEDIUM | Filename regex |
| Orphan file | Not imported/referenced | LOW | Import analysis |
| Abandoned experiment | In `/experiments/` with no commits | MEDIUM | Path + git log |

### Stale File Detection

```python
# Pseudocode for stale detection
def is_stale(file_path):
    # Check last modification
    last_commit_date = git_log_date(file_path)
    if (today - last_commit_date).days > 90:
        return True, "No commits in 90+ days"

    # Check naming patterns
    stale_patterns = [
        r'^wip[-_]',
        r'^draft[-_]',
        r'^old[-_]',
        r'[-_]backup$',
        r'\.bak$',
        r'\.old$',
        r' \(\d+\)$',  # "file (1).txt"
    ]
    for pattern in stale_patterns:
        if re.match(pattern, file_path.name):
            return True, f"Matches stale pattern: {pattern}"

    return False, None
```

### Merged Branch Detection

```bash
# Find merged branches still present
git branch --merged main | grep -v "^\*" | grep -v "main"

# Find branches with no recent commits
git for-each-ref --sort=-committerdate refs/heads/ \
  --format='%(refname:short) %(committerdate:relative)' | \
  grep -E "(months|year)"
```

### Orphan File Detection

Low confidence - requires import analysis:

1. Scan all import statements
2. Build dependency graph
3. Find files with no incoming edges
4. Exclude entry points (main.py, __init__.py, etc.)

## Category 3: Organizational Mess (Low Priority)

### Action: Rename, move, or flag for review

### Detection Patterns

| Issue | Detection Method | Example | Confidence |
|-------|------------------|---------|------------|
| **Duplicates** | Content hash match | Two identical config files | HIGH |
| **Near-duplicates** | >90% content similarity | `config.py` vs `config_old.py` | MEDIUM |
| **Misplaced files** | Extension in wrong dir | `.py` in `docs/` | MEDIUM |
| **Naming inconsistency** | Mixed conventions | `MyFile.py` + `other_file.py` | HIGH |
| **Temp files** | Pattern match | `*.tmp`, `*~`, `*.swp` | HIGH |
| **Backup files** | Pattern match | `*.bak`, `*.backup`, `*~` | HIGH |
| **Version suffixes** | Pattern match | `file_v2.py`, `file_final.py` | MEDIUM |
| **Conflicting names** | Case difference | `Config.py` vs `config.py` | HIGH |

### Duplicate Detection

```python
import hashlib

def find_duplicates(directory):
    hash_map = {}
    duplicates = []

    for file_path in directory.rglob("*"):
        if file_path.is_file():
            file_hash = hashlib.md5(file_path.read_bytes()).hexdigest()
            if file_hash in hash_map:
                duplicates.append((file_path, hash_map[file_hash]))
            else:
                hash_map[file_hash] = file_path

    return duplicates
```

### Naming Inconsistency Detection

```python
from collections import Counter

def detect_naming_inconsistency(directory):
    patterns = {
        'snake_case': r'^[a-z][a-z0-9_]*$',
        'kebab-case': r'^[a-z][a-z0-9-]*$',
        'camelCase': r'^[a-z][a-zA-Z0-9]*$',
        'PascalCase': r'^[A-Z][a-zA-Z0-9]*$',
    }

    file_patterns = Counter()
    for file_path in directory.glob("*.py"):
        name = file_path.stem
        for pattern_name, regex in patterns.items():
            if re.match(regex, name):
                file_patterns[pattern_name] += 1
                break

    # Flag if no dominant pattern (>80%)
    total = sum(file_patterns.values())
    dominant = file_patterns.most_common(1)[0]
    if dominant[1] / total < 0.8:
        return True, file_patterns
    return False, None
```

### Temp File Patterns

| Pattern | Typical Source | Action |
|---------|----------------|--------|
| `*.tmp` | Various | Delete |
| `*~` | Vim, Emacs | Delete |
| `*.swp`, `*.swo` | Vim | Delete |
| `#*#` | Emacs | Delete |
| `.#*` | Emacs | Delete |
| `*.bak` | Various | Review/Delete |
| `*.orig` | Merge conflicts | Delete |
| `*.rej` | Patch failures | Delete |

## Severity Scoring

### Base Scores

| Category | Base Score | Rationale |
|----------|------------|-----------|
| Generated Files | Critical (10) | Should never be tracked |
| Stale Content | Moderate (5) | May contain valuable work |
| Organizational Mess | Low (2) | Inconvenience, not critical |

### Modifiers

| Condition | Modifier | Applied To |
|-----------|----------|------------|
| Size > 100MB | +2 | Generated |
| Tracked in git | +1 | Generated |
| >180 days old | +3 | Stale |
| Merged branch | +2 | Stale |
| In root directory | +1 | Organizational |
| Multiple inconsistencies | +1 each | Organizational |

### Total Score Calculation

```python
def calculate_clutter_score(items):
    total = 0
    for item in items:
        base = BASE_SCORES[item.category]
        modifiers = sum(m for m in item.modifiers if m.applies)
        total += base + modifiers
    return total
```

### Score Interpretation

| Score | Severity | Recommendation |
|-------|----------|----------------|
| 0-10 | Clean | No action needed |
| 11-30 | Minor | Optional cleanup |
| 31-50 | Moderate | Cleanup recommended |
| 51-75 | Significant | Cleanup strongly recommended |
| 76+ | Critical | Immediate cleanup needed |

## Detection Order

1. **Generated Files First** - Quick wins, high impact
2. **Temp Files** - Easy cleanup
3. **Duplicates** - Clear action
4. **Stale Content** - Requires review
5. **Organizational Issues** - Part of larger reorganization

## Edge Cases

### False Positives to Avoid

| Pattern | Why It's Not Clutter |
|---------|---------------------|
| `__init__.py` | Required for packages |
| `__main__.py` | Entry point convention |
| `.gitkeep` | Intentional placeholder |
| `requirements.txt` | Dependency tracking |
| `Makefile` | Build automation |
| `.env.example` | Template for secrets |

### Context-Dependent

| Pattern | Clutter If | Not Clutter If |
|---------|------------|----------------|
| `data/` | Large binaries tracked | Properly gitignored |
| `*.csv` | In root, tracked | In `data/`, gitignored |
| `old-*` | No recent commits | Active development |

## Reporting Format

```markdown
## Clutter Summary

- **Total items**: 47
- **Total score**: 68 (Significant)
- **Critical**: 15 items (score: 42)
- **Moderate**: 12 items (score: 18)
- **Low**: 20 items (score: 8)

## Top Priority Items

| Path | Category | Score | Recommendation |
|------|----------|-------|----------------|
| node_modules/ | Generated | 12 | Add to .gitignore |
| __pycache__/ | Generated | 10 | Add to .gitignore |
| old-experiment/ | Stale | 8 | Archive or delete |
| config.py.bak | Mess | 3 | Delete |
```

## References

- [GitHub gitignore templates](https://github.com/github/gitignore)
- [Git best practices](https://sethrobertson.github.io/GitBestPractices/)
