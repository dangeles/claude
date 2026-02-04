# Gitignore Patterns by Project Type

## Overview

This reference provides gitignore patterns for different project types. Patterns are sourced from GitHub's official gitignore repository and community best practices.

## Python Projects

### Essential Patterns

```gitignore
# Byte-compiled / optimized / DLL files
__pycache__/
*.py[cod]
*$py.class

# C extensions
*.so

# Distribution / packaging
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
share/python-wheels/
*.egg-info/
.installed.cfg
*.egg
MANIFEST

# Virtual environments
.env
.venv
env/
venv/
ENV/
env.bak/
venv.bak/

# Installer logs
pip-log.txt
pip-delete-this-directory.txt

# Unit test / coverage
htmlcov/
.tox/
.nox/
.coverage
.coverage.*
.cache
nosetests.xml
coverage.xml
*.cover
*.py,cover
.hypothesis/
.pytest_cache/
cover/

# Type checkers
.mypy_cache/
.dmypy.json
dmypy.json
.pyre/
.pytype/

# Jupyter Notebook
.ipynb_checkpoints

# pyenv
.python-version

# Environments
.env
.env.local
*.env

# IDE
.idea/
.vscode/
*.swp
*.swo
```

### Data Science Additions

```gitignore
# Data files (if large)
*.csv
*.parquet
*.h5
*.hdf5
*.pkl
*.joblib
data/raw/*
data/processed/*
!data/raw/.gitkeep
!data/processed/.gitkeep

# Model files
models/*
!models/.gitkeep
*.onnx

# MLflow
mlruns/
mlartifacts/

# Weights & Biases
wandb/

# DVC
.dvc/cache
.dvc/tmp
```

## Node.js / JavaScript Projects

### Essential Patterns

```gitignore
# Dependencies
node_modules/
jspm_packages/

# Build outputs
dist/
build/
out/
.next/
.nuxt/
.cache/

# Logs
logs
*.log
npm-debug.log*
yarn-debug.log*
yarn-error.log*
lerna-debug.log*
.pnpm-debug.log*

# Runtime data
pids
*.pid
*.seed
*.pid.lock

# Coverage
coverage/
*.lcov
.nyc_output

# TypeScript
*.tsbuildinfo

# Optional npm cache directory
.npm

# Optional eslint cache
.eslintcache

# Parcel cache
.parcel-cache

# Yarn
.yarn/cache
.yarn/unplugged
.yarn/build-state.yml
.yarn/install-state.gz
.pnp.*

# Environment variables
.env
.env.local
.env.development.local
.env.test.local
.env.production.local

# IDE
.idea/
.vscode/
*.swp
*.swo
```

### React / Next.js Additions

```gitignore
# Next.js
.next/
out/
.vercel

# Storybook
storybook-static/

# Testing
jest-results/
playwright-report/
test-results/
```

## Rust Projects

### Essential Patterns

```gitignore
# Build artifacts
/target/

# Cargo.lock for libraries (keep for binaries)
# Cargo.lock

# IDE
.idea/
.vscode/
*.swp

# MacOS
.DS_Store

# Profiling
perf.data
perf.data.old
flamegraph.svg
```

## Go Projects

### Essential Patterns

```gitignore
# Binaries
*.exe
*.exe~
*.dll
*.so
*.dylib

# Test binary
*.test

# Output of go coverage
*.out

# Go workspace
go.work

# Dependency directories
vendor/

# IDE
.idea/
.vscode/

# Build output
/bin/
/build/
```

## General / Cross-Platform

### Operating System

```gitignore
# macOS
.DS_Store
.AppleDouble
.LSOverride
._*
.Spotlight-V100
.Trashes

# Windows
Thumbs.db
Thumbs.db:encryptable
ehthumbs.db
ehthumbs_vista.db
*.stackdump
[Dd]esktop.ini

# Linux
*~
.fuse_hidden*
.directory
.Trash-*
.nfs*
```

### IDE / Editor

```gitignore
# JetBrains (IntelliJ, PyCharm, etc.)
.idea/
*.iml
*.ipr
*.iws
out/

# VS Code
.vscode/*
!.vscode/settings.json
!.vscode/tasks.json
!.vscode/launch.json
!.vscode/extensions.json
*.code-workspace

# Vim
*.swp
*.swo
*~
.netrwhist

# Emacs
*~
\#*\#
/.emacs.desktop
/.emacs.desktop.lock
*.elc
auto-save-list
tramp

# Sublime Text
*.sublime-workspace
*.sublime-project
```

### Secrets and Credentials

```gitignore
# Environment files
.env
.env.local
.env.*.local
*.env

# Credentials
credentials.json
credentials.yaml
secrets.json
secrets.yaml
*.pem
*.key
*.crt

# AWS
.aws/credentials
aws-credentials.json

# Google Cloud
*.json
!package.json
!tsconfig.json
service-account*.json
```

### Logs and Temporary Files

```gitignore
# Logs
*.log
logs/

# Temporary files
*.tmp
*.temp
*.bak
*.backup
*~
.tmp/
tmp/
temp/

# Cache
.cache/
cache/
```

## Data Projects

### Essential Patterns

```gitignore
# Raw data (large files)
data/raw/*
!data/raw/.gitkeep
!data/raw/README.md

# Processed data
data/processed/*
!data/processed/.gitkeep

# External data
data/external/*
!data/external/.gitkeep

# Features
data/features/*
!data/features/.gitkeep

# Models
models/*
!models/.gitkeep

# Database files
*.db
*.sqlite
*.sqlite3

# Parquet, CSV, etc.
*.parquet
*.csv
*.tsv
*.xlsx
*.xls

# Pickle
*.pkl
*.pickle

# HDF5
*.h5
*.hdf5

# Feather
*.feather

# Apache Arrow
*.arrow

# Checkpoint files
checkpoints/*
!checkpoints/.gitkeep
```

### dbt Specific

```gitignore
# dbt
target/
dbt_packages/
dbt_modules/
logs/
```

### Airflow Specific

```gitignore
# Airflow
airflow.cfg
airflow.db
airflow-webserver.pid
logs/
```

## Combining Patterns

### Template for Multi-Language Project

```gitignore
# ============================================
# Operating System
# ============================================
.DS_Store
Thumbs.db

# ============================================
# IDE / Editor
# ============================================
.idea/
.vscode/
*.swp
*.swo
*~

# ============================================
# Python
# ============================================
__pycache__/
*.py[cod]
.venv/
.pytest_cache/
.mypy_cache/
*.egg-info/
dist/
build/

# ============================================
# Node.js
# ============================================
node_modules/
dist/
.next/

# ============================================
# Secrets
# ============================================
.env
.env.local
*.pem
credentials.json

# ============================================
# Logs and Temp
# ============================================
*.log
*.tmp
logs/

# ============================================
# Data (if applicable)
# ============================================
data/raw/*
data/processed/*
!data/*/.gitkeep
```

## Usage Guidelines

### When to Gitignore

| Category | Rule | Example |
|----------|------|---------|
| Generated | Always ignore | `__pycache__/`, `node_modules/` |
| Secrets | Always ignore | `.env`, `credentials.json` |
| Large data | Usually ignore | `*.parquet`, `*.csv` |
| IDE settings | Usually ignore | `.idea/`, `.vscode/` |
| OS files | Always ignore | `.DS_Store`, `Thumbs.db` |
| Logs | Usually ignore | `*.log`, `logs/` |

### When NOT to Gitignore

| File | Reason |
|------|--------|
| `Cargo.lock` (binaries) | Reproducible builds |
| `package-lock.json` | Reproducible builds |
| `poetry.lock` | Reproducible builds |
| `.gitkeep` | Preserve empty directories |
| Sample data (<1MB) | Documentation, testing |

### Adding .gitkeep

To track empty directories:
```bash
touch data/raw/.gitkeep
touch data/processed/.gitkeep
touch models/.gitkeep
```

## References

- [GitHub gitignore repository](https://github.com/github/gitignore)
- [gitignore.io](https://www.toptal.com/developers/gitignore)
