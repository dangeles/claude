# Naming Conventions: Code Projects

## Overview

This reference defines naming conventions for code projects (Python, JavaScript/TypeScript, Rust, Go, etc.). Conventions are based on community standards and industry best practices.

## Python Projects

### Files

| Element | Convention | Example | Notes |
|---------|------------|---------|-------|
| Modules | snake_case | `my_module.py` | Lowercase with underscores |
| Packages | snake_case | `my_package/` | Same as modules |
| Tests | `test_*.py` | `test_utils.py` | pytest convention |
| Scripts | snake_case | `run_migrations.py` | Descriptive verbs |

### Directories

| Element | Convention | Example |
|---------|------------|---------|
| Package dirs | snake_case | `my_package/` |
| Test dirs | `tests/` | `tests/unit/`, `tests/integration/` |
| Source root | `src/` | `src/my_package/` |
| Docs | `docs/` | `docs/api/` |

### Special Files

| File | Convention | Required |
|------|------------|----------|
| `__init__.py` | Exact name | For packages |
| `__main__.py` | Exact name | For runnable packages |
| `conftest.py` | Exact name | pytest fixtures |
| `setup.py` | Exact name | Legacy packaging |
| `pyproject.toml` | Exact name | Modern packaging |

### Python-Specific Rules

1. **Never start with numbers**: `2024_analysis.py` is invalid for imports
2. **Avoid hyphens**: Use `my_module.py`, not `my-module.py`
3. **Private modules**: Prefix with underscore `_internal.py`

## JavaScript/TypeScript Projects

### Files

| Element | Convention | Example | Notes |
|---------|------------|---------|-------|
| Components | PascalCase | `UserProfile.tsx` | React convention |
| Utilities | camelCase or kebab-case | `apiClient.ts`, `api-client.ts` | Project-dependent |
| Tests | `*.test.ts` or `*.spec.ts` | `utils.test.ts` | Jest/Vitest convention |
| Configs | lowercase | `webpack.config.js` | Tool convention |

### Directories

| Element | Convention | Example |
|---------|------------|---------|
| Components | PascalCase or kebab-case | `UserProfile/`, `user-profile/` |
| General dirs | kebab-case | `src/`, `lib/`, `utils/` |
| Test dirs | `__tests__/` or `tests/` | Jest convention |

### Special Files

| File | Convention | Required |
|------|------------|----------|
| `index.ts` | Exact name | Barrel exports |
| `package.json` | Exact name | Required |
| `.eslintrc.js` | Dotfile | Linting config |
| `tsconfig.json` | Exact name | TypeScript config |

## Rust Projects

### Files

| Element | Convention | Example |
|---------|------------|---------|
| Modules | snake_case | `my_module.rs` |
| Tests | `tests/*.rs` | `tests/integration_test.rs` |
| Binaries | snake_case | `my_cli.rs` |

### Directories

| Element | Convention | Example |
|---------|------------|---------|
| Source | `src/` | `src/lib.rs`, `src/main.rs` |
| Tests | `tests/` | Integration tests |
| Benches | `benches/` | Benchmarks |

### Special Files

| File | Convention |
|------|------------|
| `lib.rs` | Library entry |
| `main.rs` | Binary entry |
| `mod.rs` | Module directory |
| `Cargo.toml` | Package manifest |

## Go Projects

### Files

| Element | Convention | Example |
|---------|------------|---------|
| Source | snake_case | `user_service.go` |
| Tests | `*_test.go` | `user_service_test.go` |
| Generated | `*_gen.go` | `models_gen.go` |

### Directories

| Element | Convention | Example |
|---------|------------|---------|
| Commands | `cmd/` | `cmd/myapp/` |
| Internal | `internal/` | Private packages |
| Packages | `pkg/` | Public packages |

## Universal Conventions

### Configuration Files

| File | Convention |
|------|------------|
| `.gitignore` | Exact, dotfile |
| `.env` | Exact, dotfile |
| `.env.example` | Template |
| `README.md` | Exact, uppercase |
| `CHANGELOG.md` | Exact, uppercase |
| `LICENSE` | Exact, uppercase |
| `CLAUDE.md` | Exact, uppercase |

### Documentation

| Element | Convention | Example |
|---------|------------|---------|
| Docs dir | `docs/` | Always lowercase |
| API docs | `docs/api/` | Subdirectory |
| Guides | `docs/guides/` | Subdirectory |
| ADRs | `docs/adr/` | Architecture Decision Records |

### Scripts

| Element | Convention | Example |
|---------|------------|---------|
| Scripts dir | `scripts/` | Utility scripts |
| Build scripts | Descriptive | `scripts/build.sh` |
| CI scripts | `.github/` or `.gitlab/` | Platform-specific |

## Anti-Patterns to Avoid

### Critical (Must Fix)

| Pattern | Problem | Solution |
|---------|---------|----------|
| Spaces in names | Breaks shell commands | Use underscores or hyphens |
| Special characters | Encoding issues | Alphanumeric only |
| Starting with number | Breaks imports | Prefix with word |
| Mixed case on same level | Confusion | Pick one convention |

### Warning (Should Fix)

| Pattern | Problem | Solution |
|---------|---------|----------|
| Version in filename | Maintenance burden | Use git for versioning |
| Abbreviations | Unclear meaning | Spell out (or document) |
| Generic names | No context | Be descriptive |
| Very long names | Hard to read | Aim for <50 chars |

### Examples of Good vs Bad

| Bad | Good | Reason |
|-----|------|--------|
| `MyModule.py` | `my_module.py` | Python uses snake_case |
| `file (1).js` | `file_backup.js` | No spaces or special chars |
| `2024analysis.py` | `analysis_2024.py` | Don't start with number |
| `data-v2.csv` | `data_20240115.csv` | Date, not version |
| `cfg.py` | `config.py` | Spell it out |
| `stuff.py` | `utilities.py` | Be descriptive |

## Adaptive Mode

When analyzing existing projects:

1. **Sample 20+ files** to detect dominant pattern
2. **Calculate consistency**: % following dominant pattern
3. **If >80% consistent**: Enforce existing pattern
4. **If <80% consistent**: Recommend project-type default
5. **Document detected pattern** in report

## References

- [PEP 8](https://peps.python.org/pep-0008/) - Python Style Guide
- [Airbnb JavaScript Style Guide](https://github.com/airbnb/javascript)
- [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- [Effective Go](https://go.dev/doc/effective_go)
