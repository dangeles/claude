# Structure Template: Code Projects

## Overview

This reference defines directory structure templates for code projects. Templates are based on community standards including Cookiecutter, GitHub conventions, and language-specific best practices.

## Python Project

### Standard Structure

```
project/
├── src/
│   └── {package_name}/
│       ├── __init__.py
│       ├── main.py           # Entry point
│       ├── config.py         # Configuration
│       ├── models/           # Data models
│       │   └── __init__.py
│       ├── services/         # Business logic
│       │   └── __init__.py
│       └── utils/            # Utilities
│           └── __init__.py
├── tests/
│   ├── __init__.py
│   ├── conftest.py           # Shared fixtures
│   ├── unit/
│   │   └── test_*.py
│   └── integration/
│       └── test_*.py
├── docs/
│   ├── api/                  # API documentation
│   ├── guides/               # User guides
│   └── conf.py               # Sphinx config
├── scripts/
│   ├── setup.sh
│   └── deploy.sh
├── .github/
│   └── workflows/
│       └── ci.yml
├── pyproject.toml            # Project config, dependencies
├── README.md
├── CLAUDE.md
├── CHANGELOG.md
├── LICENSE
└── .gitignore
```

### Minimal Structure (Small Projects)

```
project/
├── src/
│   └── {package_name}/
│       ├── __init__.py
│       └── main.py
├── tests/
│   └── test_main.py
├── pyproject.toml
├── README.md
├── CLAUDE.md
└── .gitignore
```

### CLI Application

```
project/
├── src/
│   └── {package_name}/
│       ├── __init__.py
│       ├── cli.py            # CLI entry point
│       ├── commands/         # Subcommands
│       │   ├── __init__.py
│       │   └── run.py
│       └── core/             # Core logic
│           └── __init__.py
├── tests/
├── pyproject.toml
├── README.md
├── CLAUDE.md
└── .gitignore
```

## JavaScript/TypeScript Project

### Standard Structure (Node.js)

```
project/
├── src/
│   ├── index.ts              # Entry point
│   ├── config/
│   │   └── index.ts
│   ├── controllers/          # Request handlers
│   │   └── index.ts
│   ├── services/             # Business logic
│   │   └── index.ts
│   ├── models/               # Data models
│   │   └── index.ts
│   ├── middleware/           # Express middleware
│   │   └── index.ts
│   └── utils/
│       └── index.ts
├── tests/
│   ├── unit/
│   │   └── *.test.ts
│   └── integration/
│       └── *.test.ts
├── dist/                     # Compiled output (gitignored)
├── docs/
├── scripts/
├── .github/
│   └── workflows/
│       └── ci.yml
├── package.json
├── tsconfig.json
├── README.md
├── CLAUDE.md
└── .gitignore
```

### React Application

```
project/
├── public/
│   ├── index.html
│   └── favicon.ico
├── src/
│   ├── index.tsx             # Entry point
│   ├── App.tsx               # Root component
│   ├── components/           # React components
│   │   ├── common/           # Shared components
│   │   └── features/         # Feature-specific
│   ├── hooks/                # Custom hooks
│   ├── context/              # React context
│   ├── services/             # API services
│   ├── utils/
│   ├── types/                # TypeScript types
│   └── styles/               # Global styles
├── tests/
│   └── __tests__/
├── package.json
├── tsconfig.json
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Rust Project

### Library

```
project/
├── src/
│   ├── lib.rs                # Library entry
│   ├── error.rs              # Error types
│   └── utils.rs
├── tests/                    # Integration tests
│   └── integration_test.rs
├── benches/                  # Benchmarks
│   └── benchmark.rs
├── examples/                 # Example usage
│   └── basic.rs
├── Cargo.toml
├── README.md
├── CLAUDE.md
└── .gitignore
```

### Binary Application

```
project/
├── src/
│   ├── main.rs               # Binary entry
│   ├── lib.rs                # Library code
│   ├── cli.rs                # CLI parsing
│   └── commands/
│       └── mod.rs
├── tests/
├── Cargo.toml
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Go Project

### Standard Structure

```
project/
├── cmd/
│   └── {app_name}/
│       └── main.go           # Entry point
├── internal/                 # Private packages
│   ├── config/
│   │   └── config.go
│   ├── handlers/
│   │   └── handlers.go
│   └── models/
│       └── models.go
├── pkg/                      # Public packages
│   └── utils/
│       └── utils.go
├── api/                      # API definitions
│   └── openapi.yaml
├── scripts/
├── docs/
├── go.mod
├── go.sum
├── Makefile
├── README.md
├── CLAUDE.md
└── .gitignore
```

## Directory Purposes

### Core Directories

| Directory | Purpose | Contents |
|-----------|---------|----------|
| `src/` | Source code | Main application code |
| `tests/` | Test code | Unit, integration, e2e tests |
| `docs/` | Documentation | Guides, API docs, ADRs |
| `scripts/` | Utility scripts | Build, deploy, setup scripts |

### Supporting Directories

| Directory | Purpose | Contents |
|-----------|---------|----------|
| `examples/` | Example usage | Working examples |
| `benches/` | Benchmarks | Performance tests |
| `tools/` | Development tools | Code generators, etc. |
| `configs/` | Configuration | Environment configs |

### Generated Directories (gitignored)

| Directory | Purpose | Language |
|-----------|---------|----------|
| `dist/` | Compiled output | TypeScript/JS |
| `build/` | Build artifacts | General |
| `target/` | Compiled output | Rust |
| `__pycache__/` | Bytecode cache | Python |
| `node_modules/` | Dependencies | Node.js |
| `.venv/` | Virtual environment | Python |

## Required Files

### Universal

| File | Purpose | Required |
|------|---------|----------|
| `README.md` | Project overview | Always |
| `CLAUDE.md` | AI assistant context | Always |
| `.gitignore` | Git exclusions | Always |
| `LICENSE` | Legal | For open source |

### Language-Specific

| File | Purpose | Language |
|------|---------|----------|
| `pyproject.toml` | Package config | Python |
| `package.json` | Package config | Node.js |
| `Cargo.toml` | Package config | Rust |
| `go.mod` | Module config | Go |

## Migration Guidelines

### From Flat to Structured

```
# Before (flat)
project/
├── main.py
├── utils.py
├── test_main.py
└── README.md

# After (structured)
project/
├── src/
│   └── project/
│       ├── __init__.py
│       ├── main.py
│       └── utils.py
├── tests/
│   └── test_main.py
├── pyproject.toml
├── README.md
├── CLAUDE.md
└── .gitignore
```

### Migration Steps

1. Create target directory structure
2. Move source files to `src/{package}/`
3. Move test files to `tests/`
4. Update imports in all files
5. Create/update package config
6. Update any CI/CD paths
7. Test that everything still works

## Anti-Patterns

### Avoid These Structures

| Pattern | Problem | Solution |
|---------|---------|----------|
| Code in root | Hard to package | Use `src/` |
| Tests with source | Confused imports | Separate `tests/` |
| Deeply nested | Hard to navigate | Max 4 levels |
| God directories | >50 files | Split by feature/type |
| Mixed concerns | Code + data + docs | Separate directories |

## Scaling Patterns

### Small -> Medium

Add:
- `tests/unit/` and `tests/integration/`
- `docs/api/` and `docs/guides/`
- `scripts/` for automation

### Medium -> Large

Add:
- Feature-based organization in `src/`
- `CODEOWNERS` for ownership
- `docs/adr/` for architecture decisions
- `.github/` for CI/CD and templates

### Large -> Monorepo

Consider:
- `packages/` or `apps/` directory
- Shared `libs/` or `common/`
- Root `package.json` with workspaces
- Lerna, Nx, or Turborepo

## References

- [Cookiecutter](https://github.com/cookiecutter/cookiecutter)
- [GitHub Repository Structure](https://docs.github.com/en/repositories)
- [Python Packaging Guide](https://packaging.python.org/)
- [The Clean Architecture](https://blog.cleancoder.com/uncle-bob/2012/08/13/the-clean-architecture.html)
