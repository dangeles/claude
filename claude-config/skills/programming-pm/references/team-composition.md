# Team Composition Guide

Reference guide for selecting team composition in programming-pm orchestrated projects.

---

## Default Team (Always Required)

These skills are always included in every programming-pm project.

| Skill | Role | Phase(s) | Responsibility |
|-------|------|----------|----------------|
| programming-pm | Orchestrator | All | Coordination, quality gates, state management |
| requirements-analyst | Requirements | 1 | Scope definition, success criteria |
| systems-architect | Architecture | 3 | System design, component boundaries |
| senior-developer | Implementation | 4-5 | Production code, integration tests, code review |
| copilot | Review Support | 5 | Code review assistance |

---

## Optional Specialists

Include these specialists based on project characteristics.

### junior-developer

**Include when**:
- Project has >3 independent implementation tasks
- Tasks can be decomposed into well-scoped units
- Routine implementations following established patterns

**Example triggers**:
- "Create multiple utility functions"
- "Implement CRUD operations for several models"
- "Add unit tests for existing code"

**Do NOT include when**:
- All tasks require senior judgment
- Project is small (<1 week effort)
- Tight deadline without time for review cycles

### mathematician

**Include when** (keywords):
- "algorithm", "complexity", "O(n)", "Big-O"
- "optimization", "minimize", "maximize"
- "numerical", "approximation", "precision"
- "sort", "search", "graph", "tree"
- "matrix", "linear algebra", "eigenvalue"

**Example triggers**:
- "Design an efficient search algorithm"
- "Implement numerical integration"
- "Optimize the matching algorithm"

**Do NOT include when**:
- Using standard library algorithms without modification
- Simple operations without complexity concerns
- Statistical analysis (use statistician)

### statistician

**Include when** (keywords):
- "statistics", "statistical", "p-value", "significance"
- "Monte Carlo", "simulation", "sampling"
- "MCMC", "Markov chain", "Bayesian"
- "confidence interval", "uncertainty"
- "bootstrap", "resampling", "permutation"
- "power analysis", "sample size", "effect size"

**Example triggers**:
- "Validate the Monte Carlo simulation"
- "Design hypothesis test for A/B comparison"
- "Implement MCMC sampler with convergence diagnostics"

**Do NOT include when**:
- Simple descriptive statistics (mean, median)
- Non-statistical numerical methods
- Algorithm design (use mathematician)

---

## Team Selection Decision Tree

```
START
  │
  ├── Is this a software development project?
  │     ├── NO → Use technical-pm instead
  │     └── YES → Include default team
  │
  ├── Does project have >3 decomposable tasks?
  │     ├── NO → senior-developer handles all
  │     └── YES → Include junior-developer
  │
  ├── Does project involve algorithm design?
  │     ├── NO → Continue
  │     └── YES → Include mathematician
  │
  ├── Does project involve statistical analysis?
  │     ├── NO → Continue
  │     └── YES → Include statistician
  │
  └── Final team composition determined
```

---

## RACI Matrix

Responsibility assignment for common deliverables.

### Legend
- **R** = Responsible (does the work)
- **A** = Accountable (approves/owns outcome)
- **C** = Consulted (provides input)
- **I** = Informed (kept updated)
- **-** = Not involved

### Project Deliverables

| Deliverable | programming-pm | senior-developer | junior-developer | mathematician | statistician |
|-------------|----------------|------------------|------------------|---------------|--------------|
| Requirements approval | A | C | - | C | C |
| Pre-mortem facilitation | R/A | C | I | C | C |
| Architecture design | C | C | - | C | C |
| Algorithm specification | C | C | - | R/A | C |
| Statistical specification | C | C | - | C | R/A |
| Implementation (complex) | I | R/A | C | C | C |
| Implementation (routine) | I | A | R | - | - |
| Code review (junior) | - | R/A | I | - | - |
| Code review (senior) | I | I | - | C | C |
| Unit tests | I | A | R | - | - |
| Integration tests | A | R | C | C | C |
| PR creation | I | R | - | - | - |
| PR merge decision | A | R | - | - | - |

### Phase Involvement

| Phase | programming-pm | senior-developer | junior-developer | mathematician | statistician |
|-------|----------------|------------------|------------------|---------------|--------------|
| 1. Requirements | A | C | I | C | C |
| 2. Pre-mortem | R/A | C | I | C | C |
| 3. Architecture | C | C | I | C | C |
| 4. Implementation | I | R/A | R | R | R |
| 5. Code Review | A | R | I | C | C |
| 6. VCS Integration | A | R | I | - | - |

---

## User Override Options

Users can override automatic team selection.

### Force Include Specialist

```bash
# Include mathematician even if not detected
programming-pm --include mathematician "Implement sorting algorithm with complexity constraints"

# Include multiple specialists
programming-pm --include mathematician --include statistician "Complex simulation project"
```

### Force Exclude Specialist

```bash
# Exclude auto-detected specialist
programming-pm --exclude statistician "Data pipeline without statistical validation"

# Exclude junior-developer (senior handles all)
programming-pm --exclude junior-developer "Time-critical security fix"
```

### Minimal Team

```bash
# Only PM + senior-developer
programming-pm --minimal "Simple CRUD API"

# Minimal plus one specialist
programming-pm --minimal --include mathematician "Simple algorithm implementation"
```

---

## Team Size Guidelines

### Small Project (1-2 weeks)

- programming-pm
- senior-developer
- (Optional) 1 specialist if domain-specific

### Medium Project (2-4 weeks)

- programming-pm
- senior-developer
- junior-developer (if >3 tasks)
- (Optional) 1-2 specialists as needed

### Large Project (>4 weeks)

- programming-pm
- senior-developer
- junior-developer
- Specialists as needed
- Consider breaking into phases with checkpoints

---

## Specialist Coordination Patterns

### mathematician + senior-developer

```
mathematician → Algorithm Specification → senior-developer
                      │
                      ├── Pseudocode
                      ├── Complexity analysis
                      ├── Verification criteria
                      └── Implementation guidance
```

### statistician + senior-developer

```
statistician → Statistical Specification → senior-developer
                      │
                      ├── Method selection
                      ├── Validation criteria
                      ├── Convergence diagnostics
                      └── Implementation guidance
```

### mathematician + statistician (both needed)

```
mathematician ─┬─→ Algorithm Specification ─┬─→ senior-developer
               │                            │
statistician ──┴─→ Statistical Specification ┘

Note: Each provides separate specification
      senior-developer integrates both
      programming-pm coordinates handoffs
```

### senior-developer + junior-developer

```
senior-developer ──→ Task Decomposition ──→ junior-developer
                            │                      │
                            │                      ↓
                            │              Implementation
                            │                      │
                            └──── Code Review ←────┘
                                      │
                                      ↓
                               Integration
```

---

## Anti-Patterns to Avoid

### Over-staffing

**Problem**: Including specialists not needed for the project
**Impact**: Overhead, delays, unnecessary complexity
**Solution**: Start minimal, add specialists only when needed

### Under-staffing

**Problem**: Missing specialist for domain-specific work
**Impact**: Quality issues, incorrect implementations
**Solution**: Use keyword detection, don't skip pre-flight checks

### Direct Specialist Communication

**Problem**: Specialists communicating without going through PM
**Impact**: Lost context, inconsistent state, coordination failures
**Solution**: All communication flows through programming-pm

### Role Overreach

**Problem**: Specialist working outside their scope
**Impact**: Quality issues, role conflict
**Solution**: Clear "What X does NOT do" sections, escalate violations

---

## Team Composition Examples

### Example 1: Web API Development

**Project**: REST API for user management

**Team**:
- programming-pm (orchestration)
- senior-developer (implementation)

**Rationale**: Standard CRUD operations, no algorithm or statistics

---

### Example 2: Search Algorithm

**Project**: Implement efficient search for large dataset

**Team**:
- programming-pm (orchestration)
- mathematician (algorithm design)
- senior-developer (implementation)

**Rationale**: Complexity analysis needed, but no statistics

---

### Example 3: Monte Carlo Simulation

**Project**: Option pricing simulation with validation

**Team**:
- programming-pm (orchestration)
- mathematician (numerical methods)
- statistician (convergence validation)
- senior-developer (implementation)

**Rationale**: Both algorithm design and statistical validation needed

---

### Example 4: Data Processing Pipeline

**Project**: ETL pipeline with multiple transformations

**Team**:
- programming-pm (orchestration)
- senior-developer (complex transformations)
- junior-developer (routine transformations)

**Rationale**: Many decomposable tasks, no specialized math/stats

---

## Pre-Flight Validation

Before starting a project, programming-pm validates team composition:

```bash
# Check required skills exist
for skill in requirements-analyst systems-architect senior-developer copilot; do
  [ -f ~/.claude/skills/$skill/SKILL.md ] || echo "ABORT: Missing $skill"
done

# Check optional specialists exist (warn if missing)
for skill in mathematician statistician junior-developer; do
  [ -f ~/.claude/skills/$skill/SKILL.md ] || echo "WARN: $skill not available"
done
```

If a selected specialist is not available:
1. Warn user about limitation
2. Document fallback approach
3. Proceed with reduced capability
