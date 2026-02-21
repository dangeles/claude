# Agent Loading Guide

Maps each agent to the specific reference file **sections** it should load. This
prevents agents from loading entire 500-1300 line reference files into context,
which would exhaust the context window when multiple agents are invoked per phase.

**Principle**: Load the minimum context needed. Never load an entire 1000+ line file.
Use the Read tool with line offset and limit parameters to load only the sections listed below.

**Base path**: All reference files are in `cfd-bioreactor/references/` relative to the
skills directory. Agents receive the absolute path from the orchestrator in their
Task tool invocation.

---

## 1. cfd-mathematician Loading Map

The mathematician needs governing equations, variational forms, dimensionless
analysis, and analytical solutions for verification.

| Reference File | Section | Lines | Content | Purpose |
|---|---|---|---|---|
| physics-models.md | 1. Incompressible Navier-Stokes Equations | 12-86 | Strong form, weak form, function spaces | Flow variational formulation |
| physics-models.md | 2. Species Transport Equations | 87-172 | ADR equation, SUPG weak form | Transport variational formulation |
| physics-models.md | 3. O2 Sink Models | 173-303 | Michaelis-Menten, regularization | Reaction term analysis |
| physics-models.md | 6. Non-Dimensionalization Guide | 506-552 | Re, Pe, Da definitions and regimes | Dimensionless number computation |
| validation-benchmarks.md | Benchmark 1: Poiseuille Flow | 8-199 | Analytical solution, error norms | Verification targets |

**Total estimated lines**: ~450
**Estimated tokens**: ~3,000

**Loading instruction for orchestrator**: When invoking cfd-mathematician, include in the Task prompt:
```
Load these reference sections using the Read tool:
- physics-models.md lines 12-303 (Sections 1-3: equations and models)
- physics-models.md lines 506-552 (Section 6: dimensionless numbers)
- validation-benchmarks.md lines 8-199 (Benchmark 1: Poiseuille analytical solution)
```

---

## 2. cfd-reviewer Loading Map

The reviewer needs failure modes, mesh quality criteria, expected validation
results, and physical parameter ranges to challenge plans effectively.

| Reference File | Section | Lines | Content | Purpose |
|---|---|---|---|---|
| troubleshooting-guide.md | Stage 4: Solver Execution Errors | 269-437 | Divergence, NaN, slow convergence | Known failure pattern recognition |
| troubleshooting-guide.md | Quick Diagnostic Checklist | 534-589 | Stage-by-stage checks | Systematic review protocol |
| mesh-generation-guide.md | 5. Memory Estimation | 427-518 | DOF-to-RAM formulas, limits | Memory feasibility checks |
| mesh-generation-guide.md | 6. Mesh Quality Metrics | 519-592 | Jacobian, aspect ratio, angle | Mesh quality threshold validation |
| validation-benchmarks.md | Benchmark 1: Poiseuille Flow | 8-199 | Expected error magnitudes | Flow validation targets |
| validation-benchmarks.md | Mass and Species Conservation Check Protocol | 464-539 | Conservation formulas, tolerances | Post-solve validation criteria |
| physics-models.md | 5. Bioreactor Parameter Lookup Tables | 402-505 | Parameter ranges, physical limits | Physical plausibility checks |

**Total estimated lines**: ~600
**Estimated tokens**: ~4,000

**Loading instruction for orchestrator**: When invoking cfd-reviewer, include in the Task prompt:
```
Load these reference sections using the Read tool:
- troubleshooting-guide.md lines 269-437 (Stage 4: solver errors)
- troubleshooting-guide.md lines 534-589 (Quick Diagnostic Checklist)
- mesh-generation-guide.md lines 427-592 (Sections 5-6: memory + quality)
- validation-benchmarks.md lines 8-199 (Benchmark 1: expected results)
- validation-benchmarks.md lines 464-539 (Conservation check protocol)
- physics-models.md lines 402-505 (Section 5: parameter tables)
```

---

## 3. Orchestrator Loading Map for Code Generation

The orchestrator loads reference sections when generating code (Steps 2.4, 3.4).
These are loaded by the orchestrator itself, not passed to agents.

### Phase 1: Mesh Code Generation

| Reference File | Section | Lines | Content |
|---|---|---|---|
| mesh-generation-guide.md | 1. gmsh Python API Fundamentals | 11-59 | Kernel selection, init/finalize |
| mesh-generation-guide.md | 2. STEP File Import | 60-177 | OCC import pattern |
| mesh-generation-guide.md | 3. Physical Group Definition | 178-291 | Tag conventions |
| mesh-generation-guide.md | 4. Mesh Refinement Strategies | 292-426 | Graded, boundary layer |
| fenicsx-patterns.md | 1. Version Assertion Pattern | 16-51 | Required script header |
| fenicsx-patterns.md | 2. Import Patterns | 52-102 | Canonical imports |
| fenicsx-patterns.md | 15. Reproducibility Header | 1252-1325 | Metadata header |

### Phase 2: Flow Code Generation

| Reference File | Section | Lines | Content |
|---|---|---|---|
| fenicsx-patterns.md | 3. Function Space Definition | 103-169 | Taylor-Hood setup |
| fenicsx-patterns.md | 4. Boundary Condition Application | 170-260 | Dirichlet/Neumann/Robin |
| fenicsx-patterns.md | 5. Stokes Variational Form | 261-307 | Stokes weak form code |
| fenicsx-patterns.md | 6. Navier-Stokes Variational Form | 308-404 | NS with Newton |
| fenicsx-patterns.md | 8. Solver Configuration | 509-613 | MUMPS, GMRES, Newton params |
| fenicsx-patterns.md | 9. Solver Progress Monitor | 614-705 | Convergence logging |
| physics-models.md | 1. Incompressible Navier-Stokes Equations | 12-86 | Weak form reference |

### Phase 3: Transport Code Generation

| Reference File | Section | Lines | Content |
|---|---|---|---|
| fenicsx-patterns.md | 7. SUPG-Stabilized Transport Form | 405-508 | SUPG weak form code |
| fenicsx-patterns.md | 10. Post-Solve Quality Checks | 706-810 | Conservation checks |
| fenicsx-patterns.md | 11. Solution Output | 811-874 | XDMF/VTK patterns |
| physics-models.md | 2. Species Transport Equations | 87-172 | ADR equation reference |
| physics-models.md | 3. O2 Sink Models | 173-303 | MM reaction code |
| physics-models.md | 4. Boundary Condition Patterns | 304-401 | Robin BC for membranes |

---

## 4. Knowledge Redundancy Principle

Agent SKILL.md files reference domain concepts by name (e.g., "check SUPG stability",
"verify inf-sup condition", "assess Michaelis-Menten regularization"). They describe
**what** to analyze, not **how** -- the formulas, numerical constants, and parameter
values live exclusively in reference files.

**Rules**:
1. Agent SKILL.md files MUST NOT contain formulas, numerical constants, or parameter values.
2. All numerical and formula knowledge lives exclusively in the reference files listed above.
3. Agents load the relevant reference sections at invocation time via the Read tool.
4. **Exception**: The orchestrator's Section 13 (Code Generation Protocol) embeds safety-critical
   "ALWAYS" and "NEVER" rules inline. These rules must always be in the orchestrator's context
   because they guard against code-level errors that could produce silently wrong results.

---

## 5. Section Dependency Table

When updating a reference file section, consult this table to identify which agents
are affected and may need their loading instructions verified.

| Reference File | Section | Depended On By |
|---|---|---|
| physics-models.md | 1. Navier-Stokes | cfd-mathematician, orchestrator (Phase 2) |
| physics-models.md | 2. Species Transport | cfd-mathematician, orchestrator (Phase 3) |
| physics-models.md | 3. O2 Sink Models | cfd-mathematician, orchestrator (Phase 3) |
| physics-models.md | 4. Boundary Conditions | orchestrator (Phase 3) |
| physics-models.md | 5. Parameter Tables | cfd-reviewer |
| physics-models.md | 6. Non-Dimensionalization | cfd-mathematician |
| physics-models.md | 7. Coupling Strategy | orchestrator (Phase 2-3 sequencing) |
| validation-benchmarks.md | Benchmark 1 | cfd-mathematician, cfd-reviewer |
| validation-benchmarks.md | Conservation Protocol | cfd-reviewer |
| troubleshooting-guide.md | Stage 4: Solver Errors | cfd-reviewer |
| troubleshooting-guide.md | Quick Diagnostic | cfd-reviewer |
| mesh-generation-guide.md | 1-4. gmsh Patterns | orchestrator (Phase 1) |
| mesh-generation-guide.md | 5. Memory Estimation | cfd-reviewer |
| mesh-generation-guide.md | 6. Quality Metrics | cfd-reviewer |
| fenicsx-patterns.md | 1-2. Version + Imports | orchestrator (all phases) |
| fenicsx-patterns.md | 3-6. Function Spaces, BCs, Stokes, NS | orchestrator (Phase 2) |
| fenicsx-patterns.md | 7, 10. SUPG, Quality Checks | orchestrator (Phase 3) |
| fenicsx-patterns.md | 8-9. Solver Config, Monitor | orchestrator (Phase 2) |
| fenicsx-patterns.md | 11, 15. Output, Reproducibility | orchestrator (Phase 2-3) |
| environment-setup.md | 4. Pre-Flight Validation | orchestrator (Phase 0) |
