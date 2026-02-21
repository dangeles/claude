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

| Reference File | Section | Section Header | Approx Lines | Content | Purpose |
|---|---|---|---|---|---|
| physics-models.md | 1 | `## 1. Incompressible Navier-Stokes Equations` | ~75 | Strong form, weak form, function spaces | Flow variational formulation |
| physics-models.md | 2 | `## 2. Species Transport Equations` | ~85 | ADR equation, SUPG weak form | Transport variational formulation |
| physics-models.md | 3 | `## 3. O2 Sink Models` | ~130 | Michaelis-Menten, regularization | Reaction term analysis |
| physics-models.md | 6 | `## 6. Non-Dimensionalization Guide` | ~45 | Re, Pe, Da definitions and regimes | Dimensionless number computation |
| validation-benchmarks.md | Benchmark 1 | `## Benchmark 1: Poiseuille Flow` | ~190 | Analytical solution, error norms | Verification targets |

**Total estimated lines**: ~525
**Estimated tokens**: ~3,500

**Loading instruction for orchestrator**: When invoking cfd-mathematician, include in the Task prompt:
```
Load these reference sections using the Read tool (use absolute paths from {skill_base_path}/references/):
- physics-models.md: from "## 1. Incompressible Navier-Stokes Equations" through end of "## 3. O2 Sink Models" (~290 lines)
- physics-models.md: section "## 6. Non-Dimensionalization Guide" through next ## heading (~45 lines)
- validation-benchmarks.md: section "## Benchmark 1: Poiseuille Flow" through next ## heading (~190 lines)

To find section boundaries: Grep for the section header, note the line number.
Grep for the next ## heading to find the end. Read from start to end-1.
```

---

## 2. cfd-reviewer Loading Map

The reviewer needs failure modes, mesh quality criteria, expected validation
results, and physical parameter ranges to challenge plans effectively.

| Reference File | Section | Section Header | Approx Lines | Content | Purpose |
|---|---|---|---|---|---|
| troubleshooting-guide.md | Stage 4 | `## Stage 4: Solver Execution Errors` | ~170 | Divergence, NaN, slow convergence | Known failure pattern recognition |
| troubleshooting-guide.md | Quick Diagnostic | `## Quick Diagnostic Checklist` | ~55 | Stage-by-stage checks | Systematic review protocol |
| mesh-generation-guide.md | 5 | `## 5. Memory Estimation` | ~90 | DOF-to-RAM formulas, limits | Memory feasibility checks |
| mesh-generation-guide.md | 6 | `## 6. Mesh Quality Metrics` | ~75 | Jacobian, aspect ratio, angle | Mesh quality threshold validation |
| validation-benchmarks.md | Benchmark 1 | `## Benchmark 1: Poiseuille Flow` | ~190 | Expected error magnitudes | Flow validation targets |
| validation-benchmarks.md | Conservation | `## Mass and Species Conservation Check Protocol` | ~75 | Conservation formulas, tolerances | Post-solve validation criteria |
| physics-models.md | 5 | `## 5. Bioreactor Parameter Lookup Tables` | ~105 | Parameter ranges, physical limits | Physical plausibility checks |

**Total estimated lines**: ~760
**Estimated tokens**: ~5,000

**Loading instruction for orchestrator**: When invoking cfd-reviewer, include in the Task prompt:
```
Load these reference sections using the Read tool (use absolute paths from {skill_base_path}/references/):
- troubleshooting-guide.md: section "## Stage 4: Solver Execution Errors" through next ## heading (~170 lines)
- troubleshooting-guide.md: section "## Quick Diagnostic Checklist" through end of file (~55 lines)
- mesh-generation-guide.md: from "## 5. Memory Estimation" through end of "## 6. Mesh Quality Metrics" (~165 lines)
- validation-benchmarks.md: section "## Benchmark 1: Poiseuille Flow" through next ## heading (~190 lines)
- validation-benchmarks.md: section "## Mass and Species Conservation Check Protocol" through next ## heading (~75 lines)
- physics-models.md: section "## 5. Bioreactor Parameter Lookup Tables" through next ## heading (~105 lines)

To find section boundaries: Grep for the section header, note the line number.
Grep for the next ## heading to find the end. Read from start to end-1.
```

---

## 3. Orchestrator Loading Map for Code Generation

The orchestrator loads reference sections when generating code (Steps 2.4, 3.4).
These are loaded by the orchestrator itself, not passed to agents.

### Phase 1: Mesh Code Generation

| Reference File | Section | Section Header | Approx Lines | Content |
|---|---|---|---|---|
| mesh-generation-guide.md | 1 | `## 1. gmsh Python API Fundamentals` | ~50 | Kernel selection, init/finalize |
| mesh-generation-guide.md | 2 | `## 2. STEP File Import` | ~115 | OCC import pattern |
| mesh-generation-guide.md | 3 | `## 3. Physical Group Definition` | ~115 | Tag conventions |
| mesh-generation-guide.md | 4 | `## 4. Mesh Refinement Strategies` | ~135 | Graded, boundary layer |
| fenicsx-patterns.md | 1 | `## 1. Version Assertion Pattern` | ~35 | Required script header |
| fenicsx-patterns.md | 2 | `## 2. Import Patterns` | ~50 | Canonical imports |
| fenicsx-patterns.md | 15 | `## 15. Reproducibility Header` | ~75 | Metadata header |

### Phase 2: Flow Code Generation

| Reference File | Section | Section Header | Approx Lines | Content |
|---|---|---|---|---|
| fenicsx-patterns.md | 3 | `## 3. Function Space Definition` | ~65 | Taylor-Hood setup |
| fenicsx-patterns.md | 4 | `## 4. Boundary Condition Application` | ~90 | Dirichlet/Neumann/Robin |
| fenicsx-patterns.md | 5 | `## 5. Stokes Variational Form` | ~45 | Stokes weak form code |
| fenicsx-patterns.md | 6 | `## 6. Navier-Stokes Variational Form` | ~95 | NS with Newton |
| fenicsx-patterns.md | 8 | `## 8. Solver Configuration` | ~105 | MUMPS, GMRES, Newton params |
| fenicsx-patterns.md | 9 | `## 9. Solver Progress Monitor` | ~90 | Convergence logging |
| physics-models.md | 1 | `## 1. Incompressible Navier-Stokes Equations` | ~75 | Weak form reference |

### Phase 3: Transport Code Generation

| Reference File | Section | Section Header | Approx Lines | Content |
|---|---|---|---|---|
| fenicsx-patterns.md | 7 | `## 7. SUPG-Stabilized Transport Form` | ~105 | SUPG weak form code |
| fenicsx-patterns.md | 10 | `## 10. Post-Solve Quality Checks` | ~105 | Conservation checks |
| fenicsx-patterns.md | 11 | `## 11. Solution Output` | ~65 | XDMF/VTK patterns |
| physics-models.md | 2 | `## 2. Species Transport Equations` | ~85 | ADR equation reference |
| physics-models.md | 3 | `## 3. O2 Sink Models` | ~130 | MM reaction code |
| physics-models.md | 4 | `## 4. Boundary Condition Patterns` | ~100 | Robin BC for membranes |

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

### 4b. Authority and Known Discrepancies

**This guide is the authoritative source** for reference file section mappings.
Agent SKILL.md files may reference sections by approximate or different names.
When there is a discrepancy between an agent's internal section reference and
this guide, use the mappings in this guide.

**Known discrepancies with agent SKILL.md files** (agent files are out of scope
for editing; these are documented here for awareness):

| Agent | Agent's Reference | Correct Reference (this guide) |
|-------|-------------------|-------------------------------|
| cfd-reviewer | `troubleshooting-guide.md "Known Failure Modes" section` | `"## Stage 4: Solver Execution Errors"` |
| cfd-reviewer | `mesh-generation-guide.md Section 5 (quality criteria) + Section 6 (memory estimation)` | Section 5 is Memory Estimation; Section 6 is Mesh Quality Metrics (descriptions are swapped in agent) |
| cfd-mathematician | `physics-models.md Sections 1-3 (governing equations, variational forms, dimensionless numbers)` | Section 3 is "O2 Sink Models" (not dimensionless numbers). Dimensionless numbers are in Section 6. Agent omits Section 6 from its table. |

The orchestrator should use THIS guide's section headers when constructing
agent Task prompts, not the agent's own approximate section names.

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
