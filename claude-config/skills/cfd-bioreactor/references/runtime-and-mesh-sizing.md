# Runtime and Mesh Sizing

Operational defaults for simulation mesh sizing, expected wall-clock times, execution
strategy, and solver convergence diagnostics. Extracted from SKILL.md to keep the
orchestration layer short.

## Contents

1. Debug vs. production mesh sizing
2. Expected wall-clock times and kill thresholds
3. Execution strategy by problem type
4. Solver convergence assessment
5. Mesh convergence study decision guide

---

## 1. Debug vs. Production Mesh Sizing

Mesh sizing defaults are tied to the complexity tier and `run_purpose`. The orchestrator
determines `run_purpose` during Phase 0. Two levels exist:

- `debug`: first run, error recovery, parameter exploration. Coarsest mesh that captures
  qualitative physics.
- `production`: publication-quality results. Fine mesh with convergence study.

**Tier-aware mesh sizing defaults**:

| Tier | Run Purpose | Elements (2D) | Elements (3D) | Rationale |
|------|-------------|---------------|----------------|-----------|
| 1 | Always debug | 200-500 | N/A (2D only) | Pipeline validation; exact solution known |
| 2 | Always debug | 500-2,000 | N/A (2D only) | Qualitative physics check; fast iteration |
| 3 | Debug (first run) | N/A | 5,000-20,000 | Coarse 3D captures flow structure |
| 3 | Production | N/A | 50,000-200,000 | Quantitative accuracy for publication |
| 4 | Debug (first run) | N/A | 10,000-50,000 | Validate before committing hours |
| 4 | Production | N/A | 200,000-1,000,000+ | Publication quality with convergence study |

**Geometric feature resolution rule**: regardless of run_purpose, keep at least 3
elements across the thinnest geometric feature (membranes, thin channels, boundary
layers). If the default element count range for a tier and run_purpose cannot achieve
this, increase the element count to the minimum that satisfies the constraint and warn
the user: "Debug mesh increased from [default] to [adjusted] elements to resolve
[feature] (>= 3 elements across thinnest feature required)."

**How the orchestrator uses this**:

1. Set `run_purpose` in Phase 0 from tier and session history.
2. During mesh planning, pass `run_purpose` and the corresponding element count range to
   `cfd-mathematician` as a constraint.
3. During code generation, use these defaults unless the mathematician or reviewer
   recommends otherwise.
4. When `run_purpose == debug`, skip the Phase 4 mesh convergence study.

**Override**: the user can change run purpose at any time. State the consequence:

```
Current run purpose: debug (Tier 3, first run)
Override to production? This will increase mesh from ~10K to ~100K elements
and runtime from ~2 minutes to ~30-60 minutes.
```

---

## 2. Expected Wall-Clock Times and Kill Thresholds

Serial execution, single core. Times cover simulation execution only, not agent
discussion.

| Tier | Run Purpose | Mesh + Flow | Transport | Total Expected | Kill Threshold (3x) |
|------|-------------|-------------|-----------|----------------|---------------------|
| 1 | debug | 30s | N/A | 30s | 90s |
| 2 | debug | 1-3 min | 1-2 min | 2-5 min | 15 min |
| 3 | debug | 2-5 min | 2-5 min | 5-10 min | 30 min |
| 3 | production | 15-30 min | 10-20 min | 30-60 min | 3 hours |
| 4 | debug | 5-15 min | 5-10 min | 10-25 min | 75 min |
| 4 | production | 1-4 hours | 30 min-2 hours | 2-6 hours | 18 hours |

**Kill threshold**: wall-clock time beyond 3x the expected total for the tier and run
purpose indicates a stalled or diverging simulation.

When a run is executed via Bash, the Bash timeout kills it. For scripts the user runs
manually (3D production), embed a timeout check in the generated script:

```python
import time
_start_time = time.time()
_kill_threshold_s = {kill_threshold_seconds}

# Inside solver loop:
if time.time() - _start_time > _kill_threshold_s:
    print(f"TIMEOUT: Simulation exceeded {_kill_threshold_s}s kill threshold.")
    print("Likely causes: mesh too fine for debug, solver diverging, or")
    print("insufficient Newton continuation steps.")
    print("Try: reduce mesh size, increase continuation ramp steps,")
    print("or switch to Picard iteration.")
    sys.exit(1)
```

Report the timeout cause in the Phase 5 summary with specific recommendations.

**Stall detection heuristic** (solver progress monitoring): if the solver has not reduced
the residual by more than 10x within the last 5 minutes of wall-clock time (production
runs) or 30 seconds (debug runs), treat it as stalled. A stall means no progress at all,
which is different from slow convergence.

---

## 3. Execution Strategy by Problem Type

| Problem Type | Expected Runtime | Execution Strategy |
|---|---|---|
| 2D validation (Tier 1) | < 2 minutes | Bash (timeout=300000) |
| 2D coupled (Tier 2) | 2-10 minutes | Bash (timeout=600000) |
| 3D coarse validation | 5-15 minutes | Bash (timeout=600000) if small; Write+manual otherwise |
| 3D production (Tier 3) | 30 min - 4 hours | Write script only. User runs manually. |
| 3D with convergence (Tier 4) | 4-24 hours | Write script only. User runs with MPI. |
| Mesh convergence study | 4x single-run time | Write script only for 3D. |

When estimated runtime exceeds 10 minutes, generate the script and instruct the user to
run it manually rather than executing it via Bash. For 3D simulations, include MPI
instructions:

```bash
mpirun -np 4 python simulation.py
```

---

## 4. Solver Convergence Assessment

Monitor these during execution.

**Newton solver convergence**:

| Metric | Healthy | Warning | Diverging |
|--------|---------|---------|-----------|
| Residual norm after 5 iterations | < 1e-4 | 1e-4 to 1e-1 | > 1e-1 or increasing |
| Residual reduction per iteration | > 10x | 2x-10x | < 2x or oscillating |
| Total iterations to converge | < 10 | 10-20 | > 20 or not converging |

**Krylov solver convergence** (GMRES/CG inner solves):

| Metric | Healthy | Warning | Diverging |
|--------|---------|---------|-----------|
| Krylov iterations per Newton step | < 100 | 100-500 | > 500 |
| Krylov residual stagnation | N/A | Stalls for > 50 iterations | Residual increases |

**Action on convergence signals**:

- **Healthy**: continue.
- **Warning**: log the warning, continue, flag it in the Phase 5 summary. When
  `run_purpose == debug`, note "convergence marginal on debug mesh -- may improve on
  finer production mesh".
- **Diverging**: enter the self-correction loop. On a first attempt with a debug mesh, do
  not refine the mesh as the first fix; try solver parameter adjustments instead
  (relaxation, Picard, continuation ramp).

---

## 5. Mesh Convergence Study Decision Guide

The Phase 4 mesh convergence study is optional. Run it when:

- `run_purpose == production` and tier >= 3
- the user asks for publication-quality results
- QoI sensitivity is high: a 10% mesh coarsening changes the QoI by more than 5%

Skip it when:

- `run_purpose == debug`
- tier 1-2 (an analytical solution is available for validation instead)
- the user declines ("I just need qualitative results")

**Convergence assessment criteria**:

- QoI change between the last two refinement levels < 1%
- observed convergence rate within 0.5 of the theoretical order:
  - P1 concentration: expect p ~ 2.0, accept 1.5-2.5
  - P2 velocity: expect p ~ 3.0, accept 2.5-3.5
- an anomalously low rate (p < 1.0) suggests a singularity or insufficient boundary layer
  resolution
