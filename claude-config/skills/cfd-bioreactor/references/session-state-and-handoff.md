# Session State and Project Handoff

Session state file schema, resume behavior, and the criteria for handing a session off to
`programming-pm`. Extracted from SKILL.md to keep the orchestration layer short.

## Contents

1. Session state file schema
2. Session resume protocol
3. programming-pm handoff criteria

---

## 1. Session State File Schema

State file location: `/tmp/cfd-bioreactor-state-{session-id}.yaml`

```yaml
session:
  session_id: "20260221-143000-12345"
  session_dir: "/tmp/cfd-bioreactor-session-20260221-143000-12345"
  skill_base_path: "/Users/davidangelesalbores/repos/claude/claude-config/skills/cfd-bioreactor"
  mode: "LITE"                    # DIRECT | LITE | FULL
  tier: 2                         # 1-4
  run_purpose: "debug"            # debug | production
  handoff_declined: false         # Set to true if user declines programming-pm handoff
  current_phase: 2                # 0-5
  status: "running"               # running | paused | complete | failed

phases_completed: [0, 1]          # List of completed phase numbers
phases_valid: [0, 1]              # List of phases whose results are still valid

phase_results:
  phase_0:
    status: "complete"
    handoff_path: null
    timestamp: "2026-02-21T14:30:00Z"
  phase_1:
    status: "complete"
    handoff_path: "{session_dir}/handoffs/phase1-mesh-plan.yaml"
    reviewer_status: "APPROVED_WITH_WARNINGS"
    retries_used: 0
    timestamp: "2026-02-21T14:35:00Z"
  phase_2:
    status: "running"
    handoff_path: null
    reviewer_status: null
    retries_used: 0
    execution_retries_used: 0     # For self-correction loop

timestamps:
  started: "2026-02-21T14:30:00Z"
  last_updated: "2026-02-21T14:40:00Z"

errors: []                        # List of non-fatal errors/warnings encountered
```

---

## 2. Session Resume Protocol

On startup, look for `/tmp/cfd-bioreactor-state-*.yaml` with `status != "complete"`. If
one exists, offer to resume:

```
Found incomplete session from [timestamp].
Last completed phase: Phase [N] ([phase_name])
Mode: [DIRECT/LITE/FULL], Tier: [N]
Resume from Phase [N+1]? (yes/no/restart)
```

On resume:

- Phase 0 re-runs regardless (it is fast and catches environment changes).
- Phases 1-4 are resumable when their handoff files still exist at the recorded paths.
  Check for the files before resuming from them.

On restart, create a new session and archive the old state file.

---

## 3. programming-pm Handoff Criteria

The cfd-bioreactor orchestrator generates standalone simulation scripts. When the user's
needs grow beyond one-off scripts into maintained software, recommend handoff to
`programming-pm`.

### Handoff Trigger Conditions

| Trigger | Detection Heuristic | Example User Request |
|---------|---------------------|---------------------|
| Library packaging | User mentions "package", "library", "module", "import from", "reusable", "API" | "Make this simulation into a library I can call with different parameters" |
| Test/CI requirements | User mentions "pytest", "CI", "continuous integration", "coverage", "unit test" | "Add unit tests for the simulation pipeline" |
| Version control integration | User mentions "git repo", "branch", "PR", "version", "release" | "Set up version control for the simulation code" |
| Multi-file project growth | Session has generated >= 4 scripts AND user requests shared utilities | "Extract the mesh generation into a shared module" |
| Parameter sweep framework | User wants automated parameter sweeps with result aggregation | "Run this simulation for 20 different inlet velocities and collect results" |

### Context Filtering Rules

To avoid false-positive handoff triggers:

1. **Phase filtering**: evaluate handoff triggers only at phase boundaries (Phase 0,
   Phase 5) or during user interaction between phases, not during Phases 2-4 (active
   computation).
2. **Compound requirement**: require keywords from 2+ different trigger categories (for
   example "library" + "test") or a single unambiguous keyword with no alternate CFD
   interpretation ("package this as a library", "set up CI").
3. **Exclusion patterns**: these are not handoff triggers:
   - "test" in "test this simulation", "test run", "test different parameters" (means
     "try", not "unit test")
   - "git" in "let me check" or "let me get" (not version control)
   - "version" in "FEniCSx version", "Python version" (not software versioning)
   - "module" in "Python module not found" (import error, not packaging)
4. **Decline suppression**: if the user declines the handoff, set `handoff_declined: true`
   in session state and do not re-suggest programming-pm for the rest of the session.

### Handoff Recommendation Format

```
HANDOFF RECOMMENDATION: Software Development Coordination

Your simulation work is evolving into a software project:
- Trigger: [which trigger was detected]
- Current state: [N scripts generated, current tier/mode]

Recommendation: Hand off to programming-pm for:
- Architecture design (module structure, shared utilities)
- Testing strategy (unit tests for simulation components)
- Code review and quality gates
- Version control integration

Options:
(A) Hand off to programming-pm now (recommended)
(B) Continue with cfd-bioreactor (standalone scripts only)
(C) Complete current simulation phase, then hand off
```

### What Gets Handed Off

- All generated scripts from `{session_dir}/scripts/`
- Session state summary (tier, mode, run_purpose, phases completed)
- Physical parameters and validation results
- List of scripts and their purposes

### When Not to Hand Off

- Standard tier progression (Tier 1 -> 2 -> 3 -> 4)
- Mesh convergence studies (these are cfd-bioreactor Phase 4)
- One-off parameter variations (just re-run with different inputs)
- Visualization enhancements (cfd-bioreactor Phase 4)
- Any session where `handoff_declined` is already `true`
