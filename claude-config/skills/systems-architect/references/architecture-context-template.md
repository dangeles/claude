# Architecture Context Document - Template and Protocols

**Version**: 1.0
**Last Updated**: 2026-02-07
**Purpose**: Template and protocols for `.architecture/context.md` - persistent architectural reference for incremental development

---

## What This Document Is

The Architecture Context Document is an **authoritative reference for architectural intent**. It describes what the architecture SHOULD BE, while the code describes what it IS. When the two diverge, the code is ground truth for implementation decisions, and the divergence is reported for future context updates.

**This is NOT**:
- A substitute for reading the code
- A guarantee of accuracy (may lag behind code changes)
- A complete specification (uses `[TBD]`, `[UNKNOWN]`, `[INFERRED]` where information is incomplete)

---

## Template Structure

### Quick Reference Index

```markdown
# Architecture Context

**Template-Version**: 1.0
**Last-Updated**: {ISO8601 timestamp}
**Codebase-State**: {git commit SHA, optional}

## Quick Reference

| Module | Tier | Key Dependencies | Modification Risk |
|--------|------|------------------|-------------------|
| {module} | {0/1/2} | {comma-separated} | {Low/Medium/High} |

**Tier Definitions**:
- **Tier 0 (Foundation)**: Core utilities, data structures, no dependencies on other modules
- **Tier 1 (Core Business Logic)**: Depends on Tier 0, provides main functionality
- **Tier 2 (Application/Interface)**: Depends on Tier 0/1, user-facing or integration code
```

### Section 1: Module Interconnections and Dependencies

```markdown
## 1. Module Interconnections and Dependencies

**Dependency Graph** (MUST be a DAG - directed acyclic graph):

```
Tier 2: app/cli.py
         ↓ (imports)
Tier 1: src/processor.py, src/analyzer.py
         ↓
Tier 0: src/utils.py, src/data_structures.py
```

### Module Descriptions

#### Tier 0: Foundation Modules

**`src/utils.py`**:
- Purpose: Common utilities (file I/O, logging, configuration)
- Exports: `load_config()`, `setup_logging()`, `read_csv_safe()`
- Dependencies: None (foundation tier)
- Modification order: First (safest - no dependents within project)

**`src/data_structures.py`**:
- Purpose: Core data classes (Sample, Dataset, Result)
- Exports: `Sample`, `Dataset`, `Result` dataclasses
- Dependencies: None
- Modification order: First

#### Tier 1: Core Business Logic

**`src/processor.py`**:
- Purpose: Data processing pipeline
- Exports: `process_dataset()`, `validate_samples()`
- Dependencies: `utils`, `data_structures`
- Modification order: Second (after Tier 0)

**`src/analyzer.py`**:
- Purpose: Statistical analysis and reporting
- Exports: `analyze()`, `generate_report()`
- Dependencies: `utils`, `data_structures`, `processor`
- Modification order: Third (depends on processor)

#### Tier 2: Application Layer

**`app/cli.py`**:
- Purpose: Command-line interface
- Exports: `main()`, CLI argument parsing
- Dependencies: All Tier 0 and Tier 1 modules
- Modification order: Last (highest-risk - depends on everything)

**Informational Dependencies** (not import dependencies):
- `src/analyzer.py` expects `processor.py` to produce `Dataset` objects with `.samples` attribute
- `app/cli.py` expects `analyzer.py` to write reports to `output/` directory
```

### Section 2: Intended Usage Patterns

```markdown
## 2. Intended Usage Patterns for Each Module

### `src/utils.py`
**Intended**: Utility functions called by any module needing file I/O or logging
**NOT Intended**: Business logic, data processing, or state management
**Example**: `config = load_config('config.yaml')` in any module

### `src/processor.py`
**Intended**: Single entry point `process_dataset(dataset)` for all processing
**NOT Intended**: Direct manipulation of `Dataset` internals outside `processor.py`
**Example**: `processed = process_dataset(raw_dataset)` from CLI or tests

### `src/analyzer.py`
**Intended**: Analyze processed datasets only (assumes `process_dataset()` has run)
**NOT Intended**: Processing raw data or bypassing validation
**Example**: `report = analyze(processed_dataset)` after processing

### `app/cli.py`
**Intended**: User-facing entry point, orchestrates calls to processor → analyzer
**NOT Intended**: Direct file manipulation or business logic
**Example**: `python app/cli.py --input data.csv --output results/`
```

### Section 3: Modification Order Considerations

```markdown
## 3. Modification Order to Avoid Breaking Changes

**General Rule**: Modify in topological order (Tier 0 → Tier 1 → Tier 2).

### Scenarios and Safe Orders

**Scenario 1: Adding a new utility function**
1. Add function to `src/utils.py` (Tier 0)
2. Use function in any Tier 1 or Tier 2 module
**Risk**: Low (new code, no existing dependents)

**Scenario 2: Changing `Dataset` dataclass schema**
1. Update `src/data_structures.py` (Tier 0)
2. Update `src/processor.py` to use new schema (Tier 1)
3. Update `src/analyzer.py` if it accesses changed fields (Tier 1)
4. Update tests
5. Update `app/cli.py` if user-facing output format changes (Tier 2)
**Risk**: High (Tier 0 change affecting all layers)

**Scenario 3: Adding a new analysis function**
1. Add function to `src/analyzer.py` (Tier 1)
2. Call from `app/cli.py` if exposing to users (Tier 2)
**Risk**: Medium (Tier 1 change, check all callers)

**Breaking Change Protocol**:
- If changing a module signature, check `git grep "function_name"` for all call sites
- Update call sites in dependent modules before committing
- Run full test suite before committing
```

### Section 4: Streaming and Incremental Change Strategies

```markdown
## 4. Streaming/Incremental Change Strategies

### Strategy 1: Feature Flags for Gradual Rollout
**When**: Adding experimental features to `analyzer.py`
**How**:
- Add `--enable-new-analysis` CLI flag
- Wrap new code in `if config.enable_new_analysis:` blocks
- Default to `False` until validated
- Remove flag after stabilization

### Strategy 2: Parallel Implementations During Refactoring
**When**: Rewriting `processor.py` with new algorithm
**How**:
- Keep `process_dataset()` (old implementation) functional
- Add `process_dataset_v2()` alongside
- Switch CLI to use v2 via `--use-v2-processor` flag
- Deprecate v1 after migration period
- Delete v1 after all users migrated

### Strategy 3: Backwards-Compatible Schema Evolution
**When**: Adding fields to `Dataset` dataclass
**How**:
- Add new fields as optional (default values)
- Old code continues working (ignores new fields)
- New code uses new fields when available
- Gradual migration: Tier 0 → Tier 1 → Tier 2

### Strategy 4: Isolated Module Addition
**When**: Adding new Tier 1 module `src/validator.py`
**How**:
- Add module without modifying existing code
- Import and use only from Tier 2 (`app/cli.py`)
- Run in parallel with existing flow initially
- Integrate fully after validation
```

---

## Generation Protocols

### Phase 3 (Standard Architecture Design)

During programming-pm Phase 3, systems-architect generates `.architecture/context.md`:

1. **Read architecture handoff** from programming-pm (contains module list, dependencies, design decisions)
2. **Populate all four sections** using architecture design
3. **Assign tiers** (0/1/2) based on dependency structure
4. **Verify DAG constraint**: Ensure no circular dependencies (cycles = architectural problem)
5. **Write to `.architecture/context.md`** in project root
6. **Create directory** `.architecture/` if it does not exist

### Bootstrap Mode (Existing Codebases)

For projects without prior Architecture Context Document:

1. **Scan directory structure**: List all modules in `src/`, `modules/`, `app/`
2. **Infer dependencies** from import statements (use `git grep "^import" "^from"`)
3. **Mark uncertainties**:
   - `[TBD]`: Placeholder for information that should be added later
   - `[UNKNOWN]`: Information that couldn't be determined from static analysis
   - `[INFERRED]`: Information guessed from code structure (may be wrong)
4. **Populate Quick Reference Index** with inferred tiers
5. **Document "Known Gaps" section** at end listing all `[TBD]`, `[UNKNOWN]`, `[INFERRED]` items

**Bootstrap Quality Standard**: Incomplete but honest documentation is better than fabricated completeness.

### SIMPLE Mode

For SIMPLE workflows (small changes, single module), systems-architect generates abbreviated context:

- Only Section 1 (Module Interconnections) required
- Sections 2-4 can be `[TBD - Not needed for SIMPLE workflow]`
- Quick Reference Index required

---

## Maintenance Protocols

### When TO Update

Update `.architecture/context.md` when:
1. **New module added** to `src/` or `modules/`
2. **Module deleted or renamed**
3. **Dependency changed** (new import statement between modules)
4. **Interface changed** (function signature, dataclass schema)
5. **Developer reports discrepancy** via code handoff field `architecture_context.discrepancy_noted: true`

### When NOT to Update

Do NOT update for:
1. **Bug fixes** within a module (no interface changes)
2. **Internal refactoring** (same inputs/outputs)
3. **Documentation changes**
4. **Test additions** (unless test reveals architectural issue)

### Staleness Detection

Programming-pm Phase 5 checks for staleness:

```bash
# Get last context update timestamp
CONTEXT_UPDATED=$(git log -1 --format=%aI -- .architecture/context.md)

# Get last structural change timestamp (narrowed to reduce false positives)
STRUCTURAL_CHANGE=$(git log -1 --format=%aI -- "src/*/__init__.py" "modules/*/__init__.py")

# If STRUCTURAL_CHANGE > CONTEXT_UPDATED, context may be stale
```

**Action on Stale Detection**: Log informational note in progress file. Do NOT warn prominently (prevents warning fatigue).

### Drift Handling (Phase 5)

If developer handoff includes `architecture_context.discrepancy_noted: true`:
1. **Systems-architect updates context** with targeted file modification (<10 min)
2. **Update specific sections** that diverged (not full regeneration)
3. **Log reason for update** in context document's revision history

---

## Merge Conflict Resolution

If `.architecture/context.md` has merge conflicts:

**Protocol**: Accept base (older version), then re-run systems-architect in bootstrap mode.

**Rationale**: Document can be regenerated from codebase. Better to have slightly outdated context than manually merged nonsense.

**Steps**:
```bash
git checkout --theirs .architecture/context.md  # Accept incoming version
# OR
git checkout --ours .architecture/context.md    # Accept current version

# Then regenerate from code
# (systems-architect will run in Phase 3 on next workflow invocation)
```

---

## Special Markers

Use these markers for incomplete information:

- `[TBD]`: To Be Determined - needs human input
- `[UNKNOWN]`: Cannot be determined from available information
- `[INFERRED]`: Guessed from code structure, may be incorrect

**Developer Instructions**: When you encounter these markers, treat the code as ground truth for your implementation decisions. If you discover the correct information, note it in your code handoff for context update.

---

## Template-Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-07 | Initial template with 4 sections, bootstrap mode, SIMPLE mode support |
