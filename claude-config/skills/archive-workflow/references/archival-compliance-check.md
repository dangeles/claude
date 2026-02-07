# Archival Compliance Check Procedure

## Purpose
This document defines the reusable archival compliance check that all file-generating
workflows follow before writing output files. It is referenced by all consuming
workflows via a short integration section in their SKILL.md.

## Source of Truth Hierarchy

When determining archival guidelines, the following precedence applies (highest first):
1. `.archive-metadata.yaml` (if present and valid)
2. CLAUDE.md (fallback, only for lit-pm and programming-pm)
3. Workflow defaults (baseline, no archival enforcement)

## Prerequisites

- YAML parser: `yaml.safe_load()` (Python) or `yq` (bash). Never use `yaml.load()`.
- If neither is available, log WARNING and skip archival checks.

## Step 1: Detect Archival Metadata

```bash
REPO_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
ARCHIVE_METADATA="${REPO_ROOT}/.archive-metadata.yaml"

if [ -f "$ARCHIVE_METADATA" ]; then
  ARCHIVAL_GUIDELINES_PRESENT=true
else
  ARCHIVAL_GUIDELINES_PRESENT=false
  # No archival guidelines -- proceed with workflow defaults
fi
```

### Exception: archive-workflow Internal Invocation
If this workflow was invoked via handoff with `archival_context: "skip"`:
- Do NOT check for `.archive-metadata.yaml`
- Proceed with the invoking workflow's instructions
- This prevents circular enforcement when archive-workflow invokes specialists

### Exception: Orchestrator-Provided Context
If the handoff from an orchestrator includes an `archival_context` block:
- Use the provided context directly
- Do NOT re-read `.archive-metadata.yaml`
- The orchestrator has already resolved guidelines and user choices

## Step 2: Read and Validate

Parse using safe methods only:

```python
import yaml, os

with open(archive_metadata_path, 'r', encoding='utf-8') as f:
    data = yaml.safe_load(f)

if data is None:
    warn("Archival metadata file is empty. Treating as no guidelines.")
    return NO_GUIDELINES
```

### Schema Validation

After parsing, validate required fields:

| Field Path | Type | Required | Valid Values | If Missing/Invalid |
|------------|------|----------|-------------|-------------------|
| version | string | YES | Pattern "N.N" | Treat entire file as corrupt |
| generated.timestamp | string | YES | ISO8601 | Treat as corrupt |
| project.type | string | YES | "code", "research", "data", "mixed" | Default to "mixed" |
| naming_conventions.summary.files | string | YES | "snake_case", "kebab-case", "camelCase", "PascalCase" | Skip naming checks |
| naming_conventions.summary.directories | string | YES | "lowercase", "snake_case", "kebab-case" | Skip directory naming checks |
| naming_conventions.project_specific_rules | list | NO | List of {pattern, convention, example} | Use summary only |
| naming_conventions.anti_patterns | list | NO | List of strings | Skip anti-pattern checks |
| naming_conventions.full_reference | string | NO | Absolute file path | Use summary only |
| structure.top_level_directories | list | YES (min 1) | List of {name, purpose, enforced} | Skip directory checks |
| structure.full_reference | string | NO | Absolute file path | Use summary only |
| enforcement.mode | string | NO | "advisory", "soft-mandatory", "hard-mandatory" | Default to "advisory" |

If any required field fails validation:
- Log WARNING with the specific field and its invalid value
- Treat that section as "no guidelines" (partial application)
- Continue with valid sections

### Version Compatibility

Check the `version` field against supported version "1.0":
- If major version matches: proceed
- If major version is higher: WARN "Schema version {version} is newer than supported"
- If version is missing: treat as corrupt

### Path Validation

All paths read from the YAML must be validated:
1. Reject any path containing ".."
2. `full_reference` paths must resolve within `~/.claude/skills/archive-workflow/`
3. Apply `os.path.expanduser()` defensively on all paths
4. If a reference file does not exist, log WARNING and use summary fields only

## Step 3: Validate Proposed File Paths

For each file the workflow intends to create, check:

### Directory Check
1. File in a LISTED directory with `enforced: true` -- always validate naming
2. File in a LISTED directory with `enforced: false` -- validate with softer warning
3. File in an UNLISTED directory that ALREADY EXISTS in the repo -- no violation
4. File in an UNLISTED directory that DOES NOT EXIST -- potential violation

Only case 4 triggers the advisory.

### Naming Check
Match the file extension against `naming_conventions.project_specific_rules`:
- If a matching rule exists, verify the filename follows the convention
- If no matching rule exists, check against `naming_conventions.summary.files`

### Anti-Pattern Check
Check the filename against `naming_conventions.anti_patterns` list.

## Step 4: Advisory Resolution

### Enforcement Mode Behavior
Read `enforcement.mode` from the YAML (default: "advisory"):
- `advisory`: Present options A/B/C, never block
- `soft-mandatory`: Present options, default to archival, log override
- `hard-mandatory`: Only options A and C available (no workflow default override)

### Batch Advisory (Primary Pattern)

Before creating files in a workflow stage, enumerate ALL files, validate ALL paths,
and collect violations into a batch. Present ONE summary:

```
ARCHIVAL GUIDELINE ADVISORY ({N} violations in {M} files)
============================================================

All violations share the root cause:
  {Description of the common issue}

Affected files:
  1. {proposed_path_1} -> {archival_path_1}
  2. {proposed_path_2} -> {archival_path_2}
  [... more]

Options:
  (A) Apply archival guidelines to ALL
  (B) Keep workflow defaults for ALL
  (C) Review individually
  (D) Custom directory: [specify]

Recommendation: (A) Follow archival guidelines
```

### Single-File Advisory (Fallback)

If only one file violates, or if user chose option (C) from batch:

```
ARCHIVAL GUIDELINE ADVISORY
============================

Proposed file: {proposed_path}
Violation: {description}

Archival guideline says: {archival_recommendation}
Workflow default says: {workflow_default}

Options:
  (A) Follow archival guidelines: {archival_path}
  (B) Follow workflow default: {default_path}
  (C) Custom path: [specify]

Recommendation: (A) Follow archival guidelines
```

### Context-Dependent Behavior
- When invoked by orchestrator via Task tool (sub-agent context): Apply archival
  guidelines silently using the orchestrator's archival_context. No interactive prompts.
- When invoked standalone by user: Present advisory options interactively.

## Step 5: Record Decision

Log the user's choice in the workflow session state:

```yaml
archival_decisions:
  - proposed_path: "{path}"
    violation_type: "{type}"
    chosen_option: "A|B|C|D"
    final_path: "{resolved_path}"
    timestamp: "{ISO8601}"
```

## Producer-Consumer Contract

### Producer (archive-workflow) GUARANTEES:
1. `.archive-metadata.yaml` will be syntactically valid YAML (or not exist)
2. All fields marked REQUIRED will be present and non-null
3. All enum fields will contain values from the specified set
4. The `version` field will follow semantic versioning (MAJOR.MINOR)
5. The file will be written atomically (temp file + rename)
6. All string values will be double-quoted
7. All paths in `full_reference` will be absolute (no tildes)
8. The file will be committed to git as part of the archive-workflow commit

### Consumer RESPONSIBILITIES:
1. MUST parse with yaml.safe_load() or yq (NEVER yaml.load())
2. MUST check `version` field for compatibility
3. MUST handle missing `.archive-metadata.yaml` (proceed without checks)
4. MUST handle malformed YAML (treat as missing, log warning)
5. MUST present user with options on violation (never silently enforce)
6. MUST log compliance decisions in session state
7. MUST NOT modify `.archive-metadata.yaml`
8. MUST ignore unknown fields (forward compatibility)

### Contract VIOLATIONS:
| Violation | Handler |
|-----------|---------|
| Invalid YAML syntax | Consumer: treat as missing |
| Missing required field | Consumer: degrade that section |
| Unknown field in YAML | Consumer: ignore it |
| Version too new | Consumer: best-effort + warn |
| File manually edited | Consumer: validate, use if valid |

## Schema Versioning Policy

- MINOR versions (1.0 -> 1.1): Additive changes only. New optional fields.
  Consumers need no changes.
- MAJOR versions (1.0 -> 2.0): Breaking changes. Consumers must be updated.
  Migration window: both schemas valid for 30 days.
- Forward compatibility: Consumers MUST ignore unknown fields.
- Default-to-v1: If `version` field is missing, treat as corrupt (not as v1).

## Standard Archival Context Block (for orchestrator handoffs)

All orchestrators MUST include this block when dispatching specialists:

```yaml
archival_context:
  guidelines_present: true/false
  source: ".archive-metadata.yaml"  # or "CLAUDE.md" or "defaults"
  naming_convention: "snake_case"
  output_directory: "docs/"
  enforcement_mode: "advisory"
  user_override: null  # If user chose non-archival option, record here
```

## Deprecation Notice

### Deprecated: Direct CLAUDE.md Reading for Archival Guidelines
- **Deprecated in**: v1.0 (2026-02-06)
- **Replaced by**: `.archive-metadata.yaml` via this compliance check
- **Fallback**: CLAUDE.md reading preserved ONLY when `.archive-metadata.yaml` absent
- **Removal timeline**: CLAUDE.md fallback will be removed in v2.0 (TBD)

When CLAUDE.md fallback is triggered, log:
"Archival guidelines read from CLAUDE.md (fallback). Run archive-workflow to
generate .archive-metadata.yaml for structured guidelines."
