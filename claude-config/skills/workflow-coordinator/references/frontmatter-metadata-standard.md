# Frontmatter Metadata Standard for Handoff-Capable Workflows

**Version**: 1.0
**Schema**: Universal Handoff Schema v3.0
**Date**: 2026-02-07

---

## Important: Extension to Anthropic Standard

Claude Code's built-in skill discovery uses **only two frontmatter fields**:

- `name` -- Skill identifier (used for invocation and display)
- `description` -- Natural language description (used for relevance matching)

**All other fields documented here are custom extensions** consumed by `workflow-coordinator` for handoff discovery and coordination. They are ignored by Claude Code's built-in skill loading mechanism and do not affect normal skill invocation.

These custom fields follow YAML syntax within the standard `---` frontmatter delimiters but have no meaning outside the workflow-coordinator ecosystem.

---

## Field Definitions

### Required Fields (Anthropic Standard)

These fields are required for all skills, regardless of handoff capability.

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `name` | string | Skill identifier (kebab-case) | `programming-pm` |
| `description` | string | Natural language description for Claude's skill matching | `Use when coordinating software development...` |

### Required Fields (Handoff Discovery)

These fields are required for a skill to participate in handoff discovery.

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `handoff.accepts_from` | list[string] | Skills that can hand off to this skill. `"*"` = any. | `["*"]` or `["skill-editor"]` |
| `handoff.provides_to` | list[string] | Skills this can hand off to. `"*"` = any. | `["programming-pm", "technical-pm"]` |
| `handoff.schema_version` | string | Handoff schema version supported | `"3.0"` |

### Recommended Fields (Enhanced Discovery)

These fields improve discovery quality but are not strictly required.

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `handoff.schema_type` | string | `"universal"` or `"internal"` | `"universal"` |
| `categories` | list[string] | Primary categories for discovery | `["implementation", "architecture"]` |

### Optional Fields (Additional Metadata)

These fields provide additional context but do not affect discovery.

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `input_requirements` | list[string] | Types of input this skill accepts | `["specification", "requirements"]` |
| `output_types` | list[string] | Types of output this skill produces | `["implementation", "code", "tests"]` |

---

## Validation Rules

1. **name**: Must be a non-empty string. Should use kebab-case.
2. **description**: Must be a non-empty string. Should start with "Use when" for consistency with Claude's skill matching.
3. **handoff.accepts_from**: Must be a list of strings. Each entry is a skill name or `"*"`.
4. **handoff.provides_to**: Must be a list of strings. Each entry is a skill name or `"*"`.
5. **handoff.schema_version**: Must be a string matching `"X.Y"` format (e.g., `"3.0"`).
6. **handoff.schema_type**: Must be `"universal"` or `"internal"`.
7. **categories**: Must be a list of strings from the controlled vocabulary (see below).
8. **input_requirements**: Must be a list of strings.
9. **output_types**: Must be a list of strings.

### YAML Validation

After any frontmatter modification, validate with:

```bash
python3 -c "
import yaml, sys
content = open('path/to/SKILL.md').read()
parts = content.split('---')
if len(parts) >= 3:
    fm = yaml.safe_load(parts[1])
    if fm and 'name' in fm and 'description' in fm:
        print(f'Valid: {fm[\"name\"]}')
    else:
        print('Invalid: missing name or description', file=sys.stderr)
        sys.exit(1)
else:
    print('Invalid: no frontmatter delimiters', file=sys.stderr)
    sys.exit(1)
"
```

---

## Controlled Vocabulary for Categories

Core categories (from existing orchestrator analysis):

| Category | Description | Example Skills |
|----------|-------------|----------------|
| `implementation` | Code writing, software development | programming-pm |
| `architecture` | System design, structural decisions | programming-pm, systems-architect |
| `research` | Information gathering, literature review | technical-pm, research-pipeline |
| `analysis` | Data analysis, problem decomposition | scientific-analysis-architect, perspective-swarm |
| `creative` | Creative writing, content generation | lit-pm, brainstorming-pm |
| `verification` | Testing, validation, review | completion-verifier |
| `coordination` | Multi-workflow orchestration | workflow-coordinator |
| `skill-development` | Skill creation and modification | skill-editor |
| `project-management` | Project planning and tracking | programming-pm, technical-pm |
| `organization` | File organization, archival | archive-workflow |
| `decision-support` | Multi-perspective analysis | perspective-swarm |
| `workflow-creation` | Creating new workflows | skill-editor |
| `integration` | Cross-system integration | workflow-coordinator |
| `workflow-management` | Managing workflow lifecycles | workflow-coordinator |

Custom categories are allowed but should be documented. Prefer existing categories when applicable.

---

## Minimum Viable Frontmatter

The absolute minimum for a skill to participate in handoff discovery:

```yaml
---
name: my-skill
description: Use when doing X that requires Y

handoff:
  accepts_from:
    - "*"
  provides_to:
    - "*"
  schema_version: "3.0"
---
```

This enables:
- Discovery by workflow-coordinator
- Handoff reception from any skill
- Handoff transmission to any skill

---

## Full Frontmatter Example

A complete example with all fields:

```yaml
---
name: programming-pm
description: Use when coordinating software development projects requiring multiple specialists (architect, developers, mathematician, statistician, notebook-writer) with quality gates for archival setup, requirements, architecture, pre-mortem, code review, testing, and version control integration.

# v3.0 universal handoff metadata (see workflow-coordinator/references/frontmatter-metadata-standard.md)
handoff:
  accepts_from:
    - "*"
  provides_to:
    - senior-developer
    - junior-developer
    - systems-architect
    - requirements-analyst
    - mathematician
    - statistician
    - skill-editor
  schema_version: "3.0"
  schema_type: universal
  # Legacy fields preserved for backward compatibility
  handoff_trigger: "--handoff {payload_path}"
  requires:
    - context.original_prompt
    - context.problem_type
  optional_consumes:
    - context.synthesis_summary
    - insights.convergent

categories:
  - implementation
  - architecture
  - project-management

input_requirements:
  - specification
  - requirements
  - architecture-decision

output_types:
  - implementation
  - code
  - tests
  - documentation
---
```

---

## Migration Guide: v2.0 to v3.0

### Field Mapping Table

| v2.0 Field | v3.0 Field | Notes |
|------------|------------|-------|
| `handoff.accepts_handoff: true` | `handoff.accepts_from: ["*"]` | Boolean replaced with explicit source list |
| `handoff.handoff_categories` | `categories` (top-level) | Moved to top level for broader use |
| `handoff.handoff_description` | `description` (top-level, Anthropic standard) | Already exists as Anthropic-standard field |
| `handoff.handoff_trigger` | `handoff.handoff_trigger` | Preserved as legacy field (optional) |
| `handoff.protocol_version: "2.0"` | `handoff.schema_version: "3.0"` | Renamed for clarity |
| `handoff.requires` | `handoff.requires` | Preserved as legacy field (optional) |
| `handoff.optional_consumes` | `handoff.optional_consumes` | Preserved as legacy field (optional) |
| (none) | `handoff.provides_to` | New: explicit target list |
| (none) | `handoff.schema_type` | New: "universal" or "internal" |
| (none) | `input_requirements` | New: input type declarations |
| (none) | `output_types` | New: output type declarations |

### Migration Steps

1. Add `handoff.accepts_from` (replace `accepts_handoff: true` with explicit list)
2. Add `handoff.provides_to` (list known target skills)
3. Change `protocol_version` to `schema_version: "3.0"`
4. Add `handoff.schema_type: universal`
5. Move `handoff_categories` to top-level `categories`
6. Remove `handoff_description` (superseded by top-level `description`)
7. Remove `accepts_handoff` (superseded by `accepts_from`)
8. Optionally add `input_requirements` and `output_types`
9. Preserve `handoff_trigger`, `requires`, `optional_consumes` for backward compatibility

---

## Migration Period Precedence Rules

During the migration period (v2.0 and v3.0 coexistence):

1. **v3.0 fields take precedence** over v2.0 fields when both are present
2. If `handoff.accepts_from` exists, ignore `handoff.accepts_handoff`
3. If top-level `categories` exists, ignore `handoff.handoff_categories`
4. If `handoff.schema_version` exists, ignore `handoff.protocol_version`
5. Legacy fields (`handoff_trigger`, `requires`, `optional_consumes`) are preserved inside the `handoff` block for consumers that still read v2.0 format
6. A skill with ONLY v2.0 fields is still discoverable but will be flagged as "legacy" in discovery results

---

## Default Values for Optional Fields

| Field | Default | Rationale |
|-------|---------|-----------|
| `handoff.schema_type` | `"universal"` | Most skills use universal schema |
| `categories` | `[]` | No default categories |
| `input_requirements` | `[]` | No declared requirements |
| `output_types` | `[]` | No declared outputs |
| `handoff.handoff_trigger` | `null` | No legacy trigger |
| `handoff.requires` | `[]` | No required context fields |
| `handoff.optional_consumes` | `[]` | No optional context fields |
