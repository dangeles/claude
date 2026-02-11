# Handoff Registry Reference

**Version**: 1.0
**Schema**: Universal Handoff Schema v3.0
**Date**: 2026-02-07

---

## Overview

The handoff registry is a **lightweight metadata system** that enables workflows to:

1. Discover other workflows that can receive handoffs
2. Determine handoff compatibility (schema versions)
3. Find workflows by category (implementation, research, analysis, etc.)

**Key Principle**: Registry is **optional** -- workflows can hand off without it using frontmatter-based discovery. Registry provides enhanced discovery and centralized metadata.

---

## Registry File: Optional

The registry file (`~/.claude/handoff-registry.yaml`) is OPTIONAL. The system operates correctly with frontmatter-only discovery. If the registry file does not exist, discovery falls through to frontmatter scan silently. No error, no warning.

---

## Architecture: Hybrid Discovery

### Discovery Modes

**Mode 1: Frontmatter Scan** (Default, no registry required)

- Scan `~/.claude/skills/*/SKILL.md` for handoff metadata in YAML frontmatter
- Parse frontmatter for `handoff.accepts_from`, `handoff.provides_to`, `categories`, etc.
- Build in-memory discovery index
- Cache for session duration
- **Pros**: No central registry needed, always up-to-date
- **Cons**: Slower discovery (file I/O), no cross-machine coordination

**Mode 2: Registry File** (Optional enhancement)

- Read `~/.claude/handoff-registry.yaml`
- Pre-computed workflow metadata
- Faster than frontmatter scan
- **Pros**: Fast lookup, can include remote workflows
- **Cons**: Needs sync when workflows change

**Hybrid Approach**: Try registry first (if exists), fall back to frontmatter scan.

---

## Registry File Format

### Location

```
~/.claude/handoff-registry.yaml
```

### Structure

```yaml
version: "1.0"
updated_at: "2026-02-07T18:30:00Z"

workflows:
  - name: "programming-pm"
    version: "1.2"
    description: "Multi-phase project management workflow for software development"
    skill_path: "~/.claude/skills/programming-pm/SKILL.md"
    handoff:
      accepts_from: ["*"]
      provides_to:
        - "senior-developer"
        - "junior-developer"
        - "systems-architect"
      schema_version: "3.0"
      schema_type: "universal"
    categories:
      - "implementation"
      - "architecture"
    input_requirements:
      - "specification"
      - "requirements"
    output_types:
      - "implementation"
      - "code"
      - "tests"
    meta:
      estimated_duration: "days-to-weeks"
      complexity: "high"
      autonomous: false

  - name: "skill-editor"
    version: "1.0"
    description: "Multi-agent workflow for creating and modifying Claude Code skills"
    skill_path: "~/.claude/skills/skill-editor/SKILL.md"
    handoff:
      accepts_from: ["*"]
      provides_to:
        - "programming-pm"
        - "technical-pm"
      schema_version: "3.0"
      schema_type: "universal"
    categories:
      - "skill-development"
      - "workflow-creation"
    input_requirements:
      - "specification"
      - "skill-request"
    output_types:
      - "skill"
      - "agent-configuration"
    meta:
      estimated_duration: "hours"
      complexity: "high"
      autonomous: true

meta:
  total_workflows: 2
  last_scan: "2026-02-07T18:30:00Z"
  discovery_mode: "hybrid"
```

---

## Discovery Algorithm

### Step 1: Load Registry (if exists)

Check for `~/.claude/handoff-registry.yaml`. If it exists, parse it. If not, proceed to Step 2.

### Step 2: Frontmatter Scan (fallback or supplement)

Scan `~/.claude/skills/*/SKILL.md` for YAML frontmatter containing a `handoff` block. For each skill with handoff metadata, extract:

- `name`
- `handoff.accepts_from`
- `handoff.provides_to`
- `handoff.schema_version`
- `handoff.schema_type`
- `categories`
- `input_requirements`
- `output_types`

### Step 3: Merge and Cache

If both registry and frontmatter data exist, merge them. Frontmatter takes precedence (always up-to-date). Cache the merged result for the session duration.

### Step 4: Query

Query the discovery index by:

- **Category**: Find workflows in a specific category
- **Accepts from**: Find workflows that accept handoffs from a specific skill
- **Schema version**: Find workflows supporting a specific schema version
- **Input requirements**: Find workflows that accept specific input types

---

## Relevance Scoring

When multiple workflows match discovery criteria, rank by relevance:

### Scoring Algorithm

| Factor | Points | Description |
|--------|--------|-------------|
| Category match | 0-40 | Direct match = 40, partial = 20 |
| Schema compatibility | 0-20 | v3.0 = 20, universal support = 10 |
| Input compatibility | 0-20 | Direct match = 20, partial = 10 |
| Previous success | 0-10 | Successful handoff history = 10 |
| Autonomy match | 0-10 | Matches supervision preference = 10 |
| Loop penalty | -30 | Target already in handoff_chain |

### Selection

Present top 3-5 workflows to the user, ranked by score. Include:

- Workflow name and description
- Relevance score
- Schema compatibility status
- Estimated duration

---

## Schema Version Compatibility

### Compatibility Matrix

| Source Schema | Target Schema | Compatible? | Notes |
|--------------|---------------|-------------|-------|
| 3.0 | 3.0 | Yes | Perfect match |
| 3.0 | 2.0 | Partial | Downgrade (lose tracing, tiered context) |
| 3.0 | 1.2 | No | programming-pm internal schema (incompatible) |
| 2.0 | 3.0 | Yes | Upgrade path exists |
| 1.2 | 3.0 | Adapter needed | Requires adapter layer |

### Adapter Layer (Future)

For incompatible schemas, an adapter translates between schema versions. Adapters are deferred to v1.1.

---

## Cache Management

Discovery results are cached for the duration of a Claude session. To force re-discovery when skills are added or modified:

1. New session automatically re-scans
2. Within a session, re-invoke discovery explicitly

---

## Registry Maintenance

### Automatic Sync

To regenerate the registry from frontmatter:

1. Scan all `~/.claude/skills/*/SKILL.md` files
2. Extract handoff metadata from frontmatter
3. Merge with existing registry (preserve custom metadata)
4. Write updated registry file
5. Validate with YAML parser

### Manual Editing

Users can manually edit `~/.claude/handoff-registry.yaml` to:

- Add remote workflows (not in local `~/.claude/skills/`)
- Override discovery settings
- Add custom metadata (notes, tags)

Manual edits are preserved during sync (merge strategy).

---

## Example: Find Workflows for Implementation Handoff

Given a skill-editor session that needs to hand off implementation work:

1. Query: `find_workflows(category="implementation", accepts_from="skill-editor")`
2. Results:
   - `programming-pm` (score: 70) -- Full implementation pipeline, accepts from any
   - `technical-pm` (score: 45) -- Lightweight coordination, accepts from any
3. Present to user for selection

---

## References

- **Frontmatter Metadata Standard**: `./frontmatter-metadata-standard.md`
- **Universal Handoff Schema v3.0**: `./universal-handoff-schema-v3.0.json`
- **perspective-swarm discovery**: `../../perspective-swarm/references/workflow-discovery.md`
