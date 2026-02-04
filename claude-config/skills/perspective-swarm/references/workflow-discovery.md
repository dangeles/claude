# Workflow Discovery Algorithm

This document specifies how perspective-swarm discovers handoff-eligible workflows.

## Discovery Mechanism

**Approach**: Hybrid (file-system scan + optional registry)
- Primary: File-system scan of SKILL.md frontmatter
- Secondary: Optional registry override for priority/metadata

## Discovery Locations (Priority Order)

| Priority | Location | Scope |
|----------|----------|-------|
| 1 (highest) | `{cwd}/.claude/skills/*/SKILL.md` | Project-specific |
| 2 | `{git_root}/.claude/skills/*/SKILL.md` | Repository-wide |
| 3 (lowest) | `~/.claude/skills/*/SKILL.md` | Global user skills |

**Note**: If `cwd == git_root`, skip Priority 2 to avoid duplicates.

## Discovery Algorithm

```
FUNCTION discover_handoff_targets():
    candidates = []
    seen_skills = Set()  # Prevent duplicates from shadowing

    # Priority 1: Project-specific (cwd)
    cwd_skills_path = "{cwd}/.claude/skills/"
    IF exists(cwd_skills_path):
        candidates.extend(scan_skills_directory(cwd_skills_path, priority=1, seen_skills))

    # Priority 2: Git root (if different from cwd)
    TRY:
        git_root = run("git rev-parse --show-toplevel", timeout=5s)
        IF git_root != cwd:
            git_skills_path = "{git_root}/.claude/skills/"
            IF exists(git_skills_path):
                candidates.extend(scan_skills_directory(git_skills_path, priority=2, seen_skills))
    CATCH GitNotFound, NotARepository, Timeout:
        # Skip git root discovery, log debug message
        log_debug("Git root detection skipped: {error}")

    # Priority 3: Global skills
    global_skills_path = expand_path("~/.claude/skills/")
    IF exists(global_skills_path):
        candidates.extend(scan_skills_directory(global_skills_path, priority=3, seen_skills))

    RETURN candidates


FUNCTION scan_skills_directory(path, priority, seen_skills):
    results = []

    FOR each entry in listdir(path):
        # Skip non-directories
        IF not is_directory(entry): CONTINUE

        # Skip hidden directories
        IF entry.name.startswith("."): CONTINUE

        # Skip if already seen (shadowed by higher priority)
        IF entry.name in seen_skills:
            log_debug("Skill {entry.name} shadowed by higher-priority location")
            CONTINUE

        skill_md_path = "{entry}/SKILL.md"
        IF not exists(skill_md_path):
            log_debug("Skipping {entry.name}: No SKILL.md found")
            CONTINUE

        TRY:
            metadata = parse_frontmatter(skill_md_path)
        CATCH YAMLParseError as e:
            log_warning("Skill {entry.name} has invalid YAML frontmatter: {e}")
            CONTINUE

        # Check if handoff-eligible
        IF not metadata.handoff or not metadata.handoff.accepts_handoff:
            CONTINUE

        # Validate required fields
        validation = validate_handoff_metadata(metadata.handoff, entry.name)
        IF validation.has_errors:
            log_warning("Skill {entry.name} has invalid handoff metadata: {validation.errors}")
            CONTINUE

        seen_skills.add(entry.name)
        results.append({
            skill: metadata.name,
            scope: priority_to_scope(priority),  # "project" | "git-root" | "global"
            priority: priority,
            categories: metadata.handoff.handoff_categories,
            description: metadata.handoff.handoff_description,
            trigger: metadata.handoff.handoff_trigger or "{payload_path}",
            protocol_version: metadata.handoff.protocol_version or "2.0",
            health_check: metadata.handoff.health_check or null,
            path: skill_md_path
        })

    RETURN results
```

## Metadata Validation

### Required Fields for Handoff Eligibility

| Field | Type | Validation |
|-------|------|------------|
| `handoff.accepts_handoff` | boolean | Must be `true` |
| `handoff.handoff_categories` | array | Must have at least 1 element |
| `handoff.handoff_description` | string | Must be non-empty |

### Optional Fields (with defaults)

| Field | Default |
|-------|---------|
| `handoff.handoff_trigger` | `"{payload_path}"` |
| `handoff.protocol_version` | `"2.0"` |
| `handoff.health_check` | `null` (no health check) |
| `handoff.requires` | `[]` (no required payload fields) |
| `handoff.optional_consumes` | `[]` |

### Lenient Parsing Rules

When parsing metadata, apply lenient defaults:

1. **Missing `handoff_categories`**: Skip skill (required field)
2. **`handoff_categories` is string**: Convert to single-element array
3. **Missing `handoff_description`**: Use skill's main `description` field
4. **Missing `handoff_trigger`**: Default to `"{payload_path}"`

## Relevance Scoring Algorithm

```
FUNCTION score_relevance(synthesis, candidate):
    score = 0

    # 1. Category match to problem_type (+3 points)
    problem_type = synthesis.problem_type
    categories = candidate.categories

    category_mappings = {
        "strategic": ["research", "analysis", "architecture"],
        "analytical": ["analysis", "research", "verification"],
        "creative": ["creative", "research", "analysis"],
        "decision": ["research", "analysis", "verification"]
    }

    IF any(cat in categories for cat in category_mappings.get(problem_type, [])):
        score += 3

    # 2. Uncertainty signals (+2 points)
    uncertainty_count = len(synthesis.key_uncertainties)
    IF uncertainty_count >= 3 AND "research" in categories:
        score += 2

    # 3. Implementation signals (+2 points)
    implementation_keywords = ["build", "implement", "create", "develop", "code", "software"]
    IF any(kw in synthesis.synthesis_summary.lower() for kw in implementation_keywords):
        IF "implementation" in categories:
            score += 2

    # 4. Low convergence signals (+2 points)
    IF synthesis.convergence_level in ["low", "none"]:
        IF "research" in categories:
            score += 2

    # 5. Scope priority bonus (+1 for project-scope)
    IF candidate.scope == "project":
        score += 1

    RETURN score
```

## Tie-Breaking Rules

When multiple candidates have the same relevance score:

1. **Primary**: Higher priority scope wins (project > git-root > global)
2. **Secondary**: Alphabetical by skill name (stable, predictable ordering)

## Caching Strategy

### Session-Scoped Cache

Discovery runs **once** at Stage 1 (after framing), results cached for Stage 4.

**Cache file**: `{session_path}/available-workflows.yaml`

```yaml
discovery:
  timestamp: ISO8601
  ttl_minutes: 30        # Cache valid for 30 minutes

workflows:
  - skill: lit-pm
    scope: global
    priority: 3
    categories: [research, literature]
    description: "Comprehensive literature review (4-24 hours)"
    trigger: "--handoff {payload_path}"
    relevance_score: null  # Populated at Stage 4

  - skill: programming-pm
    scope: global
    priority: 3
    categories: [implementation]
    description: "Software implementation (2-8 hours)"
    trigger: "--handoff {payload_path}"
    relevance_score: null
```

### Cache Invalidation

Cache is invalidated if:
1. TTL (30 minutes) exceeded
2. User explicitly requests re-scan
3. Session is resumed from a checkpoint

## Error Handling

### EC-01: No Workflows Found

```
IF len(candidates) == 0:
    # Differentiate the message
    total_skills_scanned = count_all_skill_directories()

    IF total_skills_scanned == 0:
        message = "No skills installed. Install skills to enable handoffs."
    ELSE:
        message = "{total_skills_scanned} skills found, none accept handoffs.
                   To enable handoffs, add `handoff:` metadata to SKILL.md.
                   See: references/handoff-schema.md"
```

### EC-02: Malformed Metadata

```
ON parse_error:
    log_warning("Skill {name} skipped: {error_details}")
    # Continue processing other skills (best-effort discovery)
```

### EC-03: Git Root Detection Failures

```
TRY:
    git_root = run("git rev-parse --show-toplevel", timeout=5s)
CATCH GitNotFound:
    log_debug("git command not available, skipping git-root discovery")
    git_root = null
CATCH NotARepository:
    log_debug("Not in a git repository, skipping git-root discovery")
    git_root = null
CATCH Timeout:
    log_warning("Git root detection timed out after 5 seconds")
    git_root = null

# Handle submodule case
IF git_root is not null:
    TRY:
        superproject_root = run("git rev-parse --show-superproject-working-tree", timeout=5s)
        # --show-superproject-working-tree returns NON-EMPTY when IN a submodule
        # and returns EMPTY when NOT in a submodule
        IF superproject_root != "":
            # We are in a submodule, use the superproject root for discovery
            git_root = superproject_root
    CATCH:
        # Keep using git_root from --show-toplevel
        pass
```

## Performance Characteristics

| Metric | Expected | Notes |
|--------|----------|-------|
| Skills scanned | 30-50 | Typical ecosystem size |
| Time to scan | <1 second | File-system + YAML parsing |
| Memory usage | <10 MB | Metadata only, not full files |

**Acceptable** because discovery runs once per session (15-30 minute workflow).

## Debugging Discovery

To diagnose why a skill isn't appearing:

```
# Check if skill directory exists
ls ~/.claude/skills/{skill-name}/SKILL.md

# Check frontmatter syntax
head -50 ~/.claude/skills/{skill-name}/SKILL.md

# Required fields in frontmatter:
# handoff:
#   accepts_handoff: true
#   handoff_categories: [category1, category2]
#   handoff_description: "Description shown in menu"
```

## Shadow Detection

When a project-scope skill has the same name as a global skill:

```
log_info("Skill {name} at {project_path} shadows global skill at {global_path}")
```

The project-scope skill takes precedence. User is **not** warned (this is expected behavior for customization).

## Circular Handoff Prevention

The `handoff_chain` field in the handoff payload tracks the sequence of skills that have participated in a handoff chain within a session. This chain is built from the session history, not at discovery time.

### How handoff_chain Works

1. **Initial handoff**: When perspective-swarm first creates a handoff payload, `handoff_chain` is set to `["perspective-swarm"]`.

2. **Subsequent handoffs**: If the receiving skill (e.g., lit-pm) later hands off to another skill, it appends itself to the chain: `["perspective-swarm", "lit-pm"]`.

3. **Loop detection at handoff time**: Before generating a handoff payload, check if the target skill already appears in `handoff_chain`. If so, display a warning (but do not block):

```
WARNING: {target_skill} already appears in the handoff chain for this session.
         Continuing may create a circular workflow.
         Chain: perspective-swarm -> lit-pm -> {target_skill}
         Proceed anyway? [y/N]
```

### Important Clarification

The `handoff_chain` is **not** used during discovery (Stage 1). It is only relevant at handoff time (Stage 4) when a user has selected a specific target. This is because:

- Discovery happens before any handoff has occurred
- The chain is session-specific, not skill-specific
- Only the handoff payload carries chain state, not the discovery cache
