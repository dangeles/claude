# Decision Matrix

Reference tables for the git-strategy-advisor skill. All decision logic is defined here.

**Version**: 1.0

---

## Scope Classification

### Threshold Table

Boundaries are inclusive. A metric value that falls on a boundary belongs to the scope level whose range includes that value.

| Scope | files_changed | lines_changed | directories_spanned |
|-------|--------------|---------------|---------------------|
| trivial | == 1 | <= 10 | == 1 |
| minor | IN [1, 2] | IN [11, 50] | == 1 |
| moderate | IN [3, 5] | IN [51, 200] | IN [2, 3] |
| major | >= 6 | >= 201 | >= 4 |

### Aggregation Algorithm (Highest-Dimension Rule)

1. Evaluate each metric independently against the threshold table.
2. Take the MAXIMUM scope level across all three metrics.
3. **Single-metric exception**: If only one metric is elevated and the other two are at the lowest level (trivial), downgrade the result by one level and set confidence to "medium".

**Example 1**: 2 files, 180 lines, 1 directory.
- files_changed = 2 -> minor
- lines_changed = 180 -> moderate
- directories_spanned = 1 -> trivial
- Maximum = moderate. Only lines_changed is elevated; files is minor (not trivial). Exception does NOT apply. Result: **moderate**.

**Example 2**: 1 file, 55 lines, 1 directory.
- files_changed = 1 -> trivial
- lines_changed = 55 -> moderate
- directories_spanned = 1 -> trivial
- Maximum = moderate. Only lines is elevated, other two are trivial. Exception applies: downgrade to **minor**, confidence = "medium".

### Large Change Detection

If metrics significantly exceed major thresholds (>50 files or >5000 lines), classify as "major", set confidence to "medium", add warning:
```
code: "LARGE_CHANGE"
message: "Unusually large change. If this includes generated code or vendored dependencies, the recommendation may be more conservative than necessary."
```

---

## Work Type Classification

### Extension Mapping

| Type | Extensions |
|------|-----------|
| code | .py, .js, .ts, .go, .rs, .java, .sh, .rb, .c, .cpp, .h |
| documentation | .md, .rst, .txt |
| configuration | .yaml, .yml, .json, .toml, .cfg, .ini, .env |
| skill | Files under `claude-config/skills/` or `~/.claude/skills/` |

### Majority Rule

1. Classify each file by extension.
2. Count files per type.
3. If ALL files belong to one type: use that type.
4. If >= 80% of files belong to one type: use that type (note minority in warnings).
5. If no type reaches 80%: classify as "mixed".
6. **Skill type override**: Requires >50% skill-path files. If <= 50%, classify as "mixed".
7. Unrecognized extensions: classify as "code" (conservative default).

---

## Decision 1: Branch Strategy

| Scope | Work Type | On Primary Branch? | Action |
|-------|-----------|-------------------|--------|
| trivial | documentation | any | direct-commit |
| trivial | configuration | any | direct-commit |
| trivial | code | yes | direct-commit |
| trivial | code | no | use-current-branch |
| trivial | skill | any | direct-commit |
| trivial | mixed | yes | direct-commit |
| trivial | mixed | no | use-current-branch |
| minor | any | yes | create-feature-branch |
| minor | any | no | use-current-branch |
| moderate | any | yes | create-feature-branch |
| moderate | any | no | use-current-branch |
| major | any | any | create-feature-branch |

**Default row**: If no row matches, use `create-feature-branch` (conservative).

---

## Decision 2: Branch Naming

### Format

`{type}/{brief-description}`

### Type Mapping

| Source | Branch Prefix |
|--------|--------------|
| Work type: code | `feature/` |
| Work type: documentation | `docs/` |
| Work type: configuration | `config/` |
| Work type: skill | `skill/` |
| Work type: mixed | `feature/` |

### Keyword Overrides (take precedence over work-type mapping)

| Keywords in Description | Branch Prefix |
|------------------------|--------------|
| "fix", "bug", "patch", "repair" | `fix/` |
| "refactor", "restructure", "reorganize", "cleanup" | `refactor/` |

### Name Generation Rules

1. Extract 2-4 keywords from task description.
2. Join with hyphens, lowercase.
3. Truncate at last complete word within 40 characters (never cut mid-word).
4. Example: "Add parallel web search to researcher" -> `feature/parallel-web-search-researcher`

---

## Decision 3: Push Strategy

| Scope | Has Remote? | Action |
|-------|-------------|--------|
| trivial | yes | stay-local |
| trivial | no | stay-local |
| minor | yes | push-after-review |
| minor | no | stay-local |
| moderate | yes | push-now |
| moderate | no | stay-local |
| major | yes | push-now |
| major | no | push-now (advisory: "Consider setting up a remote for backup and collaboration.") |

---

## Decision 4: PR Strategy

| Scope | Work Type | Action |
|-------|-----------|--------|
| trivial | any | no-pr-needed |
| minor | documentation | no-pr-needed |
| minor | code | consider-pr |
| minor | configuration | consider-pr |
| minor | skill | consider-pr |
| minor | mixed | consider-pr |
| moderate | any | create-pr |
| major | any | create-pr |

For major scope, rationale includes: "Consider breaking into smaller PRs for better review quality (research shows review effectiveness drops significantly above 200 lines)."

---

## Consistency Rules

After computing all four decisions independently, apply these override rules in order:

1. **Direct-commit + PR conflict**: If branch.action = "direct-commit" AND pr.action IN ["create-pr", "consider-pr"], override pr.action to "no-pr-needed". Append to pr.rationale: "Direct commits do not use pull requests."

2. **PR requires push**: If pr.action = "create-pr" AND push.action = "stay-local", override push.action to "push-now". Append to push.rationale: "Creating a PR requires pushing to remote."

3. **No remote + PR**: If pr.action = "create-pr" AND has_remote = false, override pr.action to "consider-pr". Append to pr.rationale: "No remote configured. Set up a remote to enable pull requests." Override push.action to "stay-local".

4. **Main branch protection**: If branch.action = "use-current-branch" AND is_on_main = true AND scope IN ["moderate", "major"], override branch.action to "create-feature-branch". Append to branch.rationale: "Moderate/major changes should not be committed directly to the primary branch."

5. **Already committed**: If branch.action = "create-feature-branch" AND mode = "post-work" AND working tree is clean with unpushed commits on current branch, override branch.action to "use-current-branch". Append to branch.rationale: "Changes are already committed on the current branch."

Log any overrides in the `warnings` array with code "CONSISTENCY_OVERRIDE".

---

## Output Schema Versioning

**Current Version**: 1.0

- **Minor versions** (1.0 -> 1.1): Additive only. New fields may appear. Consumers MUST ignore unknown fields.
- **Major versions** (1.x -> 2.0): Breaking changes. Migration guidance provided.
