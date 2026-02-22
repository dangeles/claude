# Configuration Management Workflow

Instructions for Claude Code Agent when modifying global configuration (skills, settings, plugins) in this repository.

## When This Applies

Follow this workflow when modifying:
- Skills in `claude-config/skills/`
- Settings in `claude-config/settings.json`
- Plugin configuration in `claude-config/plugins/`

## Important Context

This repository tracks **global** `~/.claude/` configuration using a sync system:
- `claude-config/` contains synced copy of `~/.claude/`
- `sync-config.py` manages bidirectional sync
- Changes here affect your global Claude Code configuration across all projects

## Workflow

### 1. Pre-Modification Safety Check

Before making any configuration changes:

```bash
# Check git status
git status

# Ensure working directory is clean
# If uncommitted changes exist, commit them first
git add .
git commit -m "Save work before config change"

# Check sync status
./sync-config.py status
```

### 2. Create Planning Journal Entry

Document the intended change:

```bash
./sync-config.py plan --title "Description of config change"
```

Or manually create: `planning/[hostname]/[YYYY-MM-DD]-[description].md`

Fill in:
- **Objective**: What are you changing and why?
- **Changes Planned**: Specific files/modifications
- **Expected Outcome**: What should improve?

### 3. Implement Changes in Repository

**Key principle**: Modify `claude-config/` in repository first, NEVER modify `~/.claude/` directly.

This ensures changes are version-controlled before affecting live system.

**Examples**:
- Modify skill: Edit `claude-config/skills/[skill-name]/SKILL.md`
- Update settings: Edit `claude-config/settings.json`
- Add new skill: Create `claude-config/skills/[new-skill]/SKILL.md`

### 4. Quality Analysis

Before syncing to live system, verify:

**Logical Consistency:**
- [ ] No conflicting plugin dependencies
- [ ] Skill prerequisites are satisfied
- [ ] Settings don't contradict each other

**Correctness:**
- [ ] Valid JSON/YAML syntax
  ```bash
  jq . claude-config/settings.json  # Validate JSON
  ```
- [ ] Proper SKILL.md structure (YAML frontmatter + markdown sections)
- [ ] File naming follows conventions
- [ ] Examples are accurate and runnable

**Succinctness:**
- [ ] No unnecessary duplication
- [ ] Clear and concise documentation
- [ ] Follows existing patterns
- [ ] Minimal file size (avoid bloat)

### 5. Preview and Sync to System

Preview what will be synced:

```bash
./sync-config.py push --dry-run
```

Review the preview carefully. If correct:

```bash
./sync-config.py push
```

This applies changes to `~/.claude/`, making them live in your Claude Code environment.

### 6. Test Changes

Test the modification:
- **Skills**: Invoke the skill in a Claude Code session
- **Settings**: Restart Claude Code if needed, verify settings applied
- **Plugins**: Check plugin functionality

Document test results in planning journal entry.

### 7. Commit or Revert

Based on testing results:

**If changes work correctly:**

```bash
git add claude-config/ planning/
git commit -m "Brief description

Detailed explanation of what was changed and why.

See planning/[hostname]/[YYYY-MM-DD]-[description].md

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"

git push
```

Update planning journal entry:
- Mark status as "Success"
- Document actual outcome
- Note any lessons learned
- Add commit SHA to "Related Commits" section

**If changes have issues:**

Ask user: "Configuration changes have issues. Should I:
1. Fix the issues and retry?
2. Revert to last known good state?"

If reverting:

```bash
# Reset repository changes
git reset --hard HEAD

# Restore ~/.claude/ from last good state
./sync-config.py push

# Update planning journal to document what didn't work
```

## Command Quick Reference

```bash
# Pre-check
git status
./sync-config.py status

# Create planning entry
./sync-config.py plan --title "Description"

# Create planning entry and open in editor
./sync-config.py plan --title "Description" --edit

# Validate JSON
jq . claude-config/settings.json

# Preview sync
./sync-config.py push --dry-run

# Apply to ~/.claude/
./sync-config.py push

# Non-interactive (for agents/scripts)
./sync-config.py push --yes             # Push without prompts (source wins)
./sync-config.py push --yes --delete    # Push and remove orphaned files
./sync-config.py push --yes --dry-run --delete  # Preview deletions

# Commit changes
git add claude-config/ planning/
git commit -m "Description"
git push

# Revert if needed
git reset --hard HEAD
./sync-config.py push
```

## Multi-Machine Workflow

**On Machine A** (where you make changes):
1. Modify `claude-config/` in repository
2. Analyze and validate changes
3. Sync to `~/.claude/`: `./sync-config.py push`
4. Test locally
5. Commit and push to remote

**On Machine B** (receiving changes):
1. `git pull` to get latest configuration
2. `./sync-config.py push` to apply to `~/.claude/`
3. Restart Claude Code if needed
4. New configuration is now active

## Rollback Procedure

If something goes wrong:

1. **Identify last known good commit**:
   ```bash
   git log --oneline claude-config/ planning/
   ```

2. **Reset to that commit**:
   ```bash
   git reset --hard [commit-sha]
   ```

3. **Restore ~/.claude/ from that state**:
   ```bash
   ./sync-config.py push
   ```

4. **Document rollback in planning journal**:
   - Edit the failed planning entry
   - Add to "Issues" section: reason for rollback
   - Mark status as "Failed"

## Examples

### Example 1: Add Example to Existing Skill

**User request**: "Add investment calculation example to calculator skill"

**Workflow**:
1. Pre-check: `git status` -> clean
2. Planning: `./sync-config.py plan --title "Add investment example to calculator"`
3. Implementation:
   ```bash
   cat > claude-config/skills/calculator/examples/investment-example.md << 'EOF'
   # Investment Calculation Example
   [content...]
   EOF
   ```
4. Analysis: Check format matches existing examples
5. Preview: `./sync-config.py push --dry-run` -> shows calculator will update
6. Sync: `./sync-config.py push`
7. Test: Invoke calculator skill, verify example appears
8. Commit:
   ```bash
   git add claude-config/skills/calculator/ planning/
   git commit -m "Add investment calculation example to calculator skill"
   git push
   ```
9. Update planning: Mark as Success, add commit SHA

### Example 2: Create New Skill

**User request**: "Create skill for API documentation generation"

**Workflow**:
1. Pre-check: `git status` -> clean
2. Planning: `./sync-config.py plan --title "Create API documentation skill"`
3. Implementation:
   ```bash
   mkdir -p claude-config/skills/api-doc-generator/{examples,references}

   cat > claude-config/skills/api-doc-generator/SKILL.md << 'EOF'
   ---
   name: api-doc-generator
   version: 1.0.0
   description: Generate comprehensive API documentation
   ---

   # API Documentation Generator

   [skill content...]
   EOF
   ```
4. Analysis: Validate YAML frontmatter, check structure
5. Preview: `./sync-config.py push --dry-run`
6. Sync: `./sync-config.py push`
7. Test: Invoke skill in Claude Code session
8. Commit: `git add claude-config/skills/api-doc-generator/ planning/` + commit
9. Update planning: Document outcome

### Example 3: Update Settings

**User request**: "Increase max tokens globally"

**Workflow**:
1. Pre-check: `git status` -> clean
2. Planning: `./sync-config.py plan --title "Increase max tokens"`
3. Implementation:
   ```bash
   jq '.maxTokens = 200000' claude-config/settings.json > tmp.json
   mv tmp.json claude-config/settings.json
   ```
4. Analysis: Validate JSON, value is reasonable
5. Preview: `./sync-config.py push --dry-run`
6. Sync: `./sync-config.py push`
7. Test: Restart Claude Code, verify new token limit
8. Commit: `git add claude-config/settings.json planning/` + commit
9. Update planning: Success

## Quality Checklist

Before syncing changes to `~/.claude/`:

**Structure:**
- [ ] Files are in correct `claude-config/` subdirectories
- [ ] Naming conventions followed
- [ ] Directory structure matches expected patterns

**Syntax:**
- [ ] JSON validated with `jq . file.json`
- [ ] YAML validated (if using frontmatter)
- [ ] Markdown linted (if applicable)

**Logic:**
- [ ] No conflicting settings
- [ ] Dependencies are satisfied
- [ ] Prerequisites documented

**Documentation:**
- [ ] Planning journal entry created
- [ ] Objective clearly stated
- [ ] Expected outcome documented
- [ ] Changes explained in commit message

**Testing:**
- [ ] Tested in Claude Code session
- [ ] Verified expected behavior
- [ ] Documented test results in planning entry

## Integration with Existing Tools

**With sync-config.py:**
- `status` - Check divergence and detect orphaned files in ~/.claude/
- `pull` - Sync changes from ~/.claude/ to repository (for manual changes)
- `push --dry-run` - Preview what will be synced
- `push` - Apply repository changes to ~/.claude/
- `push --yes` - Apply without interactive prompts (source always wins)
- `push --delete` - Remove files in ~/.claude/ that have no repo counterpart
- `push --yes --delete` - Non-interactive push with orphan cleanup (for agents)
- `plan --title` - Create planning journal entry (prints path)
- `plan --title --edit` - Create planning journal entry and open in editor

**With Planning Journal:**
- Document all configuration changes
- Track outcomes and lessons learned
- Link to git commit SHAs
- Organize by machine hostname

**With Git:**
- Repository is source of truth
- All changes version controlled
- Commit before syncing to live system
- Descriptive messages with planning journal references

## Anti-Patterns to Avoid

-- **Don't** modify `~/.claude/` directly - always modify `claude-config/` in repository first
-- **Don't** sync without testing - always preview with `--dry-run`
-- **Don't** skip planning journal - always document why changes are made
-- **Don't** commit untested changes - always verify in Claude Code first
-- **Don't** forget to sync to ~/.claude/ after committing - changes won't take effect
-- **Don't** make changes on multiple machines simultaneously - causes conflicts
-- **Don't** use `--delete` without `--dry-run` first - always preview what will be removed
-- **Don't** use `push --yes --delete` in scripts without validating exit code - it may abort on safety threshold

## Special Cases

### Handling Sync Conflicts

If `~/.claude/` was modified outside the repository:

```bash
# Pull changes from ~/.claude/ to repository
./sync-config.py pull

# Review what changed
git diff

# If changes are good, commit them
git add claude-config/
git commit -m "Sync manual changes from ~/.claude/"

# Then proceed with your modifications
```

### Emergency Rollback

If Claude Code is broken after a configuration change:

```bash
# Immediately revert to last commit
git reset --hard HEAD^

# Restore ~/.claude/ from that state
./sync-config.py push

# Restart Claude Code

# Document what went wrong in planning journal
```

### Team Collaboration

When multiple people share this configuration:

1. Create feature branch for changes
2. Make modifications in `claude-config/`
3. Test thoroughly before merging
4. Team reviews changes via PR
5. After merge, each team member runs `git pull && ./sync-config.py push`

## Summary

**Key Principles:**
1. **Repository first**: Always modify `claude-config/`, never `~/.claude/` directly
2. **Plan always**: Create planning journal entry before changes
3. **Analyze carefully**: Quality check before syncing
4. **Preview before sync**: Use `--dry-run` to preview
5. **Test before commit**: Verify in Claude Code session
6. **Document outcomes**: Update planning journal with results
7. **Commit and sync**: Commit to git, then ensure other machines sync

**Workflow Summary:**
```
git status -> create plan -> modify claude-config/ -> analyze ->
preview sync -> sync to ~/.claude/ -> test -> commit -> update plan -> push
```

This workflow ensures all global configuration changes are:
- Version controlled
- Analyzed for quality
- Tested before deployment
- Documented for future reference
- Reversible if issues arise
