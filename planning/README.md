# Planning Journal

This directory contains a journal system for tracking Claude Code configuration changes and assessing their impact.

## Purpose

The planning journal helps you:

1. **Document configuration changes** before implementing them
2. **Track outcomes** of configuration experiments
3. **Learn from experience** by recording what worked and what didn't
4. **Maintain machine-specific context** organized by hostname
5. **Link plans to commits** for full traceability

## Structure

```
planning/
├── README.md                      # This file
├── .template.md                   # Template for new entries
└── [hostname]/                    # Machine-specific directories
    └── YYYY-MM-DD-description.md  # Timestamped planning entries
```

Each machine gets its own directory (using the hostname), allowing you to track configuration changes per-machine while maintaining a unified repository.

## Workflow

### 1. Creating a Planning Entry

Before making configuration changes:

```bash
./sync-config.py plan --title "Enable scientific skills plugin"
```

This creates a new entry with:
- Current date
- Auto-detected hostname
- Pre-filled template
- Status set to "Planned"

The file opens in your editor (using `$EDITOR` environment variable).

### 2. During Implementation

As you make changes:

- ✅ Check off items in the "Changes Planned" section
- 📝 Note any issues encountered
- 🔗 Record git commit SHAs in "Related Commits"
- 💭 Update status to "In Progress"

### 3. After Completion

Document the outcome:

```bash
./sync-config.py plan --title "Enable scientific skills plugin" --status "Success"
```

Or manually edit the entry to add:

- **Actual Outcome**: What really happened
- **Assessment**: Success/Partial/Failed
- **Improvements**: What got better
- **Issues**: Problems that emerged
- **Lessons Learned**: What you'd do differently
- **Next Steps**: Follow-up actions needed

### 4. Reviewing History

List all planning entries:

```bash
./sync-config.py plan --list
```

Filter by machine:

```bash
./sync-config.py plan --list --machine macbook-pro
```

## Template Fields

### Objective
What configuration change are you planning? What problem are you solving?

Example:
```
Enable the scientific-skills plugin to access bioinformatics tools
and database integrations for research projects.
```

### Changes Planned
Checklist of specific changes you'll make.

Example:
```
- [ ] Update settings.json to enable scientific-skills plugin
- [ ] Sync configuration to repository
- [ ] Test skills are available in Claude Code
- [ ] Update project-specific .claude-project.json if needed
```

### Expected Outcome
What improvement do you expect from these changes?

Example:
```
Will have access to 100+ scientific skills including:
- Database skills (PubMed, PDB, UniProt, etc.)
- Analysis skills (plotting, statistics, ML)
- Tool integrations (Benchling, DNAnexus, etc.)
```

### Actual Outcome
After implementation: What actually happened?

Example:
```
Successfully enabled plugin. All skills now available in Claude Code.
Tested with /pubmed-database skill - works perfectly.
```

### Assessment

**Result**: Success / Partial / Failed

**Improvements**:
List concrete improvements or benefits realized.

**Issues**:
Problems encountered during or after implementation.

**Lessons Learned**:
What would you do differently next time?

Example:
```
Result: Success

Improvements:
- Access to 170+ scientific skills
- Can now search PubMed directly from Claude Code
- Bioinformatics workflows much easier

Issues:
- Initial sync took ~5 minutes (1.1 MB of skills)
- Some skills require additional setup (API keys)

Lessons Learned:
- Always check skill documentation before first use
- Consider creating project-specific .claude-project.json to enable
  only relevant skills for each project
```

### Related Commits
Link the planning entry to git commits.

Example:
```
- abc1234: Enable scientific-skills plugin in settings.json
- def5678: Sync skills directory from ~/.claude
```

### Next Steps
Follow-up actions or adaptations needed.

Example:
```
- Set up API keys for external services (PubMed, etc.)
- Create project-specific config for bioreactor project
- Document commonly-used skills in project README
```

## Best Practices

### 1. Create Entries Before Changes
Always create a planning entry BEFORE making significant configuration changes. This helps clarify your thinking and provides a record of intent.

### 2. Be Specific
Instead of "Enable new plugin", write "Enable scientific-skills plugin for bioinformatics research".

### 3. Document Issues
If something goes wrong, document it! Failed experiments are valuable learning experiences.

### 4. Update Status
Keep the status field current:
- **Planned**: Entry created, not yet implemented
- **In Progress**: Currently implementing changes
- **Success**: Completed successfully
- **Partial**: Partially successful, some issues remain
- **Failed**: Did not achieve intended outcome

### 5. Link to Commits
Always record related git commit SHAs. This provides traceability between plans and implementation.

### 6. Review Periodically
Periodically review your planning journal to:
- Identify patterns (what works, what doesn't)
- Improve future configuration decisions
- Share lessons with collaborators

## Examples

### Good Entry Example

```markdown
# Enable Scientific Skills Plugin for Research Projects

**Date**: 2026-02-02
**Machine**: macbook-pro
**Status**: Success

## Objective

Enable the scientific-skills plugin to access bioinformatics tools,
database integrations, and analysis capabilities for ongoing research
projects (bioreactor optimization, protein design).

## Changes Planned

- [x] Update settings.json to enable scientific-skills plugin
- [x] Sync configuration to repository
- [x] Test core skills (pubmed-database, uniprot-database)
- [x] Document available skills in project README

## Expected Outcome

Access to 170+ scientific skills covering databases, tools, and
analysis capabilities. Should streamline literature review and
data analysis workflows.

## Actual Outcome

Successfully enabled plugin. All 170+ skills now available.
Tested several database skills - all working correctly.

## Assessment

**Result**: Success

**Improvements**:
- Direct PubMed searches from Claude Code
- Protein structure lookups (PDB, AlphaFold)
- Data analysis skills (pandas, scikit-learn, etc.)
- 50+ database integrations

**Issues**:
- Initial sync took 5 minutes (1.1 MB)
- Some skills need API keys (will set up separately)

**Lessons Learned**:
- Check skill documentation before first use
- Consider project-specific configs to reduce clutter
- Skills directory is large - exclude from frequent syncs

## Related Commits

- a1b2c3d: Enable scientific-skills plugin in settings.json
- e4f5g6h: Initial sync of skills directory

## Next Steps

- Set up API keys for external databases
- Create bioreactor project config
- Test hypothesis-generation skill
```

### Bad Entry Example

```markdown
# Update config

**Date**: 2026-02-02
**Machine**: macbook-pro
**Status**: Success

## Objective

Change some settings

## Changes Planned

- [ ] Update stuff

## Expected Outcome

Better

## Actual Outcome

Done

## Assessment

**Result**: Success

**Improvements**:
- Works

**Issues**:
- None

**Lessons Learned**:
- It's good
```

**Why this is bad**:
- Vague title ("Update config")
- No specific objective
- No actionable checklist
- No detail about what changed
- No traceability (no commits linked)
- No useful information for future reference

## Machine-Specific Organization

Each machine gets its own directory to allow:

1. **Machine-specific context**: Different machines may have different config needs
2. **Easy filtering**: `--machine` flag to view entries for specific machine
3. **Distributed workflow**: Multiple people/machines contributing to same repo
4. **Migration tracking**: Document configuration evolution across machines

## Integration with Git

The planning journal integrates naturally with git:

```bash
# Create plan
./sync-config.py plan --title "Enable new plugin"

# Edit plan file
# ... implement changes ...

# Sync config
./sync-config.py pull
git add claude-config/ planning/
git commit -m "Enable scientific-skills plugin

See planning/macbook-pro/2026-02-02-enable-new-plugin.md for details"

# Update plan with commit SHA
# Edit planning entry to add commit SHA

# Push to remote
git push
```

On another machine:

```bash
# Pull latest
git pull

# Review plan
cat planning/macbook-pro/2026-02-02-enable-new-plugin.md

# Apply config
./sync-config.py push
```

## Tips

### Finding Your Hostname

```bash
hostname
# or
./sync-config.py plan --list  # Shows all hostnames with entries
```

### Searching Entries

Use standard Unix tools:

```bash
# Search for keyword across all entries
grep -r "plugin" planning/

# Find entries by date
ls planning/*/2026-02-*

# List entries for specific machine
ls planning/macbook-pro/
```

### Editor Configuration

Set your preferred editor:

```bash
export EDITOR=vim
# or
export EDITOR="code --wait"  # VS Code
# or
export EDITOR=nano
```

Add to your shell profile (`~/.bashrc`, `~/.zshrc`) to persist.

## Summary

The planning journal is a lightweight system for:
- ✅ Documenting configuration changes
- ✅ Tracking outcomes and lessons learned
- ✅ Maintaining machine-specific context
- ✅ Linking plans to implementation (git commits)
- ✅ Building institutional knowledge over time

Use it consistently, and it will become an invaluable record of your configuration evolution.
