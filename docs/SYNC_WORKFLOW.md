# Claude Code Configuration Sync Workflows

Practical workflows and best practices for managing Claude Code configuration across multiple machines.

## Table of Contents

- [Common Workflows](#common-workflows)
- [Multi-Machine Scenarios](#multi-machine-scenarios)
- [Best Practices](#best-practices)
- [Advanced Scenarios](#advanced-scenarios)

## Common Workflows

### Workflow 1: Enable a New Plugin

You want to enable a new plugin and sync the configuration to other machines.

**On Machine A** (where you enable the plugin):

```bash
# 1. Create planning entry
./sync-config.py plan --title "Enable scientific-skills plugin"

# 2. Enable plugin in Claude Code
#    (Use Claude Code's plugin management interface)

# 3. Pull the new configuration
./sync-config.py pull

# 4. Review changes
git diff claude-config/

# 5. Commit and push
git add claude-config/ planning/
git commit -m "Enable scientific-skills plugin

Provides access to 170+ scientific skills for bioinformatics,
databases, and analysis tools.

See planning/macbook-pro/2026-02-02-enable-scientific-skills.md"
git push
```

**On Machine B** (applying the change):

```bash
# 1. Pull latest changes
git pull

# 2. Review planning entry
cat planning/macbook-pro/2026-02-02-enable-scientific-skills.md

# 3. Preview changes
./sync-config.py push --dry-run

# 4. Apply configuration
./sync-config.py push

# 5. Restart Claude Code (if running)
#    Plugin will download automatically on first use
```

**Total time**: ~5-10 minutes (including plugin download)

### Workflow 2: Modify Skill Configuration

You want to customize a skill or add a new custom skill.

**Steps**:

```bash
# 1. Document your plan
./sync-config.py plan --title "Customize systematic-debugging skill"

# 2. Edit skill file
nano ~/.claude/skills/superpowers:systematic-debugging/SKILL.md

# 3. Test the changes in Claude Code
#    (Start a session and test the skill)

# 4. Pull changes if satisfied
./sync-config.py pull

# 5. Review and commit
git diff claude-config/skills/
git add claude-config/skills/
git commit -m "Customize systematic-debugging skill for Python projects"
git push

# 6. Update planning entry with outcome
nano planning/$(hostname)/2026-02-02-customize-systematic-debugging.md
git add planning/
git commit -m "Update planning entry with outcome"
git push
```

### Workflow 3: Set Up New Machine

You're setting up Claude Code on a new machine and want to replicate your configuration.

**Steps**:

```bash
# 1. Install Claude Code
#    (Follow official installation instructions)

# 2. Clone configuration repository
git clone <repo-url> ~/repos/claude
cd ~/repos/claude

# 3. Install dependencies
pip install pyyaml

# 4. Preview what will be pushed
./sync-config.py push --dry-run

# 5. Apply configuration
./sync-config.py push

# 6. Restart Claude Code
#    Plugins will download automatically

# 7. Verify setup
./sync-config.py status
```

**Expected behavior**:
- Settings synced immediately
- Skills available immediately
- Plugins download on first use (~5 minutes for large plugins)

### Workflow 4: Project-Specific Configuration

You want to configure Claude Code settings for a specific project.

**Steps**:

```bash
# 1. Navigate to project
cd ~/repos/bioreactor

# 2. Create project config
cat > .claude-project.json << 'EOF'
{
  "enabledSkills": [
    "scientific-skills:pubmed-database",
    "scientific-skills:biopython",
    "superpowers:test-driven-development"
  ],
  "projectContext": "Bioreactor optimization project"
}
EOF

# 3. Sync to repository
cd ~/repos/claude
./sync-config.py pull-projects bioreactor

# 4. Commit
git add project-configs/bioreactor/
git commit -m "Add Claude Code config for bioreactor project"
git push
```

**On another machine**:

```bash
# 1. Pull latest
git pull

# 2. Restore project config
./sync-config.py push-projects bioreactor
```

### Workflow 5: Resolve Configuration Conflict

Two machines have diverged, and you need to resolve conflicts.

**Scenario**: Machine A enabled plugin X, Machine B enabled plugin Y.

**Machine A**:
```bash
# Pull changes from Machine B
git pull

# Attempt to push (will show conflicts)
./sync-config.py push
```

**Conflict prompt**:
```
Conflict detected:
  Source: ~/.claude/settings.json (Machine A's version)
  Target: ./claude-config/settings.json (Machine B's version)

[Diff shows Machine A has plugin X, repo has plugin Y]

Choose resolution:
[1] Use source (enable plugin X, disable Y)
[2] Use target (keep plugin Y, disable X)
[3] Show full diff
[4] Skip this file
[5] Abort sync
```

**Resolution options**:

**Option 1: Merge manually**
```bash
# Choose [5] to abort
# Manually merge in editor
nano ~/.claude/settings.json
nano claude-config/settings.json
# Enable BOTH plugins X and Y
# Then pull to sync
./sync-config.py pull
```

**Option 2: Use one version**
```bash
# Choose [1] or [2] depending on which config you prefer
# This will overwrite the other version
```

**Best practice**: If both changes are intentional, merge manually to keep both.

## Multi-Machine Scenarios

### Scenario 1: Work + Personal Machines

**Setup**: Two machines with different plugin needs.

**Approach**: Use project-specific configs

**Work machine** (macbook-pro-work):
```bash
# Enable work-related plugins only
# settings.json: minimal global plugins
# Use project configs for work-specific needs

# Example: Data science project
cd ~/repos/data-project
cat > .claude-project.json << 'EOF'
{
  "enabledSkills": [
    "scientific-skills:pandas",
    "scientific-skills:scikit-learn",
    "scientific-skills:matplotlib"
  ]
}
EOF
```

**Personal machine** (macbook-pro-personal):
```bash
# Enable personal project plugins
# settings.json: different plugins

# Sync only the project configs you need
./sync-config.py push-projects data-project
```

**Benefit**: Global settings can differ, but project-specific configs are shared.

### Scenario 2: Linux + macOS Machines

**Challenge**: Path differences, potential line ending issues.

**Solution**:

1. **Configure git** on both machines:
   ```bash
   git config core.autocrlf input  # Consistent line endings
   ```

2. **Use relative paths** in configs (Claude Code handles this)

3. **Test both directions**:
   ```bash
   # On Linux
   ./sync-config.py pull
   git push

   # On macOS
   git pull
   ./sync-config.py push
   ```

4. **Document machine-specific issues** in planning journal:
   ```bash
   ./sync-config.py plan --title "Cross-platform sync testing"
   ```

### Scenario 3: Laptop + Remote Server

**Setup**: Laptop for development, remote server for heavy computation.

**Approach**: Selective sync

**Laptop** (full configuration):
```bash
# Enable all plugins
./sync-config.py pull
git push
```

**Remote server** (minimal configuration):
```bash
# Clone repo
git clone <repo-url> ~/repos/claude

# Selectively push only needed parts
# Option 1: Use dry-run to see what would sync
./sync-config.py push --dry-run

# Option 2: Manually copy specific files
cp claude-config/settings.json ~/.claude/settings.json

# Option 3: Create server-specific config
# (Modify sync.config.yaml with different sync rules)
```

### Scenario 4: Team Collaboration

**Setup**: Multiple team members sharing configurations.

**Approach**: Branch-based workflow

**Team member A**:
```bash
# Create feature branch for config change
git checkout -b add-scientific-skills

# Make changes
./sync-config.py plan --title "Add scientific skills for team"
./sync-config.py pull

# Commit and push
git add claude-config/ planning/
git commit -m "Add scientific-skills plugin"
git push origin add-scientific-skills

# Create PR for team review
```

**Team member B**:
```bash
# Review PR
gh pr view 123

# Test config change
git fetch origin add-scientific-skills
git checkout add-scientific-skills
./sync-config.py push --dry-run
./sync-config.py push

# Test in Claude Code
# Provide feedback on PR
```

**After merge**:
```bash
# All team members update
git pull origin main
./sync-config.py push
```

## Best Practices

### 1. Always Use Planning Journal

**Why**: Documents intent and outcomes, builds institutional knowledge.

**When**: Before ANY configuration change.

**Example**:
```bash
# Instead of:
./sync-config.py pull
git commit -m "Update config"

# Do:
./sync-config.py plan --title "Enable scientific-skills plugin"
# ... make changes ...
./sync-config.py pull
git commit -m "Enable scientific-skills plugin

See planning/macbook-pro/2026-02-02-enable-scientific-skills.md"
# ... update planning entry with outcomes ...
```

### 2. Preview Before Executing

**Why**: Avoid accidents, understand what will change.

**Pattern**:
```bash
# Always preview first
./sync-config.py pull --dry-run
./sync-config.py status

# Then execute
./sync-config.py pull
```

### 3. Commit Configuration Changes Separately

**Why**: Makes git history clear and reviewable.

**Good**:
```bash
git add claude-config/
git commit -m "Enable scientific-skills plugin"

git add planning/
git commit -m "Add planning entry for plugin change"
```

**Bad**:
```bash
git add .
git commit -m "Various changes"
```

### 4. Use Descriptive Commit Messages

**Why**: Makes history searchable and understandable.

**Good**:
```
Enable scientific-skills plugin

Provides access to:
- 50+ database skills (PubMed, UniProt, PDB, etc.)
- Analysis skills (pandas, scikit-learn, matplotlib)
- Tool integrations (Benchling, DNAnexus)

See planning/macbook-pro/2026-02-02-enable-scientific-skills.md
for detailed assessment.
```

**Bad**:
```
Update config
```

### 5. Sync Regularly

**Why**: Avoid large divergences, minimize conflicts.

**Frequency**:
- After enabling/disabling plugins: immediate
- After modifying skills: immediate
- Periodic check: weekly

**Command**:
```bash
# Check for changes
./sync-config.py status

# If changes exist
./sync-config.py pull
git add claude-config/
git commit -m "Sync configuration changes"
git push
```

### 6. Test on New Machine Before Syncing Back

**Why**: Ensure configuration works before propagating changes.

**Pattern**:
```bash
# On new machine, after push
./sync-config.py push

# Test thoroughly
# ... use Claude Code ...

# Only after confirming it works
./sync-config.py pull  # If you made changes
```

### 7. Keep Global Config Minimal

**Why**: Reduces machine-specific divergence.

**Strategy**:
- Enable only essential plugins globally
- Use project-specific configs for specialized needs
- Document exceptions in planning journal

**Example global settings.json**:
```json
{
  "enabledPlugins": {
    "superpowers@claude-plugins-official": true
  }
}
```

**Example project-specific .claude-project.json**:
```json
{
  "enabledSkills": [
    "scientific-skills:pubmed-database",
    "scientific-skills:biopython"
  ]
}
```

### 8. Backup Before Major Changes

**Why**: Easy rollback if something goes wrong.

**Automatic**: Backups are enabled by default in `sync.config.yaml`

**Manual backup**:
```bash
# Before major change
cp -r ~/.claude ~/.claude.backup.$(date +%Y%m%d)

# After confirming change works
rm -r ~/.claude.backup.YYYYMMDD
```

### 9. Document Machine-Specific Quirks

**Why**: Helps future you and collaborators.

**Where**: Planning journal

**Example**:
```markdown
# Machine-Specific Notes

## macbook-pro-work
- Has VPN that blocks some database access
- Requires proxy settings for external APIs
- Uses different Python version (3.11 vs 3.9)

## linux-server
- Minimal plugins (no GUI-based skills)
- Headless setup (no image generation)
- Network-restricted (no external database access)
```

### 10. Review Planning Journal Periodically

**Why**: Learn from past experiences, improve practices.

**Frequency**: Monthly

**Process**:
```bash
# List recent entries
./sync-config.py plan --list

# Review entries
cat planning/$(hostname)/2026-02-*.md

# Extract lessons learned
# Update this document or team wiki
```

## Advanced Scenarios

### Scenario: Selective Skill Sync

You want to sync only specific skills, not all.

**Approach**: Modify `sync.config.yaml`

```yaml
sync_rules:
  always:
    - path: "settings.json"
      description: "Plugin flags"

    # Sync only specific skill directories
    - path: "skills/superpowers:test-driven-development/"
      description: "TDD skill"

    - path: "skills/scientific-skills:pubmed-database/"
      description: "PubMed database skill"
```

**Note**: This requires modifying the sync script to handle individual skill directories.

### Scenario: Multiple Configuration Profiles

You want different configurations for different contexts (work, personal, research).

**Approach**: Use git branches

```bash
# Create work profile
git checkout -b config/work
./sync-config.py pull  # Pull work config
git commit -m "Work configuration profile"
git push origin config/work

# Create research profile
git checkout -b config/research
# ... modify settings for research ...
./sync-config.py pull
git commit -m "Research configuration profile"
git push origin config/research

# Switch profiles
git checkout config/work
./sync-config.py push  # Apply work config

git checkout config/research
./sync-config.py push  # Apply research config
```

### Scenario: Automated Sync

You want to automatically sync configuration on a schedule.

**Approach**: Use cron or systemd timer

**Cron example** (sync every hour):
```bash
# Add to crontab
0 * * * * cd ~/repos/claude && ./sync-config.py pull && git add claude-config/ && git commit -m "Auto-sync $(date)" && git push
```

**Systemd timer example**:
```ini
# ~/.config/systemd/user/claude-config-sync.service
[Unit]
Description=Claude Code Configuration Sync

[Service]
Type=oneshot
WorkingDirectory=/home/user/repos/claude
ExecStart=/home/user/repos/claude/sync-config.py pull
ExecStartPost=/usr/bin/git add claude-config/
ExecStartPost=/usr/bin/git commit -m "Auto-sync"
ExecStartPost=/usr/bin/git push

# ~/.config/systemd/user/claude-config-sync.timer
[Unit]
Description=Claude Code Configuration Sync Timer

[Timer]
OnCalendar=hourly
Persistent=true

[Install]
WantedBy=timers.target
```

**Enable**:
```bash
systemctl --user enable claude-config-sync.timer
systemctl --user start claude-config-sync.timer
```

### Scenario: CI/CD Integration

You want to validate configuration changes in CI.

**Approach**: GitHub Actions

**.github/workflows/validate-config.yml**:
```yaml
name: Validate Configuration

on:
  pull_request:
    paths:
      - 'claude-config/**'
      - 'project-configs/**'

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'

      - name: Install dependencies
        run: pip install pyyaml

      - name: Validate config
        run: |
          # Check YAML syntax
          python -c "import yaml; yaml.safe_load(open('sync.config.yaml'))"

          # Check JSON syntax
          if [ -f claude-config/settings.json ]; then
            python -c "import json; json.load(open('claude-config/settings.json'))"
          fi

      - name: Check for excluded files
        run: |
          # Ensure plans/ directory is not committed
          if [ -d claude-config/plans/ ]; then
            echo "Error: plans/ directory should not be committed"
            exit 1
          fi

      - name: Dry-run push
        run: ./sync-config.py push --dry-run
```

### Scenario: Configuration Diff Tool

You want to quickly compare configurations across machines.

**Approach**: Create custom script

**config-diff.sh**:
```bash
#!/bin/bash

# Compare configuration between machines

MACHINE_A=$1
MACHINE_B=$2

if [ -z "$MACHINE_A" ] || [ -z "$MACHINE_B" ]; then
    echo "Usage: $0 <machine-a-branch> <machine-b-branch>"
    exit 1
fi

echo "Comparing configurations:"
echo "  Machine A: $MACHINE_A"
echo "  Machine B: $MACHINE_B"
echo

# Compare settings.json
echo "=== settings.json ==="
git diff $MACHINE_A:claude-config/settings.json $MACHINE_B:claude-config/settings.json

# Compare installed plugins
echo "=== installed_plugins.json ==="
git diff $MACHINE_A:claude-config/plugins/installed_plugins.json $MACHINE_B:claude-config/plugins/installed_plugins.json

# List skill differences
echo "=== Skills ==="
comm -3 \
  <(git ls-tree -r --name-only $MACHINE_A:claude-config/skills/ | sort) \
  <(git ls-tree -r --name-only $MACHINE_B:claude-config/skills/ | sort)
```

**Usage**:
```bash
./config-diff.sh config/work config/research
```

## Summary

Key takeaways:

1. **Always preview** with `--dry-run` and `status`
2. **Use planning journal** for ALL configuration changes
3. **Commit atomically** with descriptive messages
4. **Sync regularly** to avoid large divergences
5. **Test before propagating** changes to other machines
6. **Keep global config minimal**, use project-specific configs
7. **Document quirks** in planning journal
8. **Review periodically** to learn and improve

Following these workflows and practices will help you maintain consistent, reproducible Claude Code configuration across all your machines.
