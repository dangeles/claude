# Claude Code Configuration Guide

Comprehensive guide to managing Claude Code configuration across multiple machines using this repository.

## Table of Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Setup](#setup)
- [File Reference](#file-reference)
- [Usage](#usage)
- [Troubleshooting](#troubleshooting)

## Overview

This repository provides a bidirectional sync system for Claude Code configuration, enabling:

1. **Version control** for Claude Code settings, skills, and plugin configurations
2. **Multi-machine sync** for consistent configuration across devices
3. **Project-specific configs** tracked separately per project
4. **Planning journal** to document and assess configuration changes
5. **Safe conflict resolution** with interactive prompts and backups

### What's Tracked

**User-wide configuration** (`~/.claude/`):
- `settings.json` - Plugin enable/disable flags
- `skills/` - ALL skills (including plugin-provided skills)
- `plugins/installed_plugins.json` - Plugin versions for reproducibility

**Project-specific configuration**:
- `.claude-project.json` - Project-specific settings
- `.claude/` - Project-specific skill overrides

### What's Excluded

Machine-specific and ephemeral data:
- `settings.local.json` - Machine-specific permissions
- `history.jsonl` - Session history
- `plans/` - Claude's ephemeral plans (NOT synced)
- `debug/`, `cache/`, `downloads/` - Temporary data
- `projects/`, `tasks/`, `todos/` - Session-specific data
- `plugins/cache/` - Plugin binaries (re-downloaded per machine)

## Architecture

### Directory Structure

```
~/repos/claude/
├── sync-config.py              # Main sync script
├── sync.config.yaml            # Sync configuration
│
├── claude-config/              # User-wide config (synced)
│   ├── settings.json           # Plugin flags
│   ├── skills/                 # ALL skills
│   │   └── [skill-name]/
│   │       ├── SKILL.md
│   │       ├── examples/
│   │       └── references/
│   ├── plugins/
│   │   └── installed_plugins.json
│   └── .backups/               # Automatic backups (not in git)
│
├── project-configs/            # Project-specific configs
│   ├── bioreactor/
│   │   └── .claude-project.json
│   └── [other-projects]/
│
├── planning/                   # Configuration journal
│   ├── [hostname]/
│   │   └── YYYY-MM-DD-description.md
│   ├── README.md
│   └── .template.md
│
└── docs/
    ├── CLAUDE_CONFIG_GUIDE.md  # This file
    └── SYNC_WORKFLOW.md        # Workflow examples
```

### Sync Flow

**Pull operation** (`.claude/` → repo):
```
~/.claude/settings.json → claude-config/settings.json
~/.claude/skills/       → claude-config/skills/
~/.claude/plugins/installed_plugins.json → claude-config/plugins/installed_plugins.json
```

**Push operation** (repo → `.claude/`):
```
claude-config/settings.json → ~/.claude/settings.json
claude-config/skills/       → ~/.claude/skills/
claude-config/plugins/installed_plugins.json → ~/.claude/plugins/installed_plugins.json
```

## Setup

### Initial Setup (Machine A)

1. **Clone or navigate to this repository**:
   ```bash
   cd ~/repos/claude
   ```

2. **Install dependencies**:
   ```bash
   pip install pyyaml
   ```

3. **Verify Claude Code directory exists**:
   ```bash
   ls ~/.claude
   ```

4. **Preview what will be synced**:
   ```bash
   ./sync-config.py pull --dry-run
   ```

5. **Pull current configuration**:
   ```bash
   ./sync-config.py pull
   ```

6. **Review synced files**:
   ```bash
   ls claude-config/
   cat claude-config/settings.json
   ```

7. **Commit to repository**:
   ```bash
   git add claude-config/ planning/
   git commit -m "Initial Claude Code configuration sync"
   git push
   ```

### Setup on Additional Machine (Machine B)

1. **Clone repository**:
   ```bash
   git clone <repo-url> ~/repos/claude
   cd ~/repos/claude
   ```

2. **Install dependencies**:
   ```bash
   pip install pyyaml
   ```

3. **Preview what will be pushed**:
   ```bash
   ./sync-config.py push --dry-run
   ```

4. **Push configuration to ~/.claude/**:
   ```bash
   ./sync-config.py push
   ```

5. **Verify configuration**:
   ```bash
   ls ~/.claude
   cat ~/.claude/settings.json
   ```

6. **Restart Claude Code** (if running) to pick up changes

## File Reference

### User-Wide Configuration Files

#### `settings.json`

Plugin enable/disable flags and global settings.

**Location**: `~/.claude/settings.json`

**Example**:
```json
{
  "enabledPlugins": {
    "superpowers@claude-plugins-official": true,
    "scientific-skills@claude-scientific-skills": true
  },
  "theme": "dark",
  "maxTokens": 200000
}
```

**What to track**: YES - This controls which plugins are active

**Why**: Ensures consistent plugin availability across machines

#### `skills/`

Directory containing all skill definitions, including plugin-provided skills.

**Location**: `~/.claude/skills/`

**Structure**:
```
skills/
├── superpowers:using-superpowers/
│   ├── SKILL.md
│   └── examples/
├── scientific-skills:pubmed-database/
│   ├── SKILL.md
│   └── references/
└── [other-skills]/
```

**What to track**: YES - ALL skills

**Why**: Skills define Claude Code's capabilities. Tracking all skills (including plugin-provided ones) ensures consistent behavior across machines.

**Note**: This directory can be large (~1.1 MB with scientific-skills plugin). Git handles this efficiently with compression.

#### `plugins/installed_plugins.json`

Plugin version manifest for reproducibility.

**Location**: `~/.claude/plugins/installed_plugins.json`

**Example**:
```json
{
  "plugins": [
    {
      "name": "superpowers",
      "version": "1.2.3",
      "author": "claude-plugins-official"
    },
    {
      "name": "scientific-skills",
      "version": "2.0.1",
      "author": "claude-scientific-skills"
    }
  ]
}
```

**What to track**: YES

**Why**: Ensures same plugin versions across machines, avoiding version-related issues

#### `settings.local.json`

Machine-specific settings (file permissions, local paths).

**Location**: `~/.claude/settings.local.json`

**What to track**: NO (excluded)

**Why**: Contains machine-specific permissions that shouldn't be synced

#### `plans/`

Claude's ephemeral planning files created during plan mode.

**Location**: `~/.claude/plans/`

**What to track**: NO (excluded)

**Why**: These are temporary files created during Claude Code sessions. They are session-specific and not meant to be shared across machines.

**Important**: This is different from `planning/` in the repository, which is for YOUR planning journal entries about config changes.

### Project-Specific Configuration

#### `.claude-project.json`

Project-specific Claude Code settings.

**Location**: `[project-root]/.claude-project.json`

**Example**:
```json
{
  "enabledSkills": [
    "scientific-skills:pubmed-database",
    "scientific-skills:uniprot-database",
    "superpowers:test-driven-development"
  ],
  "disabledSkills": [
    "scientific-skills:matlab"
  ],
  "projectContext": "Bioreactor optimization project"
}
```

**What to track**: YES (in `project-configs/`)

**Why**: Project-specific settings should be versioned with the project

#### `.claude/` (project-specific)

Project-specific skill overrides or custom skills.

**Location**: `[project-root]/.claude/`

**What to track**: YES (in `project-configs/`)

**Why**: Custom skills or project-specific modifications should be versioned

## Usage

### Basic Operations

#### Pull Configuration (Source → Repo)

Pull current configuration from `~/.claude/` into the repository:

```bash
./sync-config.py pull
```

This copies:
- `~/.claude/settings.json` → `claude-config/settings.json`
- `~/.claude/skills/` → `claude-config/skills/`
- `~/.claude/plugins/installed_plugins.json` → `claude-config/plugins/installed_plugins.json`

**When to use**:
- After enabling/disabling plugins
- After installing new plugins
- After modifying skills
- Before committing configuration changes

#### Push Configuration (Repo → Source)

Push repository configuration to `~/.claude/`:

```bash
./sync-config.py push
```

This copies:
- `claude-config/settings.json` → `~/.claude/settings.json`
- `claude-config/skills/` → `~/.claude/skills/`
- `claude-config/plugins/installed_plugins.json` → `~/.claude/plugins/installed_plugins.json`

**When to use**:
- Setting up new machine
- Restoring configuration from repository
- Applying changes from another machine

**Safety**: Always prompts for confirmation before modifying `~/.claude/`

#### Check Status

Show differences between source and target:

```bash
./sync-config.py status
```

Output:
```
[INFO] Comparing configuration...

User-wide configuration:
  Source: /Users/user/.claude
  Target: /Users/user/repos/claude/claude-config

Changes detected:
  • settings.json: modified
  • skills/: 3 new, 1 modified
```

**When to use**:
- Before pulling or pushing
- To check if configuration has diverged
- During troubleshooting

#### Preview Changes (Dry Run)

Preview what would happen without making changes:

```bash
./sync-config.py pull --dry-run
./sync-config.py push --dry-run
```

Output:
```
[INFO] Would copy: settings.json -> claude-config/settings.json
[INFO] Would copy: skills/superpowers:test-driven-development/SKILL.md -> ...
```

**When to use**:
- Before any sync operation (recommended)
- When unsure about changes
- For verification

### Project Configuration Operations

#### Pull Project Configs

Discover and sync all project-specific configurations:

```bash
./sync-config.py pull-projects
```

This:
1. Searches `~/repos/*` for projects with `.claude-project.json` or `.claude/`
2. Copies configs to `project-configs/[project-name]/`

Output:
```
[INFO] Discovering projects with Claude Code configuration...
[INFO] Found 2 project(s):
  - bioreactor (/Users/user/repos/bioreactor)
  - protein-design (/Users/user/repos/protein-design)

[INFO] Syncing project: bioreactor
[SUCCESS] Copied: .claude-project.json
```

#### Pull Specific Project

Sync only a specific project:

```bash
./sync-config.py pull-projects bioreactor
```

#### Push Project Configs

Restore project configurations from repository:

```bash
./sync-config.py push-projects
```

This:
1. Finds project configurations in `project-configs/`
2. Locates corresponding projects in search paths
3. Copies configs back to projects

**Note**: If a project doesn't exist on the machine, it's skipped gracefully.

#### Push Specific Project

Restore config for a specific project:

```bash
./sync-config.py push-projects bioreactor
```

### Planning Journal

#### Create Planning Entry

Before making configuration changes:

```bash
./sync-config.py plan --title "Enable scientific skills plugin"
```

This:
1. Creates `planning/[hostname]/YYYY-MM-DD-enable-scientific-skills-plugin.md`
2. Fills template with date and hostname
3. Opens in your editor (`$EDITOR`)

See [planning/README.md](../planning/README.md) for detailed planning workflow.

#### List Planning Entries

Show all planning journal entries:

```bash
./sync-config.py plan --list
```

Filter by machine:

```bash
./sync-config.py plan --list --machine macbook-pro
```

### Conflict Resolution

When files differ between source and target, you'll see an interactive prompt:

```
Conflict detected:
  Source: ~/.claude/settings.json
    Modified: 2026-02-02 18:30, 256 bytes
  Target: ./claude-config/settings.json
    Modified: 2026-02-01 10:00, 243 bytes

Diff preview:
--- target/settings.json
+++ source/settings.json
@@ -2,6 +2,7 @@
   "enabledPlugins": {
     "superpowers@claude-plugins-official": true,
-   "scientific-skills@claude-scientific-skills": false
+   "scientific-skills@claude-scientific-skills": true,
+   "new-plugin@author": true
   }
 }

Choose resolution:
[1] Use source (overwrite target)
[2] Use target (keep current)
[3] Show full diff
[4] Skip this file
[5] Abort sync

Choice:
```

**Options**:
1. **Use source**: Copy source file to target (overwrite)
2. **Use target**: Keep target file unchanged
3. **Show full diff**: Display complete diff (if diff tool available)
4. **Skip this file**: Skip and continue with next file
5. **Abort sync**: Stop entire sync operation

**Backups**: Before overwriting, a backup is automatically created in `claude-config/.backups/`

### Advanced Options

#### Verbose Output

Show detailed information about operations:

```bash
./sync-config.py pull --verbose
```

Output includes:
- Files skipped (identical)
- Files excluded
- Detailed error messages

#### Disable Backups

Skip backup creation (not recommended):

```bash
./sync-config.py push --no-backup
```

#### Custom Config File

Use a different configuration file:

```bash
./sync-config.py pull --config custom.config.yaml
```

## Troubleshooting

### Problem: "Source directory does not exist"

**Error**:
```
[ERROR] Source directory does not exist: /Users/user/.claude
```

**Cause**: Claude Code hasn't been run yet, or configuration directory is in a different location.

**Solution**:
1. Run Claude Code at least once to initialize configuration
2. Verify location: `ls ~/.claude`
3. Check if Claude Code uses a different config directory

### Problem: Conflicts on every pull

**Symptom**: Every pull operation shows conflicts, even for unchanged files.

**Cause**: Line ending differences (CRLF vs LF) or timestamp-only changes.

**Solution**:
1. Ensure consistent line endings:
   ```bash
   git config core.autocrlf input
   ```

2. Check for timestamp issues:
   ```bash
   ./sync-config.py status --verbose
   ```

### Problem: Large skills directory slows down sync

**Symptom**: First sync takes several minutes.

**Cause**: Skills directory is ~1.1 MB with plugins like scientific-skills.

**Solution**:
- This is expected on first sync
- Subsequent syncs only copy changed files (much faster)
- Git compresses this efficiently
- Consider using `.gitattributes` for LFS if needed

### Problem: Plugin binaries not working on new machine

**Symptom**: Plugins fail to load after push.

**Cause**: Plugin binaries (in `~/.claude/plugins/cache/`) are machine-specific and not synced.

**Solution**:
1. Plugins will automatically re-download on first run
2. Wait for download to complete
3. Restart Claude Code

### Problem: Permission denied when pushing

**Error**:
```
[ERROR] Permission denied: /Users/user/.claude/settings.json
```

**Cause**: File permissions issue or Claude Code is running.

**Solution**:
1. Close Claude Code
2. Check file permissions: `ls -la ~/.claude/settings.json`
3. Fix permissions if needed: `chmod 644 ~/.claude/settings.json`

### Problem: Planning entry won't open in editor

**Symptom**: Planning entry created but doesn't open in editor.

**Cause**: `$EDITOR` environment variable not set or editor not in PATH.

**Solution**:
1. Set editor:
   ```bash
   export EDITOR=nano  # or vim, code, etc.
   ```

2. Add to shell profile (`~/.bashrc` or `~/.zshrc`):
   ```bash
   echo 'export EDITOR=nano' >> ~/.zshrc
   ```

3. Manually open file:
   ```bash
   nano planning/$(hostname)/YYYY-MM-DD-description.md
   ```

### Problem: Project configs not discovered

**Symptom**: `pull-projects` finds no projects.

**Cause**: Projects not in configured search paths, or no `.claude-project.json` files.

**Solution**:
1. Check search paths in `sync.config.yaml`:
   ```yaml
   project_configs:
     search_paths:
       - "~/repos/*"
   ```

2. Verify projects have config files:
   ```bash
   find ~/repos -name ".claude-project.json"
   ```

3. Add search paths if needed:
   ```yaml
   project_configs:
     search_paths:
       - "~/repos/*"
       - "~/projects/*"
       - "~/work/*"
   ```

### Problem: Sync aborted accidentally

**Symptom**: Chose "abort" during conflict resolution, want to retry.

**Solution**:
- Backups are already created (if enabled)
- Simply run the sync command again
- Choose different resolution option

### Problem: Lost local changes after push

**Symptom**: Accidentally overwrote local changes with repository version.

**Solution**:
1. Check backups:
   ```bash
   ls claude-config/.backups/
   ```

2. Restore from backup:
   ```bash
   cp claude-config/.backups/settings.json.TIMESTAMP.bak ~/.claude/settings.json
   ```

3. If no backup, check git history:
   ```bash
   git log claude-config/settings.json
   git show COMMIT:claude-config/settings.json
   ```

### Problem: Python dependencies missing

**Error**:
```
ModuleNotFoundError: No module named 'yaml'
```

**Solution**:
```bash
pip install pyyaml
```

Or use a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install pyyaml
```

### Getting Help

If you encounter issues not covered here:

1. Check verbose output: `--verbose` flag
2. Check sync status: `./sync-config.py status`
3. Review configuration: `cat sync.config.yaml`
4. Check file permissions: `ls -la ~/.claude`
5. Verify git status: `git status`

## Summary

This system provides:
- ✅ Version-controlled Claude Code configuration
- ✅ Multi-machine synchronization
- ✅ Project-specific configuration tracking
- ✅ Planning journal for configuration changes
- ✅ Safe conflict resolution with backups
- ✅ Reproducible setup across machines

For workflow examples and best practices, see [SYNC_WORKFLOW.md](./SYNC_WORKFLOW.md).

For planning journal usage, see [planning/README.md](../planning/README.md).
