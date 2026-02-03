# Claude Code Repository

A repository for Claude Code projects and experiments, with a comprehensive configuration management system for syncing Claude Code settings, skills, and plugins across multiple machines.

## Quick Start

### Using the Configuration Sync Tool

```bash
# Pull current configuration from ~/.claude to repository
./sync-config.py pull

# Push configuration from repository to ~/.claude
./sync-config.py push

# Check status (show differences)
./sync-config.py status

# Sync project-specific configs
./sync-config.py pull-projects
./sync-config.py push-projects

# Create planning journal entry
./sync-config.py plan --title "Description of config change"
```

## Configuration Management

This repository includes a bidirectional sync system for managing Claude Code configuration across multiple machines.

### What's Tracked

**User-wide configuration** (`~/.claude/`):
- ✅ `settings.json` - Plugin enable/disable flags
- ✅ `skills/` - ALL skills (including plugin-provided)
- ✅ `plugins/installed_plugins.json` - Plugin versions

**Project-specific configuration**:
- ✅ `.claude-project.json` - Project-specific settings
- ✅ `.claude/` - Project-specific skill overrides

**Planning journal**:
- ✅ Configuration change planning and assessment
- ✅ Organized by machine hostname
- ✅ Linked to git commits

### What's Excluded

Machine-specific and ephemeral data:
- ❌ `settings.local.json` - Machine-specific permissions
- ❌ `history.jsonl` - Session history
- ❌ `plans/` - Claude's ephemeral plans (session-specific)
- ❌ `debug/`, `cache/`, `downloads/` - Temporary data
- ❌ `plugins/cache/` - Plugin binaries (re-downloaded per machine)

### Features

- **Bidirectional sync** between `~/.claude/` and this repository
- **Project-specific configs** tracked separately per project
- **Planning journal** to document and assess configuration changes
- **Conflict resolution** with interactive prompts and diff preview
- **Automatic backups** before overwriting files
- **Dry-run mode** to preview changes before executing
- **Multi-machine support** with machine-specific planning entries

### Setup

1. **Install dependencies**:
   ```bash
   pip install pyyaml
   ```

2. **Pull current configuration**:
   ```bash
   ./sync-config.py pull
   ```

3. **Commit to repository**:
   ```bash
   git add claude-config/ planning/
   git commit -m "Initial Claude Code configuration sync"
   git push
   ```

4. **On other machines**:
   ```bash
   git clone <repo-url> ~/repos/claude
   cd ~/repos/claude
   ./sync-config.py push
   ```

### Common Commands

```bash
# Preview changes (recommended before any operation)
./sync-config.py pull --dry-run
./sync-config.py push --dry-run

# Show differences between source and target
./sync-config.py status

# Sync with verbose output
./sync-config.py pull --verbose

# Pull specific project config
./sync-config.py pull-projects bioreactor

# List planning journal entries
./sync-config.py plan --list
./sync-config.py plan --list --machine macbook-pro
```

## Documentation

- **[Configuration Guide](docs/CLAUDE_CONFIG_GUIDE.md)** - Comprehensive setup and usage guide
- **[Sync Workflow](docs/SYNC_WORKFLOW.md)** - Workflow examples and best practices
- **[Planning Journal Guide](planning/README.md)** - Using the planning journal system

## Repository Structure

```
~/repos/claude/
├── README.md                       # This file
├── .gitignore                      # Git exclusions
│
├── sync-config.py                  # Main sync script
├── sync.config.yaml                # Sync configuration
│
├── claude-config/                  # User-wide configuration (synced)
│   ├── settings.json               # Plugin enable/disable flags
│   ├── skills/                     # ALL skills
│   └── plugins/
│       └── installed_plugins.json  # Plugin versions
│
├── project-configs/                # Project-specific configs
│   └── [project-name]/
│       └── .claude-project.json
│
├── planning/                       # Configuration journal
│   ├── [hostname]/
│   │   └── YYYY-MM-DD-description.md
│   ├── README.md
│   └── .template.md
│
└── docs/
    ├── CLAUDE_CONFIG_GUIDE.md      # Comprehensive guide
    └── SYNC_WORKFLOW.md            # Workflow examples
```

## Workflows

### Enable a New Plugin

```bash
# 1. Create planning entry
./sync-config.py plan --title "Enable scientific-skills plugin"

# 2. Enable plugin in Claude Code

# 3. Pull configuration
./sync-config.py pull

# 4. Commit and push
git add claude-config/ planning/
git commit -m "Enable scientific-skills plugin"
git push
```

### Set Up New Machine

```bash
# 1. Clone repository
git clone <repo-url> ~/repos/claude
cd ~/repos/claude

# 2. Push configuration to ~/.claude
./sync-config.py push

# 3. Restart Claude Code
# Plugins will download automatically
```

### Sync Project Configuration

```bash
# Pull from project
./sync-config.py pull-projects bioreactor

# Commit
git add project-configs/
git commit -m "Sync bioreactor project config"
git push

# On another machine
git pull
./sync-config.py push-projects bioreactor
```

## Best Practices

1. **Always preview** with `--dry-run` before executing
2. **Use planning journal** for all configuration changes
3. **Commit atomically** with descriptive messages
4. **Sync regularly** to avoid large divergences
5. **Test on new machine** before syncing back
6. **Keep global config minimal**, use project-specific configs

## Troubleshooting

See [CLAUDE_CONFIG_GUIDE.md](docs/CLAUDE_CONFIG_GUIDE.md#troubleshooting) for detailed troubleshooting guide.

Common issues:
- **"Source directory does not exist"**: Run Claude Code at least once to initialize `~/.claude/`
- **Plugin binaries not working**: Plugins auto-download on first use (may take a few minutes)
- **Permission denied**: Close Claude Code before pushing configuration
- **Conflicts on every pull**: Check git line ending settings (`git config core.autocrlf input`)

## License

MIT
