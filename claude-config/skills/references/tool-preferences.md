# Tool Preferences: Package and Environment Management

This document defines the preferred tools for package installation and environment
management across all skills. All skills MUST follow these preferences when
generating installation commands, environment setup instructions, or dependency
management workflows.

---

## Primary Tool: micromamba

**micromamba** is the preferred package and environment manager for all non-Docker,
non-GUI contexts.

### When to use micromamba

- Installing packages from conda-forge (the primary channel)
- Creating and managing Python environments
- Installing scientific computing packages (numpy, scipy, pandas, FEniCSx, etc.)
- Installing CLI tools (yq, jq, jupytext, etc.)
- Any `conda` command equivalent -- use `micromamba` instead

### Commands

```bash
# Environment management:
micromamba create -n ENV_NAME python=VERSION
micromamba activate ENV_NAME
micromamba deactivate
micromamba env create -f environment.yml
micromamba env export > environment.yml
micromamba env remove -n ENV_NAME
micromamba env list

# Package management:
micromamba install PACKAGE
micromamba install -c conda-forge PACKAGE
micromamba search PACKAGE
```

### Claude Code Bash Sessions

> **Note**: Claude Code's Bash tool does not inherit the user's shell initialization files (`.zshrc`, `.bashrc`, etc.). A bare `micromamba activate ENV_NAME` will fail silently in agent-executed Bash sessions.

**Canonical pattern for Claude Code Bash tool sessions:**

```bash
# Initialize shell hook + activate + run commands — all on one line:
eval "$(micromamba shell hook --shell=zsh)" 2>/dev/null; micromamba activate ENV_NAME 2>/dev/null; YOUR_COMMANDS_HERE
```

**Key details:**
- `--shell=zsh`: Matches macOS/Darwin with zsh. Change to `--shell=bash` for bash users.
- `2>/dev/null` on both: Suppresses harmless warnings; graceful degradation if micromamba is absent.
- Semicolons (not `&&`): The `eval` may print warnings but still succeed; `;` is more forgiving.
- Everything on one line: The Bash tool does not persist shell state between calls.

**When required:**
- Any time a skill or agent runs `micromamba activate` via the Claude Code Bash tool.
- When adapting multi-step human-facing bash blocks for agent execution.

**When NOT required:**
- Human-facing setup instructions run in the user's terminal (shell already initialized).
- Docker contexts (micromamba typically not installed).
- Virtualenv/venv workflows (no micromamba).

---

### Channel configuration

- **Primary channel**: conda-forge
- **Bioinformatics channel**: bioconda (when needed)
- **Channel priority**: strict

Note: `conda-forge` is the channel name. It is NOT a reference to the `conda` tool.
Never change `conda-forge` to `micromamba-forge`.

---

## Secondary Tool: pip

**pip** is acceptable in the following contexts:

### When to use pip

- PyPI-only packages not available on conda-forge (e.g., scanpy, pydeseq2)
- Inside virtualenv/venv workflows (where micromamba is not the environment manager)
- `pip install -r requirements.txt` for generic requirements files
- `pip freeze > requirements.txt` for exporting pip-installed packages
- `!pip install` in Jupyter notebooks (magic commands)
- Inside Docker containers (where micromamba is not installed)

### When NOT to use pip

- For packages available on conda-forge (use micromamba instead)
- As the primary environment manager (use micromamba instead)
- For packages like numpy, pandas, jupyter, matplotlib, ipykernel, pyvista, psutil
  (all available on conda-forge)

---

## Acceptable: brew (macOS only, GUI/cask installs only)

**brew** (Homebrew) is acceptable ONLY for macOS GUI application installs via cask:

```bash
# Acceptable:
brew install --cask mactex

# NOT acceptable (use micromamba instead):
# brew install yq
# brew install jq
```

### Detection paths

`/opt/homebrew/bin/` paths in detection/PATH checking code are acceptable and should
NOT be changed. These are system detection paths, not installation instructions.

---

## Docker contexts

Inside Dockerfiles and Docker-related code blocks, use the appropriate tools for
the container environment:

- `apt-get install` for system packages
- `pip install` for Python packages
- Do NOT suggest micromamba inside Docker containers unless the Dockerfile
  explicitly installs it

---

## Summary Table

| Context | Preferred Tool | Acceptable Alternatives |
|---------|---------------|------------------------|
| conda-forge packages | micromamba | -- |
| PyPI-only packages | pip | -- |
| Environment creation | micromamba | -- |
| macOS GUI apps | brew --cask | -- |
| Docker containers | apt-get + pip | -- |
| Virtualenv workflows | pip | -- |
| Jupyter magic commands | !pip | -- |
| Generic requirements.txt | pip | -- |
