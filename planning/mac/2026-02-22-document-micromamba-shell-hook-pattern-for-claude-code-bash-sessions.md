# Document micromamba shell hook pattern for Claude Code Bash sessions

**Date**: 2026-02-22
**Machine**: mac
**Status**: Complete

## Objective

Claude Code's Bash tool does not inherit the user's shell initialization files (.zshrc, .bashrc), causing bare `micromamba activate ENV_NAME` calls to fail silently. Document the canonical workaround pattern (`eval "$(micromamba shell hook --shell=zsh)"`) in the shared tool-preferences reference and update notebook-debugger's environment-management reference where agents execute micromamba activate.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add "Claude Code Bash Sessions" subsection to tool-preferences.md with canonical shell hook pattern
- [x] Add "Claude Code Agent Usage" section to notebook-debugger environment-management.md
- [x] Update 3 agent-executable bash blocks in environment-management.md troubleshooting/registering sections

## Expected Outcome

Skills that need to activate micromamba environments in Claude Code sessions will have a documented, working pattern. Reduces debugging time for environment-related failures. Single source of truth in the shared reference file.

## Actual Outcome

All changes implemented and synced. Five edits across 2 files. All 53 skills pass YAML validation smoke test. Live files verified identical to repo copies.

## Assessment

**Result**: Success

**Improvements**:
- Canonical shell hook pattern now documented in shared reference (tool-preferences.md)
- notebook-debugger environment-management reference updated with 4 shell hook instances
- Human-facing vs agent-executable instructions clearly distinguished

**Issues**:
- Pre-existing sync divergence in plugins/installed_plugins.json (unrelated, auto-resolved)
- Pre-existing 1 new/1 deleted skill difference in sync status (unrelated, not from this change)

**Lessons Learned**:
- sync-config.py --dry-run still prompts interactively on conflicts; use --yes for non-interactive operation

## Related Commits

- [pending]: docs(references): document micromamba shell hook pattern for Claude Code Bash sessions

## Next Steps

- Consider updating cfd-bioreactor references if agents there also execute bare micromamba activate (currently out of scope as those are human-facing setup instructions)
- Monitor for any other skills that encounter silent micromamba activate failures
