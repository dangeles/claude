# Replace brew/conda/pip with micromamba across all skills

**Date**: 2026-02-21
**Machine**: mac
**Status**: Complete

## Objective

Replace all brew install, conda commands, and pip install (for conda-forge packages) with micromamba equivalents across all skill files. User's preferred package/environment manager is micromamba, not Homebrew or Anaconda/Miniconda.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow (via skill-editor 4-phase process)
- [x] Replace brew install -> micromamba install (except cask installs)
- [x] Replace all conda commands -> micromamba equivalents
- [x] Replace pip install -> micromamba install for conda-forge packages
- [x] Update prose, headers, and comments to reference micromamba
- [x] Rewrite Miniconda WSL2 installation section to micromamba
- [x] Handle env update/export behavioral differences with notes
- [x] Create tool-preferences.md reference file for future drift prevention
- [x] Sync to ~/.claude/

## Expected Outcome

All skill workflows use micromamba commands. Future skills can reference tool-preferences.md for correct tool choices.

## Actual Outcome

21 files modified across 9 skills, 1 new reference file created. 245 insertions, 199 deletions. All validation checks passed (YAML, no conda-forge corruption, cask installs preserved, Docker contexts untouched).

## Assessment

**Result**: Success

**Improvements**:
- All workflows now reference micromamba consistently
- conda env export includes compound command with pip freeze (addresses micromamba limitation)
- conda env update shows both incremental and full recreation options
- New tool-preferences.md prevents future drift

**Issues**:
- Pre-existing sync conflicts on settings.json and plugins (unrelated, skipped)

**Lessons Learned**:
- micromamba env update and env export have behavioral differences from conda -- need notes
- pip extras syntax (e.g., pyvista[jupyter]) not supported by micromamba
- shell init required for micromamba activate

## Related Commits

- [pending]: refactor(skills): replace brew/conda/pip with micromamba across 9 skills

## Next Steps

- Resolve pre-existing settings.json sync conflict
- Monitor for drift in new skills (tool-preferences.md should help)
