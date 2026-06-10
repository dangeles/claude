# Add test suite, CI, and dependency hardening for sync-config.py

**Date**: 2026-06-09
**Machine**: mac
**Status**: Success

## Objective

`sync-config.py` (1278 lines) is the single point of failure for this repo —
it copies, deletes (`--delete`), and gates pushes of live `~/.claude/` config.
It had zero automated tests, no CI, and no pinned dependencies. A regression in
its exclusion matching, orphan deletion, or frontmatter gate could silently
corrupt config across machines. This change adds a safety net.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add `validate_json` gate (settings.json + installed_plugins.json) so the
      JSON-syntax gate documented in CLAUDE.md is actually enforced in code,
      wired into both the `validate` command and `push_config`
- [x] Add `tests/` pytest suite covering the risky pure functions:
      `is_excluded`, `_parse_frontmatter`, `_hint_unquoted_colon_in_value`,
      `normalize_path_for_comparison`, `validate_frontmatter`, `validate_json`,
      and `find_orphans` orphan-detection safety
- [x] Add `requirements.txt` pinning pyyaml + pytest
- [x] Add `.github/workflows/ci.yml` running `validate` + pytest on PRs/pushes
- [x] Cross-reference the duplicated frontmatter logic (sync-config.py gate vs.
      the zero-dep `skill-frontmatter-validator.py` hook) and give both a shared
      fixture set they must agree on (tests/test_frontmatter_contract.py)
- [x] Add `.pre-commit-config.yaml` running `validate` + pytest

## Expected Outcome

The load-bearing sync tooling gains regression coverage; the documented
JSON/frontmatter gates run automatically at commit time (pre-commit) and on PRs
(CI), not just manually at push time. New-machine setup becomes reproducible via
pinned deps.

## Actual Outcome

All seven items landed. 26 pytest tests pass; `./sync-config.py validate` now
reports "Validation passed (frontmatter + JSON)". The only synced change (a
docstring/cross-reference in `skill-frontmatter-validator.py`) was pushed to
`~/.claude/` via `push --yes`; everything else (tests/, CI, requirements.txt,
pre-commit) is repo-only tooling and not synced.

Notable: writing the `find_orphans` test surfaced that exclusions there are
matched relative to the sync-rule base (e.g. `local-only/`), not the live root
(`skills/local-only/`). Behavior is correct; the test now documents it.

## Assessment

**Result**: Success

**Improvements**:
- `sync-config.py` now has a 26-test regression net over its destructive paths
  (exclusion matching, orphan detection, both push gates).
- The CLAUDE.md-documented "JSON syntax gate" is now actually enforced in code,
  not just a manual `jq` step.
- The two frontmatter validators (push gate + zero-dep hook) are pinned together
  by a contract test over shared fixtures, so they can't silently diverge.
- New-machine setup is reproducible via pinned `requirements.txt`.
- Gates now run at three layers: pre-commit, CI, and push.

**Issues**:
- None blocking. CI's `validate` step relies on `source_dir` (~/.claude) not
  being required for the read-only validation — confirmed the constructor only
  creates `target_dir`, so it runs cleanly on a fresh runner.

**Lessons Learned**:
- The hyphenated `sync-config.py` filename forces importlib-by-path loading in
  tests; conftest centralizes that so future test files stay clean.

## Related Commits

- e60931e: test(sync): add test suite, CI, and dependency hardening for sync-config.py

## Next Steps

- Optional: `pre-commit install` locally to activate the commit-time gate.
- Optional: expand tests to cover copy_tree/backup rotation if those paths churn.
