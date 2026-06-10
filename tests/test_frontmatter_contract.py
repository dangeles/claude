"""Contract test pinning the two frontmatter validators together.

The frontmatter rule is implemented twice, on purpose:
  - sync-config.py validate_frontmatter / _parse_frontmatter (push gate, uses PyYAML)
  - claude-config/hooks/skill-frontmatter-validator.py _validate_frontmatter
    (PreToolUse hook, deliberately zero-dep so it stays fast/portable)

Two implementations drift. This test feeds BOTH the same shared fixture set
(tests/fixtures/frontmatter/good_*.md vs bad_*.md) and asserts they agree on
the valid/invalid verdict, so a change to one that diverges from the other
fails CI.
"""
from pathlib import Path

from conftest import REPO_ROOT, load_module_by_path

HOOK = load_module_by_path(
    "skill_frontmatter_validator",
    REPO_ROOT / "claude-config" / "hooks" / "skill-frontmatter-validator.py",
)


def _sync_verdict(sync_module, path: Path) -> bool:
    """Reproduce the sync-config push gate's verdict for a single file."""
    sync = sync_module.ConfigSync.__new__(sync_module.ConfigSync)
    sync.verbose = False
    data, err = sync._parse_frontmatter(path)
    if err is not None or data is None:
        return False
    required = sync_module.REQUIRED_FRONTMATTER_FIELDS
    return all(data.get(k) not in (None, "") for k in required)


def _hook_verdict(path: Path) -> bool:
    ok, _ = HOOK._validate_frontmatter(path.read_text())
    return ok


def test_both_validators_agree_on_fixtures(sync_module, frontmatter_fixtures):
    assert frontmatter_fixtures, "no frontmatter fixtures discovered"
    mismatches = []
    for path, expected_valid in frontmatter_fixtures:
        sync_ok = _sync_verdict(sync_module, path)
        hook_ok = _hook_verdict(path)
        if not (sync_ok == hook_ok == expected_valid):
            mismatches.append(
                f"{path.name}: expected={expected_valid} "
                f"sync={sync_ok} hook={hook_ok}"
            )
    assert not mismatches, "validator disagreement:\n" + "\n".join(mismatches)
