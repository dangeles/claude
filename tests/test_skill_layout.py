"""Layout tests for claude-config/skills/.

These guard properties that a case-insensitive filesystem hides. macOS sets
core.ignorecase=true, so a skill committed as `skill.md` still resolves to
`SKILL.md` locally and looks fine — while on Linux (CI, or any fresh clone on a
case-sensitive filesystem) the same skill is invisible to both the
sync-config.py push gate, which globs `*/SKILL.md`, and to skill discovery.

notebook-writer shipped that way from the repo's first commit and went unnoticed
until 2026-07-26; CI was validating 54 of 55 skills.

The checks run against `git ls-files` rather than the working tree, because the
working tree is exactly where the case discrepancy is invisible.
"""
import subprocess
from pathlib import Path

from conftest import REPO_ROOT

SKILLS_PREFIX = "claude-config/skills/"
# Files that live directly in skills/ without being a skill of their own.
NON_SKILL_ENTRIES = {"CHANGELOG.md", "README.md", "references"}


def _tracked_paths() -> list[str]:
    out = subprocess.run(
        ["git", "ls-files", "-z", SKILLS_PREFIX],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return [p for p in out.split("\0") if p]


def _skill_dirs(tracked: list[str]) -> set[str]:
    dirs = set()
    for path in tracked:
        rest = path[len(SKILLS_PREFIX):]
        head = rest.split("/", 1)[0]
        if head in NON_SKILL_ENTRIES or "/" not in rest:
            continue
        dirs.add(head)
    return dirs


def test_every_skill_has_uppercase_skill_md():
    """A skill whose entry point is not exactly SKILL.md is invisible on Linux."""
    tracked = _tracked_paths()
    assert tracked, "git ls-files returned nothing for " + SKILLS_PREFIX
    entry_points = {p for p in tracked if p.rsplit("/", 1)[-1].lower() == "skill.md"}

    wrong_case = sorted(p for p in entry_points if not p.endswith("/SKILL.md"))
    assert not wrong_case, (
        "skill entry point(s) tracked with wrong filename case; these are silently "
        "skipped on case-sensitive filesystems: " + ", ".join(wrong_case)
    )

    have = {p[len(SKILLS_PREFIX):].split("/", 1)[0] for p in entry_points}
    missing = sorted(_skill_dirs(tracked) - have)
    assert not missing, "skill director(ies) with no tracked SKILL.md: " + ", ".join(missing)


def test_no_paths_differing_only_by_case():
    """Two paths differing only by case collide on a case-insensitive checkout."""
    by_lower: dict[str, list[str]] = {}
    for path in _tracked_paths():
        by_lower.setdefault(path.lower(), []).append(path)
    collisions = {k: v for k, v in by_lower.items() if len(set(v)) > 1}
    assert not collisions, "tracked paths differing only by case: " + repr(collisions)


def test_tracked_skill_entry_points_exist_on_disk():
    """Catches an index/worktree mismatch that core.ignorecase would otherwise mask."""
    missing = [
        p
        for p in _tracked_paths()
        if p.endswith("/SKILL.md") and not (Path(REPO_ROOT) / p).is_file()
    ]
    assert not missing, "tracked SKILL.md not present on disk: " + ", ".join(missing)
