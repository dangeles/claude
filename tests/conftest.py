"""Shared pytest fixtures for the sync-config.py test suite.

sync-config.py and the hook scripts have hyphenated filenames that are not
importable as normal modules, so we load them by path via importlib.
"""
import importlib.util
import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent


def load_module_by_path(name: str, path: Path):
    """Load a module from an arbitrary file path (handles hyphenated names)."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader, f"could not create import spec for {path}"
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def sync_module():
    """The sync-config.py module, loaded once per session."""
    return load_module_by_path("sync_config", REPO_ROOT / "sync-config.py")


@pytest.fixture
def make_sync(sync_module, tmp_path):
    """Factory building an isolated ConfigSync (dry_run) backed by tmp dirs.

    Returns a callable: make_sync(exclusions=[...], source=Path, target=Path).
    Defaults create empty source/target dirs under tmp_path so destructive
    paths never touch the real ~/.claude or the repo.
    """
    def _make(exclusions=None, source=None, target=None):
        source = source or (tmp_path / "live")
        target = target or (tmp_path / "repo" / "claude-config")
        source.mkdir(parents=True, exist_ok=True)
        target.mkdir(parents=True, exist_ok=True)

        config = {
            "version": "1.0",
            "source_dir": str(source),
            "target_dir": "./claude-config",
            "sync_rules": {
                "always": [
                    {"path": "settings.json", "description": "flags"},
                    {"path": "skills/", "description": "skills"},
                    {"path": "agents/", "description": "agents"},
                ]
            },
            "exclusions": exclusions if exclusions is not None else [],
            "conflict_resolution": {"default": "prompt"},
            "backup": {"enabled": False, "location": "./.backups"},
        }
        # repo_root is config_path.parent, so target_dir resolves under it.
        config_path = target.parent / "sync.config.yaml"
        config_path.write_text(yaml.safe_dump(config))
        return sync_module.ConfigSync(str(config_path), dry_run=True)

    return _make


FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "frontmatter"


@pytest.fixture(scope="session")
def frontmatter_fixtures():
    """Shared good/bad frontmatter samples used by BOTH the sync-config gate
    and the skill-frontmatter-validator hook, so the two implementations are
    held to the same contract.

    Returns list of (path, is_valid) tuples.
    """
    cases = []
    for f in sorted(FIXTURE_DIR.glob("*.md")):
        cases.append((f, f.name.startswith("good_")))
    return cases
