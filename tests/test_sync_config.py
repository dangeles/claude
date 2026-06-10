"""Regression tests for the risky logic in sync-config.py.

Focus areas (the paths where a bug silently corrupts live config):
  - exclusion matching (is_excluded)
  - case normalization for set comparison (normalize_path_for_comparison)
  - the unquoted-colon frontmatter hint (_hint_unquoted_colon_in_value)
  - frontmatter parsing (_parse_frontmatter)
  - the push gates (validate_frontmatter, validate_json)
  - orphan detection (find_orphans) used by --delete
"""
import json


# --------------------------------------------------------------------------
# is_excluded
# --------------------------------------------------------------------------

class TestIsExcluded:
    def test_exact_file_match(self, make_sync, tmp_path):
        sync = make_sync(exclusions=["settings.local.json"])
        base = tmp_path / "base"
        base.mkdir()
        assert sync.is_excluded(base / "settings.local.json", base)
        assert not sync.is_excluded(base / "settings.json", base)

    def test_directory_prefix_match(self, make_sync, tmp_path):
        sync = make_sync(exclusions=["cache/"])
        base = tmp_path / "base"
        base.mkdir()
        assert sync.is_excluded(base / "cache" / "a" / "b.bin", base)
        # A sibling that merely shares the prefix string must NOT be excluded.
        assert not sync.is_excluded(base / "cache-keep" / "b.bin", base)

    def test_glob_workspace_pattern(self, make_sync, tmp_path):
        sync = make_sync(exclusions=["*-workspace/"])
        base = tmp_path / "base"
        base.mkdir()
        assert sync.is_excluded(
            base / "skills" / "patent-review-workspace" / "scratch.md", base
        )
        assert not sync.is_excluded(base / "skills" / "patent-review" / "SKILL.md", base)

    def test_path_outside_base_not_excluded(self, make_sync, tmp_path):
        sync = make_sync(exclusions=["cache/"])
        base = tmp_path / "base"
        base.mkdir()
        outside = tmp_path / "elsewhere" / "cache" / "x"
        assert not sync.is_excluded(outside, base)


# --------------------------------------------------------------------------
# normalize_path_for_comparison
# --------------------------------------------------------------------------

class TestNormalizePath:
    def test_matches_platform(self, make_sync, sync_module):
        from pathlib import Path
        sync = make_sync()
        result = sync.normalize_path_for_comparison(Path("Skills/MySkill.md"))
        if sync_module._IS_DARWIN:
            assert result == "skills/myskill.md"
        else:
            assert result == "Skills/MySkill.md"


# --------------------------------------------------------------------------
# _hint_unquoted_colon_in_value (static)
# --------------------------------------------------------------------------

class TestColonHint:
    def hint(self, sync_module, block):
        return sync_module.ConfigSync._hint_unquoted_colon_in_value(block)

    def test_detects_unquoted_colon(self, sync_module):
        msg = self.hint(sync_module, "description: Use when: something happens")
        assert msg is not None and "description" in msg

    def test_ignores_double_quoted_value(self, sync_module):
        assert self.hint(sync_module, 'description: "Use when: ok"') is None

    def test_ignores_single_quoted_value(self, sync_module):
        assert self.hint(sync_module, "description: 'Use when: ok'") is None

    def test_ignores_block_scalar(self, sync_module):
        assert self.hint(sync_module, "description: |\n  Use when: ok") is None

    def test_ignores_comment_and_list_items(self, sync_module):
        assert self.hint(sync_module, "# a comment: with colon") is None
        assert self.hint(sync_module, "- item: with colon") is None

    def test_clean_value_returns_none(self, sync_module):
        assert self.hint(sync_module, "name: my-skill\ndescription: a clean one") is None


# --------------------------------------------------------------------------
# _parse_frontmatter
# --------------------------------------------------------------------------

class TestParseFrontmatter:
    def write(self, tmp_path, text):
        p = tmp_path / "doc.md"
        p.write_text(text)
        return p

    def test_valid(self, make_sync, tmp_path):
        sync = make_sync()
        p = self.write(tmp_path, "---\nname: x\ndescription: y\n---\nbody\n")
        data, err = sync._parse_frontmatter(p)
        assert err is None
        assert data == {"name": "x", "description": "y"}

    def test_missing_opening(self, make_sync, tmp_path):
        sync = make_sync()
        p = self.write(tmp_path, "no frontmatter here\n")
        data, err = sync._parse_frontmatter(p)
        assert data is None and "missing frontmatter" in err

    def test_unclosed(self, make_sync, tmp_path):
        sync = make_sync()
        p = self.write(tmp_path, "---\nname: x\nbody but no close\n")
        data, err = sync._parse_frontmatter(p)
        assert data is None and "no closing" in err

    def test_not_a_mapping(self, make_sync, tmp_path):
        sync = make_sync()
        p = self.write(tmp_path, "---\n- a\n- b\n---\n")
        data, err = sync._parse_frontmatter(p)
        assert data is None and "not a YAML mapping" in err

    def test_invalid_yaml_surfaces_colon_hint(self, make_sync, tmp_path):
        sync = make_sync()
        # An unquoted colon-space that makes PyYAML choke; the hint should fire.
        p = self.write(tmp_path, "---\ndescription: Use when: x is: y\n---\n")
        data, err = sync._parse_frontmatter(p)
        assert data is None
        assert "invalid YAML" in err


# --------------------------------------------------------------------------
# validate_frontmatter (push gate)
# --------------------------------------------------------------------------

class TestValidateFrontmatter:
    def test_passes_with_good_skill(self, make_sync):
        sync = make_sync()
        skill = sync.target_dir / "skills" / "good" / "SKILL.md"
        skill.parent.mkdir(parents=True)
        skill.write_text("---\nname: good\ndescription: fine\n---\n")
        assert sync.validate_frontmatter() is True

    def test_fails_with_missing_field(self, make_sync):
        sync = make_sync()
        skill = sync.target_dir / "skills" / "bad" / "SKILL.md"
        skill.parent.mkdir(parents=True)
        skill.write_text("---\nname: bad\n---\n")
        assert sync.validate_frontmatter() is False

    def test_fails_with_bad_agent(self, make_sync):
        sync = make_sync()
        agent = sync.target_dir / "agents" / "broken.md"
        agent.parent.mkdir(parents=True)
        agent.write_text("no frontmatter at all\n")
        assert sync.validate_frontmatter() is False


# --------------------------------------------------------------------------
# validate_json (push gate)
# --------------------------------------------------------------------------

class TestValidateJson:
    def test_passes_with_valid_json(self, make_sync):
        sync = make_sync()
        (sync.target_dir / "settings.json").write_text(json.dumps({"a": 1}))
        assert sync.validate_json() is True

    def test_fails_with_malformed_json(self, make_sync):
        sync = make_sync()
        (sync.target_dir / "settings.json").write_text('{"a": 1,}')  # trailing comma
        assert sync.validate_json() is False

    def test_no_json_targets_passes(self, make_sync):
        sync = make_sync()
        # No settings.json written; gate should not fail on absence.
        assert sync.validate_json() is True


# --------------------------------------------------------------------------
# find_orphans (used by --delete)
# --------------------------------------------------------------------------

class TestFindOrphans:
    def _setup(self, sync):
        repo_skills = sync.target_dir / "skills"
        live_skills = sync.source_dir / "skills"
        repo_skills.mkdir(parents=True)
        live_skills.mkdir(parents=True)
        return repo_skills, live_skills

    def test_detects_orphan(self, make_sync):
        sync = make_sync()
        repo, live = self._setup(sync)
        (repo / "keep.md").write_text("x")
        (live / "keep.md").write_text("x")
        (live / "orphan.md").write_text("x")  # in live, not in repo
        orphans = sync.find_orphans("skills/")
        names = {rel.name for _, rel in orphans}
        assert names == {"orphan.md"}

    def test_no_orphans_when_in_sync(self, make_sync):
        sync = make_sync()
        repo, live = self._setup(sync)
        (repo / "a.md").write_text("x")
        (live / "a.md").write_text("x")
        assert sync.find_orphans("skills/") == []

    def test_excluded_files_not_orphaned(self, make_sync):
        # Exclusions in find_orphans are matched relative to the sync-rule base
        # (skills/), so the pattern is "local-only/" not "skills/local-only/".
        sync = make_sync(exclusions=["local-only/"])
        repo, live = self._setup(sync)
        (repo / "a.md").write_text("x")
        (live / "a.md").write_text("x")
        excluded = live / "local-only"
        excluded.mkdir()
        (excluded / "scratch.md").write_text("x")  # excluded -> not an orphan
        assert sync.find_orphans("skills/") == []
