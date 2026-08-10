"""Behavioral tests for the session-scope-guard PreToolUse hook.

The guard blocks staging files Claude didn't touch this session. It is the one
hook here that blocks rather than warns, and it shells out to git, so the cost of
a false positive is high: it stops a legitimate commit and the user has to
disable a safety check to get past it.

Two regressions are pinned here, both observed in real use:

  - `git add some-dir/` was refused when the directory held untracked files.
    `git status --porcelain` collapses an untracked directory to a single
    `?? some-dir/` line, so expansion produced the directory name, which is
    never in the session set. Fixed with `-uall`.
  - The session set accumulated non-paths (`=1`, `/dev/null`) because the
    redirect regex matched comparisons inside quoted strings. The set is printed
    in the block message, so this noise was user-visible.

Also pinned: the guard must still block what it is for, and must fail OPEN
whenever the session set cannot be determined.
"""
import json
import subprocess

import pytest

from conftest import REPO_ROOT, load_module_by_path

HOOK = load_module_by_path(
    "session_scope_guard",
    REPO_ROOT / "claude-config" / "hooks" / "session-scope-guard.py",
)


def _git(repo, *args):
    subprocess.check_call(
        ["git", "-C", str(repo), *args],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A real git repo with one committed file, as the process cwd.

    The hook calls `git rev-parse` and `git status` against the cwd, so these
    tests need genuine git rather than a mock — the collapsed-untracked-directory
    behaviour that caused the bug is git's, and a mock would have hidden it.

    Several tests request this fixture only for its chdir side effect and never
    reference it; that is deliberate, not a leftover parameter.
    """
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    _git(r, "config", "user.email", "t@example.com")
    _git(r, "config", "user.name", "tester")
    (r / "seed.txt").write_text("seed\n")
    _git(r, "add", "seed.txt")
    _git(r, "commit", "-qm", "seed")
    monkeypatch.chdir(r)
    return r


def _transcript(tmp_path, paths, name="t.jsonl"):
    """A transcript whose Write events name `paths`."""
    p = tmp_path / name
    p.write_text("".join(
        json.dumps({"message": {"content": [
            {"type": "tool_use", "name": "Write", "input": {"file_path": x}}
        ]}}) + "\n" for x in paths
    ))
    return str(p)


def _bash_transcript(tmp_path, commands, name="b.jsonl"):
    """A transcript whose Bash events run `commands`."""
    p = tmp_path / name
    p.write_text("".join(
        json.dumps({"message": {"content": [
            {"type": "tool_use", "name": "Bash", "input": {"command": c}}
        ]}}) + "\n" for c in commands
    ))
    return str(p)


def _run(command, transcript, capsys, monkeypatch, description=""):
    payload = {
        "tool_name": "Bash",
        "tool_input": {"command": command, "description": description},
        "transcript_path": transcript,
    }
    text = json.dumps(payload)
    monkeypatch.delenv("SESSION_SCOPE_BYPASS", raising=False)
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
    code = HOOK.main()
    return code, capsys.readouterr().err


# ----- regression 1: directory arguments -----


def test_expand_untracked_directory_lists_files(repo):
    (repo / "tools").mkdir()
    (repo / "tools" / "x.py").write_text("x\n")
    assert HOOK.expand_path_for_session_check("tools/", str(repo)) == ["tools/x.py"]


def test_expand_untracked_directory_without_slash(repo):
    (repo / "tools").mkdir()
    (repo / "tools" / "x.py").write_text("x\n")
    assert HOOK.expand_path_for_session_check("tools", str(repo)) == ["tools/x.py"]


def test_expand_lists_multiple_untracked_files(repo):
    (repo / "pkg").mkdir()
    (repo / "pkg" / "a.py").write_text("a\n")
    (repo / "pkg" / "b.py").write_text("b\n")
    assert sorted(HOOK.expand_path_for_session_check("pkg/", str(repo))) == [
        "pkg/a.py", "pkg/b.py",
    ]


def test_expand_modified_tracked_file_in_directory(repo):
    (repo / "src").mkdir()
    (repo / "src" / "m.py").write_text("one\n")
    _git(repo, "add", "src/m.py")
    _git(repo, "commit", "-qm", "add m")
    (repo / "src" / "m.py").write_text("two\n")
    assert HOOK.expand_path_for_session_check("src/", str(repo)) == ["src/m.py"]


def test_expand_clean_directory_yields_nothing(repo):
    (repo / "clean").mkdir()
    (repo / "clean" / "c.py").write_text("c\n")
    _git(repo, "add", "clean/c.py")
    _git(repo, "commit", "-qm", "add c")
    assert HOOK.expand_path_for_session_check("clean/", str(repo)) == []


def test_expand_leaves_plain_file_alone(repo):
    (repo / "f.py").write_text("f\n")
    assert HOOK.expand_path_for_session_check("f.py", str(repo)) == ["f.py"]


def test_git_add_directory_of_untracked_session_files_allowed(
    repo, tmp_path, capsys, monkeypatch
):
    """The exact failure seen in practice: `git add tools/` with a new file."""
    (repo / "tools").mkdir()
    (repo / "tools" / "token-report.py").write_text("x\n")
    t = _transcript(tmp_path, ["tools/token-report.py"])
    code, err = _run("git add tools/", t, capsys, monkeypatch)
    assert code == 0, err


def test_git_add_directory_still_blocks_foreign_file(repo, tmp_path, capsys, monkeypatch):
    """Expansion must not become a loophole: a non-session file under the
    directory has to keep blocking."""
    (repo / "tools").mkdir()
    (repo / "tools" / "mine.py").write_text("mine\n")
    (repo / "tools" / "theirs.py").write_text("theirs\n")
    t = _transcript(tmp_path, ["tools/mine.py"])
    code, err = _run("git add tools/", t, capsys, monkeypatch)
    assert code == 2
    assert "tools/theirs.py" in err


def test_git_add_clean_directory_is_allowed(repo, tmp_path, capsys, monkeypatch):
    (repo / "clean").mkdir()
    (repo / "clean" / "c.py").write_text("c\n")
    _git(repo, "add", "clean/c.py")
    _git(repo, "commit", "-qm", "add c")
    t = _transcript(tmp_path, ["other.py"])
    code, err = _run("git add clean/", t, capsys, monkeypatch)
    assert code == 0, err


# ----- regression 2: session-set noise -----


@pytest.mark.parametrize("path", [
    "/dev/null", "/dev/stdout", "/dev/stderr", "/dev/fd/2",
    "=1", "&2", "<in", ">out", "",
])
def test_is_stageable_rejects_non_paths(path):
    assert not HOOK.is_stageable(path, "/repo")


def test_is_stageable_rejects_outside_repo():
    assert not HOOK.is_stageable("/tmp/scratch.py", "/repo")


def test_is_stageable_accepts_relative_paths():
    assert HOOK.is_stageable("claude-config/settings.json", "/repo")
    assert HOOK.is_stageable("f.py", "/repo")


def test_is_stageable_keeps_absolutes_without_repo_root():
    assert HOOK.is_stageable("/tmp/scratch.py", None)


def test_comparison_in_quoted_string_is_not_a_redirect():
    """`echo "count >=1"` produced a path named `=1` before the fix."""
    assert HOOK._extract_bash_mutation_paths('echo "sessions with >=1 skill"') == []


def test_fd_duplication_is_not_a_redirect():
    assert HOOK._extract_bash_mutation_paths("cmd >&2") == []
    assert HOOK._extract_bash_mutation_paths("cmd 2>&1") == []


def test_real_redirects_are_still_detected():
    assert "out.txt" in HOOK._extract_bash_mutation_paths("cmd > out.txt")
    assert "log.txt" in HOOK._extract_bash_mutation_paths("cmd >> log.txt")


def test_devnull_is_filtered_from_session_set(repo, tmp_path):
    t = _bash_transcript(tmp_path, ["jq -e . f.json >/dev/null", "echo hi > real.txt"])
    files = HOOK.derive_session_files(t, str(repo))
    assert "real.txt" in files
    assert not any("/dev/" in f for f in files)


def test_session_set_excludes_out_of_repo_writes(repo, tmp_path):
    t = _transcript(tmp_path, ["/tmp/helper.py", "in_repo.py"])
    files = HOOK.derive_session_files(t, str(repo))
    assert files == {"in_repo.py"}


def test_block_message_lists_only_stageable_paths(repo, tmp_path, capsys, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    t = _bash_transcript(tmp_path, ["echo hi > mine.py", "ls >/dev/null"])
    code, err = _run("git add foreign.py", t, capsys, monkeypatch)
    assert code == 2
    assert "/dev/null" not in err
    assert "mine.py" in err


# ----- preserved behaviour -----


def test_bulk_stage_is_refused(repo, tmp_path, capsys, monkeypatch):
    t = _transcript(tmp_path, ["a.py"])
    for cmd in ("git add -A", "git add .", "git add --all", "git add -u"):
        code, err = _run(cmd, t, capsys, monkeypatch)
        assert code == 2, cmd
        assert "bulk-stage" in err


def test_out_of_session_file_is_refused(repo, tmp_path, capsys, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    t = _transcript(tmp_path, ["mine.py"])
    code, err = _run("git add foreign.py", t, capsys, monkeypatch)
    assert code == 2
    assert "foreign.py" in err


def test_session_file_is_allowed(repo, tmp_path, capsys, monkeypatch):
    (repo / "mine.py").write_text("x\n")
    t = _transcript(tmp_path, ["mine.py"])
    code, err = _run("git add mine.py", t, capsys, monkeypatch)
    assert code == 0, err


def test_fails_open_with_empty_session_set(repo, tmp_path, capsys, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    t = _transcript(tmp_path, [], name="empty.jsonl")
    code, _ = _run("git add foreign.py", t, capsys, monkeypatch)
    assert code == 0


def test_fails_open_with_missing_transcript(repo, capsys, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    code, _ = _run("git add foreign.py", "/nonexistent.jsonl", capsys, monkeypatch)
    assert code == 0


def test_description_bypass(repo, tmp_path, capsys, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    t = _transcript(tmp_path, ["mine.py"])
    code, _ = _run(
        "git add foreign.py", t, capsys, monkeypatch,
        description="intentional SKIP_SESSION_SCOPE_GUARD",
    )
    assert code == 0


def test_env_bypass(repo, tmp_path, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    t = _transcript(tmp_path, ["mine.py"])
    payload = {
        "tool_name": "Bash",
        "tool_input": {"command": "git add foreign.py"},
        "transcript_path": t,
    }
    text = json.dumps(payload)
    monkeypatch.setenv("SESSION_SCOPE_BYPASS", "1")
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
    assert HOOK.main() == 0


def test_non_git_commands_ignored(repo, tmp_path, capsys, monkeypatch):
    t = _transcript(tmp_path, ["mine.py"])
    code, err = _run("ls -la", t, capsys, monkeypatch)
    assert code == 0
    assert err == ""


def test_malformed_stdin_does_not_raise(monkeypatch):
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: "{oops")}))
    assert HOOK.main() == 0


def test_git_commit_blocks_foreign_staged_path(repo, tmp_path, capsys, monkeypatch):
    (repo / "foreign.py").write_text("x\n")
    _git(repo, "add", "foreign.py")
    t = _transcript(tmp_path, ["mine.py"])
    code, err = _run("git commit -m x", t, capsys, monkeypatch)
    assert code == 2
    assert "foreign.py" in err
    assert "git reset HEAD" in err
