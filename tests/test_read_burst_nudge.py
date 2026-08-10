"""Behavioral tests for the read-burst-nudge PostToolUse hook.

The hook nudges the main thread toward delegating broad file sweeps. It is
advisory: it must never block, and it must fire at most once per session or it
becomes noise the user learns to ignore.

Properties pinned here:
  - it stays silent below the threshold and fires exactly once at it
  - a session that has already dispatched a subagent is never nudged
  - subagents are exempt, because a subagent reading many files is the
    recommended pattern
  - the common path does no transcript I/O (the counter is O(1)); the
    transcript is walked at most once per session
"""
import json

from conftest import REPO_ROOT, load_module_by_path

HOOK = load_module_by_path(
    "read_burst_nudge",
    REPO_ROOT / "claude-config" / "hooks" / "read-burst-nudge.py",
)


def _run(payload, capsys, monkeypatch, tmp_path, stdin_text=None):
    """Invoke main() with a payload; return (exit_code, stderr)."""
    monkeypatch.delenv("CLAUDE_AGENT_TYPE", raising=False)
    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path))
    text = stdin_text if stdin_text is not None else json.dumps(payload)
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
    code = HOOK.main()
    return code, capsys.readouterr().err


def _read_payload(session="s1", transcript=""):
    return {
        "tool_name": "Read",
        "tool_input": {"file_path": "/x/y.py"},
        "session_id": session,
        "transcript_path": transcript,
    }


def _pump(n, capsys, monkeypatch, tmp_path, session="s1", transcript=""):
    """Send n read events; return the list of stderr strings."""
    out = []
    for _ in range(n):
        code, err = _run(_read_payload(session, transcript), capsys, monkeypatch, tmp_path)
        assert code == 0, "hook must never block"
        out.append(err)
    return out


def _transcript_with_dispatch(tmp_path):
    p = tmp_path / "t.jsonl"
    p.write_text(json.dumps({
        "message": {"content": [{"type": "tool_use", "name": "Agent", "input": {}}]}
    }) + "\n")
    return str(p)


def test_silent_below_threshold(capsys, monkeypatch, tmp_path):
    errs = _pump(HOOK.THRESHOLD - 1, capsys, monkeypatch, tmp_path)
    assert not any("read-burst-nudge" in e for e in errs)


def test_fires_at_threshold(capsys, monkeypatch, tmp_path):
    errs = _pump(HOOK.THRESHOLD, capsys, monkeypatch, tmp_path)
    assert sum("read-burst-nudge" in e for e in errs) == 1


def test_fires_only_once_per_session(capsys, monkeypatch, tmp_path):
    errs = _pump(HOOK.THRESHOLD + 8, capsys, monkeypatch, tmp_path)
    assert sum("read-burst-nudge" in e for e in errs) == 1


def test_no_nudge_when_already_dispatched(capsys, monkeypatch, tmp_path):
    t = _transcript_with_dispatch(tmp_path)
    errs = _pump(HOOK.THRESHOLD + 3, capsys, monkeypatch, tmp_path, transcript=t)
    assert not any("read-burst-nudge" in e for e in errs)


def test_subagent_is_exempt(capsys, monkeypatch, tmp_path):
    payload = _read_payload()
    payload["agent_type"] = "general-purpose"
    for _ in range(HOOK.THRESHOLD + 3):
        code, err = _run(payload, capsys, monkeypatch, tmp_path)
        assert code == 0
        assert "read-burst-nudge" not in err


def test_subagent_env_signal_is_exempt(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("CLAUDE_AGENT_TYPE", "Explore")
    for _ in range(HOOK.THRESHOLD + 3):
        text = json.dumps(_read_payload())
        monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path))
        monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
        assert HOOK.main() == 0
        assert "read-burst-nudge" not in capsys.readouterr().err


def test_noqa_suppresses(capsys, monkeypatch, tmp_path):
    payload = _read_payload()
    payload["tool_input"]["description"] = "sweep noqa-read-burst"
    for _ in range(HOOK.THRESHOLD + 3):
        code, err = _run(payload, capsys, monkeypatch, tmp_path)
        assert code == 0
        assert "read-burst-nudge" not in err


def test_other_tools_ignored(capsys, monkeypatch, tmp_path):
    payload = _read_payload()
    payload["tool_name"] = "Bash"
    for _ in range(HOOK.THRESHOLD + 3):
        code, err = _run(payload, capsys, monkeypatch, tmp_path)
        assert code == 0
        assert "read-burst-nudge" not in err


def test_grep_and_glob_count_toward_threshold(capsys, monkeypatch, tmp_path):
    fired = 0
    for i in range(HOOK.THRESHOLD):
        payload = _read_payload()
        payload["tool_name"] = ("Grep", "Glob", "Read")[i % 3]
        code, err = _run(payload, capsys, monkeypatch, tmp_path)
        assert code == 0
        fired += "read-burst-nudge" in err
    assert fired == 1


def test_sessions_are_tracked_independently(capsys, monkeypatch, tmp_path):
    _pump(HOOK.THRESHOLD - 1, capsys, monkeypatch, tmp_path, session="alpha")
    errs = _pump(1, capsys, monkeypatch, tmp_path, session="beta")
    assert not any("read-burst-nudge" in e for e in errs)


def test_malformed_stdin_does_not_raise(capsys, monkeypatch, tmp_path):
    code, err = _run(None, capsys, monkeypatch, tmp_path, stdin_text="{not json")
    assert code == 0
    assert "read-burst-nudge" not in err


def test_common_path_does_no_transcript_io(capsys, monkeypatch, tmp_path):
    """Below the threshold the hook must not walk the transcript.

    This is the perf contract: PostToolUse on Read fires thousands of times per
    session, so the hot path has to stay O(1).
    """
    calls = []
    monkeypatch.setattr(HOOK, "_has_dispatched", lambda p: calls.append(p) or False)
    _pump(HOOK.THRESHOLD - 1, capsys, monkeypatch, tmp_path)
    assert calls == []
    _pump(1, capsys, monkeypatch, tmp_path)
    assert len(calls) == 1, "transcript walked exactly once, at the threshold"
