"""Behavioral tests for the context-budget-warn Stop hook.

The hook reports session context size once per band. It is advisory: it must
never block, and a band must not re-fire, or every turn past 150K becomes noise.

Properties pinned here:
  - bands fire once, in ascending order, and never regress
  - the size is read from the tail of the transcript, so a huge transcript with
    the usage record far from the start still works
  - a missing, empty, or usage-free transcript is silent rather than a crash
"""
import json

from conftest import REPO_ROOT, load_module_by_path

HOOK = load_module_by_path(
    "context_budget_warn",
    REPO_ROOT / "claude-config" / "hooks" / "context-budget-warn.py",
)


def _transcript(tmp_path, *sizes, pad_lines=0, name="t.jsonl"):
    """Write a transcript whose usage records carry the given cache-read sizes."""
    p = tmp_path / name
    lines = []
    for _ in range(pad_lines):
        lines.append(json.dumps({"message": {"content": [{"type": "text", "text": "x" * 400}]}}))
    for s in sizes:
        lines.append(json.dumps({
            "message": {"usage": {
                "cache_read_input_tokens": s,
                "cache_creation_input_tokens": 0,
                "input_tokens": 0,
            }}
        }))
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def _run(capsys, monkeypatch, tmp_path, transcript, session="s1"):
    payload = {"session_id": session, "transcript_path": transcript}
    text = json.dumps(payload)
    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
    code = HOOK.main()
    assert code == 0, "hook must never block"
    return capsys.readouterr().err


def test_silent_below_first_band(capsys, monkeypatch, tmp_path):
    t = _transcript(tmp_path, 90_000)
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, t)


def test_warns_at_first_band(capsys, monkeypatch, tmp_path):
    t = _transcript(tmp_path, 160_000)
    err = _run(capsys, monkeypatch, tmp_path, t)
    assert "context-budget-warn" in err
    assert "160K" in err


def test_band_does_not_refire(capsys, monkeypatch, tmp_path):
    t = _transcript(tmp_path, 160_000)
    assert "context-budget-warn" in _run(capsys, monkeypatch, tmp_path, t)
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, t)


def test_escalates_to_next_band(capsys, monkeypatch, tmp_path):
    t1 = _transcript(tmp_path, 160_000, name="a.jsonl")
    assert "context-budget-warn" in _run(capsys, monkeypatch, tmp_path, t1)
    t2 = _transcript(tmp_path, 260_000, name="b.jsonl")
    err = _run(capsys, monkeypatch, tmp_path, t2)
    assert "context-budget-warn" in err
    assert "260K" in err


def test_does_not_warn_again_when_context_shrinks(capsys, monkeypatch, tmp_path):
    t1 = _transcript(tmp_path, 260_000, name="a.jsonl")
    assert "context-budget-warn" in _run(capsys, monkeypatch, tmp_path, t1)
    t2 = _transcript(tmp_path, 160_000, name="b.jsonl")
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, t2)


def test_uses_last_usage_record(capsys, monkeypatch, tmp_path):
    t = _transcript(tmp_path, 400_000, 100_000)
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, t)


def test_sums_usage_components(capsys, monkeypatch, tmp_path):
    p = tmp_path / "sum.jsonl"
    p.write_text(json.dumps({
        "message": {"usage": {
            "cache_read_input_tokens": 100_000,
            "cache_creation_input_tokens": 55_000,
            "input_tokens": 2_000,
        }}
    }) + "\n")
    err = _run(capsys, monkeypatch, tmp_path, str(p))
    assert "context-budget-warn" in err
    assert "157K" in err


def test_finds_usage_beyond_tail_window(capsys, monkeypatch, tmp_path):
    """A large transcript still works: the tail window must cover the last record."""
    t = _transcript(tmp_path, 300_000, pad_lines=2000)
    err = _run(capsys, monkeypatch, tmp_path, t)
    assert "context-budget-warn" in err
    assert "300K" in err


def test_missing_transcript_is_silent(capsys, monkeypatch, tmp_path):
    assert "context-budget-warn" not in _run(
        capsys, monkeypatch, tmp_path, str(tmp_path / "nope.jsonl")
    )


def test_empty_transcript_is_silent(capsys, monkeypatch, tmp_path):
    p = tmp_path / "empty.jsonl"
    p.write_text("")
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, str(p))


def test_transcript_without_usage_is_silent(capsys, monkeypatch, tmp_path):
    p = tmp_path / "nousage.jsonl"
    p.write_text(json.dumps({"message": {"content": [{"type": "text", "text": "hi"}]}}) + "\n")
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, str(p))


def test_malformed_lines_are_skipped(capsys, monkeypatch, tmp_path):
    p = tmp_path / "mixed.jsonl"
    p.write_text(
        "{not json\n"
        + json.dumps({"message": {"usage": {"cache_read_input_tokens": 160_000}}}) + "\n"
        + "{cache_read_input_tokens: broken\n"
    )
    assert "context-budget-warn" in _run(capsys, monkeypatch, tmp_path, str(p))


def test_muted_session_is_silent(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path))
    HOOK._save(HOOK._state_path("s1"), {"muted": True})
    t = _transcript(tmp_path, 400_000)
    assert "context-budget-warn" not in _run(capsys, monkeypatch, tmp_path, t)


def test_malformed_stdin_does_not_raise(capsys, monkeypatch, tmp_path):
    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path))
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: "{oops")}))
    assert HOOK.main() == 0
    assert "context-budget-warn" not in capsys.readouterr().err
