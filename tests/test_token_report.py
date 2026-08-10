"""Tests for tools/token-report.py.

The tool's job is to be trustworthy about a number decisions get made on
(cache read per request), so the aggregation and the transcript parsing are
pinned here. Malformed transcripts must degrade rather than crash — they are
live session files and can be truncated mid-write.
"""
import json

from conftest import REPO_ROOT, load_module_by_path

TOOL = load_module_by_path("token_report", REPO_ROOT / "tools" / "token-report.py")


def _write(tmp_path, name, events):
    d = tmp_path / "proj"
    d.mkdir(exist_ok=True)
    p = d / name
    p.write_text("".join(json.dumps(e) + "\n" for e in events))
    return p


def _usage(cache_read=0, output=0, creation=0, inp=0):
    return {"message": {"usage": {
        "cache_read_input_tokens": cache_read,
        "cache_creation_input_tokens": creation,
        "output_tokens": output,
        "input_tokens": inp,
    }}}


def _tool_use(name, **inputs):
    return {"message": {"content": [{"type": "tool_use", "name": name, "input": inputs}]}}


def test_counts_requests_and_totals(tmp_path):
    p = _write(tmp_path, "a.jsonl", [_usage(1000, 10), _usage(2000, 20)])
    s = TOOL.scan(p)
    assert s.requests == 2
    assert s.cache_read == 3000
    assert s.totals["output_tokens"] == 30
    assert s.per_request == 1500


def test_per_request_is_zero_without_requests(tmp_path):
    p = _write(tmp_path, "b.jsonl", [_tool_use("Read", file_path="/x")])
    assert TOOL.scan(p).per_request == 0


def test_records_skill_and_dispatch_names(tmp_path):
    p = _write(tmp_path, "c.jsonl", [
        _usage(100),
        _tool_use("Skill", skill="papers-library"),
        _tool_use("Skill", skill="papers-library"),
        _tool_use("Agent", subagent_type="Explore"),
        _tool_use("Task", subagent_type="general-purpose"),
        _tool_use("Read", file_path="/x"),
    ])
    s = TOOL.scan(p)
    assert s.skills == {"papers-library": 2}
    assert s.dispatches == {"Explore": 1, "general-purpose": 1}


def test_malformed_lines_are_skipped(tmp_path):
    d = tmp_path / "proj"
    d.mkdir(exist_ok=True)
    p = d / "d.jsonl"
    p.write_text("{broken\n" + json.dumps(_usage(500)) + "\n\n{also broken\n")
    s = TOOL.scan(p)
    assert s.requests == 1
    assert s.cache_read == 500


def test_missing_file_yields_no_requests(tmp_path):
    assert TOOL.scan(tmp_path / "absent.jsonl").requests == 0


def test_summary_aggregates_across_sessions(tmp_path):
    _write(tmp_path, "e.jsonl", [_usage(1000, 100), _tool_use("Skill", skill="editor")])
    _write(tmp_path, "f.jsonl", [_usage(3000, 100)])
    sessions = TOOL.collect(tmp_path / "proj", None)
    summary = TOOL.build_summary(sessions)
    assert summary["sessions"] == 2
    assert summary["requests"] == 2
    assert summary["cache_read_per_request"] == 2000
    assert summary["output_per_request"] == 100
    assert summary["pct_sessions_with_skill"] == 50.0
    assert summary["skill_calls"] == 1


def test_collect_skips_transcripts_without_usage(tmp_path):
    _write(tmp_path, "g.jsonl", [_tool_use("Read", file_path="/x")])
    assert TOOL.collect(tmp_path / "proj", None) == []


def test_delta_direction():
    assert "better" in TOOL._delta(50, 100)
    assert "worse" in TOOL._delta(150, 100)
    assert "better" in TOOL._delta(150, 100, lower_is_better=False)
    assert "flat" in TOOL._delta(100, 100)
    assert TOOL._delta(10, 0) == ""


def test_main_reports_and_exits_zero(tmp_path, capsys):
    _write(tmp_path, "h.jsonl", [_usage(1000, 100)])
    assert TOOL.main(["--projects-dir", str(tmp_path / "proj")]) == 0
    out = capsys.readouterr().out
    assert "per request" in out
    assert "cache read" in out


def test_main_json_mode_includes_baseline(tmp_path, capsys):
    _write(tmp_path, "i.jsonl", [_usage(1000, 100)])
    assert TOOL.main(["--projects-dir", str(tmp_path / "proj"), "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["baseline"]["cache_read_per_request"] == 163573
    assert payload["current"]["requests"] == 1


def test_main_rejects_bad_since(tmp_path, capsys):
    _write(tmp_path, "j.jsonl", [_usage(1000)])
    code = TOOL.main(["--projects-dir", str(tmp_path / "proj"), "--since", "11-08-2026"])
    assert code == 1
    assert "YYYY-MM-DD" in capsys.readouterr().err


def test_main_errors_on_missing_dir(tmp_path, capsys):
    assert TOOL.main(["--projects-dir", str(tmp_path / "nope")]) == 1
    assert "no transcript directory" in capsys.readouterr().err
