"""Behavioral tests for the paper-context-guard PreToolUse hook.

The hook nudges the main thread toward reading full-text papers in a
subagent instead of inline. It is advisory: it must never block, and it
must never fire on ordinary file reads, or it becomes noise the user
learns to ignore.

Two properties are easy to get wrong and are pinned here:
  - host matching is suffix-anchored, so `biorxiv.org.attacker.com` is
    not treated as biorxiv.org
  - subagent invocations are exempt, because a subagent reading the
    paper is exactly the pattern the hook recommends
"""
import json

from conftest import REPO_ROOT, load_module_by_path

HOOK = load_module_by_path(
    "paper_context_guard",
    REPO_ROOT / "claude-config" / "hooks" / "paper-context-guard.py",
)


def _run(payload, capsys, monkeypatch, stdin_text=None):
    """Invoke the hook's main() with a payload; return (exit_code, stderr)."""
    monkeypatch.delenv("CLAUDE_AGENT_TYPE", raising=False)
    text = stdin_text if stdin_text is not None else json.dumps(payload)
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
    code = HOOK.main()
    return code, capsys.readouterr().err


def _warns(payload, capsys, monkeypatch) -> bool:
    code, err = _run(payload, capsys, monkeypatch)
    assert code == 0, "hook must never block"
    return "paper-context-guard" in err


def test_warns_on_local_pdf(capsys, monkeypatch):
    assert _warns(
        {"tool_name": "Read", "tool_input": {"file_path": "/papers/smith-2024.pdf"}},
        capsys, monkeypatch,
    )


def test_warns_on_publisher_url(capsys, monkeypatch):
    for url in (
        "https://www.biorxiv.org/content/10.1101/2024.01.01.573000v1.full",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC1234567/",
        "https://www.nature.com/articles/s41586-024-00000-0",
        "https://doi.org/10.1038/s41586-024-00000-0",
    ):
        assert _warns(
            {"tool_name": "WebFetch", "tool_input": {"url": url}}, capsys, monkeypatch
        ), url


def test_silent_on_ordinary_reads(capsys, monkeypatch):
    for path in ("/repo/README.md", "/repo/notes.txt", "/repo/data.csv"):
        assert not _warns(
            {"tool_name": "Read", "tool_input": {"file_path": path}}, capsys, monkeypatch
        ), path


def test_silent_on_non_paper_urls(capsys, monkeypatch):
    for url in ("https://github.com/foo/bar", "https://news.ycombinator.com/item?id=1"):
        assert not _warns(
            {"tool_name": "WebFetch", "tool_input": {"url": url}}, capsys, monkeypatch
        ), url


def test_host_match_is_suffix_anchored(capsys, monkeypatch):
    """A lookalike host must not inherit a real publisher's match."""
    for url in (
        "https://biorxiv.org.attacker.com/paper",
        "https://notnature.com/articles/x",
        "https://example.com/?ref=nature.com",
    ):
        assert not _warns(
            {"tool_name": "WebFetch", "tool_input": {"url": url}}, capsys, monkeypatch
        ), url


def test_subagent_is_exempt(capsys, monkeypatch):
    payload = {
        "tool_name": "Read",
        "tool_input": {"file_path": "/papers/smith-2024.pdf"},
        "is_subagent": True,
    }
    assert not _warns(payload, capsys, monkeypatch)


def test_subagent_env_signal_is_exempt(capsys, monkeypatch):
    monkeypatch.setenv("CLAUDE_AGENT_TYPE", "general-purpose")
    text = json.dumps(
        {"tool_name": "Read", "tool_input": {"file_path": "/papers/x.pdf"}}
    )
    monkeypatch.setattr("sys.stdin", type("S", (), {"read": staticmethod(lambda: text)}))
    assert HOOK.main() == 0
    assert "paper-context-guard" not in capsys.readouterr().err


def test_noqa_suppresses(capsys, monkeypatch):
    payload = {
        "tool_name": "Read",
        "tool_input": {
            "file_path": "/papers/x.pdf",
            "description": "read it inline noqa-paper-context",
        },
    }
    assert not _warns(payload, capsys, monkeypatch)


def test_other_tools_ignored(capsys, monkeypatch):
    assert not _warns(
        {"tool_name": "Bash", "tool_input": {"command": "cat paper.pdf"}},
        capsys, monkeypatch,
    )


def test_malformed_stdin_does_not_raise(capsys, monkeypatch):
    code, err = _run(None, capsys, monkeypatch, stdin_text="not json{")
    assert code == 0
    assert "paper-context-guard" not in err
