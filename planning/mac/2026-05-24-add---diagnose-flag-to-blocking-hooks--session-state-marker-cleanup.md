# Add --diagnose flag to blocking hooks + session-state marker cleanup

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

Two follow-ups carried over from the `_hook_lib.py` refactor entry:

- **B.** Add a `--diagnose` flag to every blocking hook so the user can probe "what would this hook do given this command?" without crafting JSON by hand AND without the hook actually interfering with a real tool call. Diagnose mode prints what would happen to stdout, prefixed with `[DIAGNOSE] would block:`, and always exits 0.

- **C.** Periodic cleanup of `session-start-context`'s marker files. The hook writes a per-session marker under `$CLAUDE_PROJECT_DIR/.claude/session-state/<session_id>-start.marker` so it fires only once per session. These accumulate forever without intervention. Reap markers older than 30 days at the start of every session — self-maintaining, no external cron.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add `is_diagnose_mode()` and `report(message, *, block)` helpers to `_hook_lib.py`. `report` returns the right exit code in either mode and routes output to stdout (diagnose) or stderr (normal).
- [x] Refactor blocking hooks to use `report()`:
  - `secret-leak-guard.py`
  - `rm-rf-sanity-guard.py`
  - `skill-frontmatter-validator.py`
  - `session-scope-guard.py` (one inline `block()` helper updated; 3 call sites untouched)
  - `claude-config-planning-required.py`
- [x] Add `_cleanup_old_markers(marker_dir)` to `session-start-context.py`; call it at the top of `main()` before the existence check
- [x] Smoke-test diagnose mode for all 4 blocking PreToolUse hooks + normal-mode regression
- [x] Smoke-test marker cleanup with a 7-month-old file + a fresh file
- [x] `./sync-config.py push` and verify live

## Implementation

### `--diagnose` flag (via `_hook_lib.report`)

`_hook_lib.py` gains two ~10-line helpers:

```python
def is_diagnose_mode() -> bool:
    return "--diagnose" in sys.argv[1:]

def report(message: str, *, block: bool = True) -> int:
    if is_diagnose_mode():
        prefix = "[DIAGNOSE] would block:\n" if block else "[DIAGNOSE] would warn:\n"
        print(prefix + message)
        return 0
    print(message, file=sys.stderr)
    return 2 if block else 0
```

Each blocking hook's `block(...)` / inline `print + return 2` was replaced with `return report(message, block=True)`. The user can now invoke any hook with `--diagnose` appended and pipe in a synthetic payload to see what it would do:

```bash
echo '{"tool_name":"Bash","tool_input":{"command":"rm -rf ~"}}' | \
  ~/.claude/hooks/rm-rf-sanity-guard.py --diagnose
# [DIAGNOSE] would block:
# [rm-rf-sanity-guard] Refusing destructive command — dangerous target(s):
#   - path: '~'  reason: $HOME exactly
# exit 0
```

`session-scope-guard.py` already had a `block()` helper; it was rewritten to delegate to `report()` so all 3 call sites picked up diagnose mode for free. The other 4 hooks had inline `print + return 2` patterns updated case-by-case.

### Marker cleanup

`session-start-context.py` gains:

```python
MARKER_TTL_DAYS = 30

def _cleanup_old_markers(marker_dir: Path) -> None:
    cutoff = time.time() - MARKER_TTL_DAYS * 86400
    if not marker_dir.exists():
        return
    for entry in marker_dir.iterdir():
        if not entry.name.endswith("-start.marker"):
            continue
        try:
            if entry.stat().st_mtime < cutoff:
                entry.unlink()
        except OSError:
            continue
```

Called at the top of `main()` before the existence-check on the current session's marker. Cost: one `iterdir` + one `stat` per existing marker. With <100 markers (typical for active use), well under 10ms.

## Expected Outcome

- Diagnose mode: any blocking hook can be probed without affecting real tool calls. Useful for debugging false positives and writing future test cases.
- Marker dir stays bounded automatically — no external sweeper needed.
- All existing behavior unchanged in normal mode (the `report()` helper produces identical stderr output and exit codes when `--diagnose` is absent).

## Actual Outcome

All 11 files (`_hook_lib.py` + 10 hooks) compile cleanly. Diagnose-mode smoke tests:

| Hook | Test | Output |
|---|---|---|
| `rm-rf-sanity-guard --diagnose` | `rm -rf ~` | `[DIAGNOSE] would block:` + reason; exit 0 |
| `secret-leak-guard --diagnose` | fake Anthropic key | `[DIAGNOSE] would block:` + match list; exit 0 |
| `session-scope-guard --diagnose` | `git add -A` | `[DIAGNOSE] would block:` + bulk-stage reason; exit 0 |
| `skill-frontmatter-validator --diagnose` | SKILL.md missing description | `[DIAGNOSE] would block:` + reason; exit 0 |
| (normal sanity) `rm-rf-sanity-guard` | `rm -rf ~` | normal block, exit 2 |

Marker cleanup test:

```
before: fresh-session-start.marker (today), old-session-start.marker (Oct 2024)
after:  fresh-session-start.marker (today), test-new-start.marker (today)
```

Old marker reaped, fresh marker kept, new session's marker created. Correct.

## Assessment

**Result**: Success

**Improvements**:
- The `report()` helper centralizes the block/warn/diagnose logic so future blocking hooks get diagnose mode for free.
- The user (or me, during future debugging) can pipe a synthetic payload into any hook with `--diagnose` to see what it would do, without modifying the hook code or temporarily disabling it.
- `session-state` dir is now self-bounding; no manual sweep needed.

**Issues**:
- None observed. The `report()` signature includes a `block=False` mode for informational hooks; none of the existing informational hooks (`scratch-redirect-suggest`, `ipynb-output-warn`, `stale-todo-tracker`, `claude-md-improver-suggest`) use it yet — they still print directly to stderr. Refactoring those is optional and low-value; deferred.

**Lessons Learned**:
- `CLAUDE_PROJECT_DIR=$X echo ... | python` only sets `CLAUDE_PROJECT_DIR` for `echo`, not for the python downstream of the pipe. Use `echo ... | CLAUDE_PROJECT_DIR=$X python` instead — env vars before the command apply to that command only.
- A "diagnose mode" via `--diagnose` flag is the cheapest way to make a blocking hook debuggable. No need for `--dry-run` style schemes that duplicate the decision logic; just route the same output differently.

## Related Commits

- `6712064`: feat(hooks): add --diagnose flag (via _hook_lib.report) + session-state marker cleanup

## Next Steps

- Optionally extend the `report()` pattern to the informational hooks (`scratch-redirect-suggest`, `ipynb-output-warn`, etc.) so they too support `--diagnose`. Low priority — informational hooks don't have false-positive debugging needs since they never block.
- If hooks accumulate even more verbose decision-making, consider a `--verbose` flag that prints a trace of every decision step. For now `--diagnose` (which prints the final verdict + reason) covers the common case.
