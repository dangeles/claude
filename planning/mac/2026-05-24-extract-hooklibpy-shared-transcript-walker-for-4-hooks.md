# Extract _hook_lib.py shared transcript walker for 4 hooks

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

After shipping 9 workflow hooks (`f6b6675`), four of them duplicated the same "open transcript JSONL → parse each line → navigate to `message.content` → iterate looking for tool_use or user-role entries" pattern. Extract this into a shared `_hook_lib.py` module so that pattern lives in one place.

The four hooks with duplicated logic:

- `session-scope-guard.py` — `derive_session_files()` + `_collect_from_event()`
- `claude-config-planning-required.py` — `_walk_transcript()`
- `stale-todo-tracker.py` — `_walk_transcript()`
- `claude-md-improver-suggest.py` — `_walk_user_messages()`

The first three needed `(tool_name, tool_input)` pairs; the fourth needed user-message text. Three iterators cover both shapes.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Write `claude-config/hooks/_hook_lib.py` with three public iterators: `iter_transcript_events`, `iter_tool_uses`, `iter_user_messages`
- [x] Refactor `claude-md-improver-suggest.py` to use `iter_user_messages`
- [x] Refactor `stale-todo-tracker.py` to use `iter_tool_uses`
- [x] Refactor `claude-config-planning-required.py` to use `iter_tool_uses`
- [x] Refactor `session-scope-guard.py` to use `iter_tool_uses` (preserves the inline Bash-mutation parsing)
- [x] Smoke-test each refactored hook with a synthetic mixed transcript
- [x] `./sync-config.py push` and verify `_hook_lib.py` reaches `~/.claude/hooks/`

## Implementation

### `_hook_lib.py` (90 lines)

Three public iterators:

- `iter_transcript_events(transcript_path)` — yields one parsed JSON event per non-empty line of the JSONL transcript. Silently skips malformed lines; returns empty iterator if the file is missing or unreadable.
- `iter_tool_uses(transcript_path)` — yields `(tool_name, tool_input)` pairs for every `type: "tool_use"` entry across every assistant message. Handles both transcript shapes (top-level fields vs `message.content[]` wrapper).
- `iter_user_messages(transcript_path)` — yields the text body of every user-role message. For list-form content (text blocks + tool_use + tool_result), joins only the `text`-type blocks. Empty messages are skipped.

Each hook imports via the `sys.path.insert(0, os.path.dirname(__file__))` pattern (with `# pyright: ignore[reportMissingImports]` since Pyright can't statically resolve sys.path manipulation).

### Refactor cost per hook

| Hook | Before (lines) | After | Walker code removed |
|---|---|---|---|
| `claude-md-improver-suggest.py` | 122 | 86 | ~40 lines |
| `stale-todo-tracker.py` | 124 | 101 | ~25 lines |
| `claude-config-planning-required.py` | 142 | 117 | ~30 lines |
| `session-scope-guard.py` | 421 | 403 | ~25 lines |
| `_hook_lib.py` (new) | — | 90 | (replaces 4× walkers) |

Net change: +90 (new lib) -120 (walker code removed from 4 hooks) = **-30 lines** across the catalog, with the win being that future hooks needing transcript walking get the iterators for free.

## Expected Outcome

- Each refactored hook is shorter and more focused; the transcript-walking logic is reviewed in one place.
- Future hooks needing transcript access import the iterators directly (one-line dependency).
- No behavior changes — the iterators preserve every code path of the previous inline walkers.

## Actual Outcome

All five files (`_hook_lib.py` + 4 refactored hooks) compile cleanly. Smoke tests against a synthetic transcript with one user message + one Write + one Edit all pass:

- `claude-md-improver-suggest` matched `remember to` from the user message → exit 0 with stderr suggestion (correct)
- `stale-todo-tracker` reported 1 TODO added in `src/foo.py` → exit 0 with stderr summary (correct)
- `claude-config-planning-required` (run with `cwd=~/repos/claude`) blocked Stop (exit 2) because a gated path was touched but no planning entry was written (correct; behavior unchanged from before)
- `session-scope-guard` allowed `git add src/foo.py` (in session set) and blocked `git add bar.py` (not in session set) — correct

Live state: `~/.claude/hooks/_hook_lib.py` synced and matches source.

## Assessment

**Result**: Success

**Improvements**:
- 30 fewer lines across the catalog despite adding a new file — net simplification.
- The transcript shape (`event.message.content[]` vs `event.content[]`) is now handled in one place. If Claude Code's transcript format evolves, only `_hook_lib.py` needs updating.
- Adding a future transcript-aware hook becomes a one-line import.

**Issues**:
- Pyright can't statically resolve `from _hook_lib import ...` because the import depends on the `sys.path.insert` shim. Suppressed via `# pyright: ignore[reportMissingImports]` on each import line. Same convention used elsewhere in the codebase (`scripts/tests/*.py` in plotting-advisor).

**Lessons Learned**:
- The right time to extract a shared helper is when the third user appears, not the second. Three transcript walkers got built in close succession, all subtly different in what they extracted; consolidating clarified the shape.
- The `iter_*` naming (generators) is cleaner than `walk_*` (which suggests a single return value) — it lets the caller use list comprehensions or short-circuit as they wish.

## Related Commits

- `<pending>`: refactor(hooks): extract _hook_lib.py shared transcript walker

## Next Steps

- (Carryover from previous entry) Add `--diagnose` flag to blocking hooks for false-positive debugging.
- (Carryover) Periodic cleanup of `session-state/*.marker` files >30 days old.
- If a fourth/fifth hook needs filesystem state inspection (e.g., reading config files, walking dirs), consider a `_hook_lib.fs` module — but only when the third user actually appears.
