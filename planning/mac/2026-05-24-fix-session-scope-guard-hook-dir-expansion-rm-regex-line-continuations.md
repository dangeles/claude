# Fix session-scope-guard hook: dir expansion, rm regex, line continuations

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

Three gaps in the new `session-scope-guard.py` PreToolUse hook surfaced during dogfooding immediately after it shipped (`dee3a86`). All three caused legitimate session commits to fail the guard:

1. **Directory args don't expand** — `git add some/dir/` checked the literal `some/dir/` string against the session set; it should expand to the actual files git would stage under that directory and check each.
2. **`rm` regex only captures the first path arg** — `rm -rf a b c` recorded only `a` in the session set, so file-level deletes inside `b/` and `c/` look out-of-session at commit time.
3. **Backslash line-continuations leak as path tokens** — `git add foo \\\n bar` had `\\` parsed as a stray path that the hook then refused.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Strip backslash-newline continuations before any path parsing
- [x] Replace per-regex single-capture for `rm`/`mv`/`cp`/`tee` with a "capture args segment + shlex split" approach so multi-arg commands record every path
- [x] Add `expand_path_for_session_check()` that uses `git status --porcelain -- <path>` to enumerate the files git would actually stage when given a directory; the per-file check then runs against that expanded list
- [x] Smoke-test all three fixes with synthetic stdin payloads + a tmp git repo
- [x] `./sync-config.py push` and verify live hook matches source

## Implementation

### Fix #1 — backslash line-continuations

Introduced `_strip_line_continuations(command)` which collapses `\\\s*\n` to a single space. Called from both `extract_git_add_paths` and `_extract_bash_mutation_paths` before any tokenization. Also added `re.DOTALL` to the `git add` capture regex so multi-line commands are captured fully (not stopped at the first newline).

Added a defensive filter — even if a stray `\` token survives, it's dropped from the path list explicitly.

### Fix #2 — `rm`/`mv`/`cp`/`tee` multi-arg capture

Replaced per-pattern single-group capture with a `(regex, mode)` table where `mode` controls how the captured args segment is interpreted:

| mode | meaning | example |
|---|---|---|
| `all` | every non-flag token is a mutated path | `rm -rf a b c` → `[a, b, c]` |
| `last` | only the last token is a mutated path (dest of mv/cp) | `mv src1 src2 dest/` → `[dest/]` |
| `non-flag` | every non-flag token (tee writes to all listed paths) | `tee -a out1 out2` → `[out1, out2]` |
| `first` | only the first non-flag token (sed -i, > redirect) | `sed -i 's/x/y/' file.txt` → `[file.txt]` |

The args segment is captured via `(.+?)(?:&&\|;\|\|\|\|\|$)` (non-greedy, terminating at command-chaining boundaries) and split with shlex.

### Fix #3 — directory expansion in `git add`

New `expand_path_for_session_check(path, repo_root)` function:

- If `path` ends with `/` OR `os.path.isdir(path)` resolves locally, treat it as a directory argument.
- Otherwise, return `[path]` unchanged.
- For directory args, run `git status --porcelain -- <dir>` to enumerate dirty files git would actually stage. Parse the porcelain output (XY status + space + path, with rename support via `X -> Y`) and return the full list.
- If the porcelain call fails (subprocess error, not a git repo), fall back to `[path]` so the upstream check still runs on the literal directory path (matches the pre-fix behavior).

Called from the per-path `git add` check loop: each raw path is expanded, then each expanded file is canonicalized and checked against the session set.

### Bonus fix — porcelain leading-space stripping

While testing #3, found a 1-char path truncation bug: `out.strip()` on the porcelain output was eating the leading space that `git status --porcelain` emits for unstaged-modified files (` M somedir/file1.txt`). Path then started at index 3 of `M somedir/file1.txt` instead of ` M somedir/file1.txt` — producing `omedir/file1.txt`. Fixed by using `out.rstrip("\n")` instead of `out.strip()`.

## Expected Outcome

- `git add somedir/` (directory) checks each file git would actually stage instead of the literal directory string.
- `rm -rf a b c` records all three paths in the session set so commits that include file-level deletes inside any of those dirs are recognized as in-session.
- `git add foo \\\n bar` parses cleanly, no spurious `\\` token.
- Existing protections (bulk-stage refusal, commit-time out-of-session check, fail-OPEN on uncertainty) unchanged.

## Actual Outcome

Three direct smoke tests + one regression test for the bonus fix, all pass:

| # | Scenario | Expected | Actual |
|---|---|---|---|
| 1 | `git add a.md \\\n b.md` with both files in session | exit 0 | exit 0 |
| 2 | `_extract_bash_mutation_paths("rm -rf a b c")` | `[a, b, c]` | `[a, b, c]` |
| 3a | `git add somedir/` with file1 in session, file2 untracked | exit 2 listing only `somedir/file2.txt` | matches |
| 3b | `git add somedir/` with both files in session | exit 0 | exit 0 |
| Bonus | porcelain output's leading space preserved | `somedir/file1.txt` (full) | `somedir/file1.txt` |

Live state: `~/.claude/hooks/session-scope-guard.py` matches source (`diff -q` confirms).

## Assessment

**Result**: Success

**Improvements**:
- The hook is now safe to leave on for normal multi-file commits — directory args, multi-line `git add`, and `rm -rf a b c` patterns all work correctly without requiring `SKIP_SESSION_SCOPE_GUARD`.
- The `(regex, mode)` table for Bash mutations is more extensible than the per-pattern approach. Adding a new mutation tool (e.g., `truncate -s 0 file`, `chmod +x file`) is one row.
- The directory-expansion path uses git itself as the source of truth ("what files would `git add <dir>` actually stage?"), which is more robust than any in-process heuristic.

**Issues**:
- None observed. The porcelain-output parsing depends on the format being `XY path` (3-char prefix), which is stable across git versions but worth re-checking if git's porcelain format changes.

**Lessons Learned**:
- `.strip()` on multi-line output is a footgun when the format relies on leading whitespace per-line. Prefer `.rstrip("\n")` or process line-by-line.
- Dogfooding a new hook against a real commit is the fastest way to catch gaps — the three issues showed up in <60 seconds of trying to commit the hook itself, all in one error message.
- Replacing per-regex capture with "regex + post-processing mode" makes the table-driven design more honest about what each Bash command actually does (rm mutates all args; mv mutates the dest; sed -i mutates the file).

## Related Commits

- `740204a`: fix(session-scope-guard): expand dir args, capture multi-arg rm/mv/cp, strip line continuations

## Next Steps

- If transcript format changes (Claude Code version bump), update the `_collect_from_event` walker.
- The `expand_path_for_session_check` could also be applied to `git rm <path>` (currently not covered) — but `git rm` is rare and the bulk-stage refusal already catches the `git rm -r somedir/` style.
- Consider adding a `--diagnose` mode on the hook that prints what it would do given an explicit command, useful for debugging without a real Bash invocation.
