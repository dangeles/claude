# Post-audit follow-ups (orphans, scope-guard subagents, version field, drift guard)

**Date**: 2026-06-09
**Machine**: mac
**Status**: Implemented

## Objective

Act on the five post-audit recommendations plus a superpowers-currency check.

## Changes

- [x] **#1 Orphan cleanup.** `push --delete` removed 4 confirmed-stale live-only files
  (git-strategy-advisor examples/refs dropped in dee3a86; workflow-coordinator
  `universal-handoff-schema-v3.0.json` rename leftover). Also renamed 3 live files to
  lowercase to match the repo (`SOURCES.md`→`sources.md`, `END-TO-END-TEST-PLAN.md`,
  `TEST-RESULTS.md`) — these case-only artifacts had been showing as perpetual phantom
  "3 new, 3 deleted" drift on a case-insensitive filesystem.
- [x] **#2 session-scope-guard ↔ subagents.** `derive_session_files()` now walks
  subagent transcripts (`<session-dir>/**/subagents/*.jsonl`) in addition to the main
  transcript, so files edited by dispatched agents (Task tool) count as session-modified
  instead of tripping the guard. Verified on this session's real transcripts: subagent-only
  edits go False→True; main-loop edits stay True; 6 transcripts walked.
- [x] **#3 Triggering check.** Dispatched a router subagent over the reworded P4
  descriptions and the two deprecations (see Actual Outcome).
- [x] **#4 version-field policy.** Decided: keep `version` OPTIONAL and meaningful (bump on
  notable rewrites; CHANGELOG is the authority). Did NOT mass-drop it — that would destroy
  CHANGELOG-referenced version history (skill-editor v3.0, scientific-analysis-architect
  2.0.0, etc.). Tightened the README frontmatter-conventions row to remove the ambiguity.
- [x] **#5 Drift guard.** Added `_drift_warning()` to `session-start-context.py`: when cwd
  is the config repo (self-scoped by presence of `sync-config.py`), it appends a one-line
  warning to the session-context block if `~/.claude` diverges from `claude-config/`.
  Best-effort, 15s timeout, fail-silent. Verified: warns on real drift, None elsewhere.
- [x] **Superpowers currency.** Refreshed the `claude-plugins-official` marketplace and ran
  `claude plugin update superpowers` → already at latest (5.1.0). The installed-vs-manifest
  sha mismatch was a red herring (resolver pins the 5.1.0 tag). No update needed.

## Actual Outcome

Hooks + README synced to `~/.claude`; status clean afterward. Both hooks compile and were
unit-tested against real data. The scope-guard fix removes the friction that forced a
`SKIP_SESSION_SCOPE_GUARD` bypass on the audit commit.

## Assessment

**Result**: Success

**Lessons Learned**:
- A case-insensitive filesystem turns a lowercase-rename commit into permanent phantom
  drift until the live dirents are also re-cased; worth doing as part of any rename pass.
- Subagent-driven edits are invisible to transcript-scoped guards unless the guard walks
  subagent transcripts — a general gotcha for any hook that reasons over "what changed."

## Related Commits

- 968cb59: chore(hooks): scope-guard walks subagent transcripts + session-start drift guard

## Next Steps

- Monitor the drift guard for false positives / latency in daily use.
