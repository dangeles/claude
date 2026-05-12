# tier-1-cleanup-claude-md-and-hostname-normalize

**Date**: 2026-05-12
**Machine**: mac
**Status**: Success

## Objective

Tier 1 of a 3-tier improvement pass on the Claude config repo. Three independent low-risk fixes:

1. **A1**: `.claude/settings.local.json` had accumulated 328 entries with ~65% garbage (shell fragments parsed as permissions, `__NEW_LINE_*` parser artifacts, one-shot path tests, 99 WebFetch domains, etc.) — bloating permission lookups and visually obscuring the actual policy.
2. **A2**: No root `CLAUDE.md`, so Claude doesn't auto-load the 7-step config workflow. The 430-line `CONFIG_MANAGEMENT.md` only loads on demand.
3. **A3**: Hostname drift — `planning/mac/` (42 entries, canonical) and `planning/mac.lan/` (3 stray entries) caused by `hostname` (FQDN) vs `hostname -s` (short) at different times.

## Changes Planned

- [x] A1: Cleanup `.claude/settings.local.json` from 328 → ~75 entries
- [x] A2: Write root `CLAUDE.md` (~60 lines)
- [x] A3: `git mv` 3 entries from `planning/mac.lan/` → `planning/mac/`, remove empty dir

## Expected Outcome

- Permission file readable at a glance, no malformed entries
- Future sessions onboard from `CLAUDE.md` without reading 430-line CONFIG_MANAGEMENT
- Planning journal consolidated under one hostname

## Actual Outcome

**A1**: 328 → 83 entries (75% reduction).
- Dropped: 57 shell fragments (`Bash(do)`, `Bash(then)`, `Bash(while read file)`, etc.); 99 WebFetch domains; 7 `__NEW_LINE_*` synthetic markers; 18 specific for-loop fragments; 19 specific path tests; 7 `/tmp/skill-editor-session/*` paths; 2 hardcoded multi-line commit messages; assorted one-off echos.
- Consolidated to wildcards: `Bash(cp|mv|rm|awk|diff|shasum|pyright|curl:*)` replacing many specific entries; `Bash(./sync-config.py:*)` replacing 5 variants; conda subpath Reads subsumed under `Read(//Users/davidangelesalbores/**)`.
- Re-added: 8 curated WebFetch domains (docs.anthropic.com, github.com, arxiv.org, pmc.ncbi.nlm.nih.gov, www.nature.com, en.wikipedia.org, www.anthropic.com, gist.github.com).
- Backup retained at `.claude/settings.local.json.bak.20260512-094158`.
- Note: `.claude/` is gitignored, so this change is local to this machine only.

**A2**: `CLAUDE.md` written at repo root, 60 lines. Sections: core invariant, 7-step workflow, push-time gates, file layout table, planning conventions, sync exclusions, skill precedence notes, anti-patterns.

**A3**: 3 entries moved via `git mv` (history preserved), `planning/mac.lan/` removed. Canonical hostname is now `hostname -s` (`mac`), documented in CLAUDE.md.

## Assessment

**Result**: Success

**Improvements**:
- Permission file 75% smaller, fully grep-able, no parser artifacts
- Top-level CLAUDE.md available to all future sessions without explicit Read
- Single hostname dir, no future drift if CLAUDE.md convention is followed

**Issues**:
- None. Settings.local.json is local-only so no sync push needed for A1.

**Lessons Learned**:
- Permission accumulation appears to come from interactive grants of specific commands rather than wildcards. Future grants should prefer `Bash(<tool>:*)` over specific argument patterns.
- The `__NEW_LINE_*` artifacts suggest a Claude Code bug in permission serialization — worth filing if it recurs.

## Related Commits

- [pending]: tier-1 cleanup commit (settings.local.json local-only; commit covers CLAUDE.md + planning move)

## Next Steps

Proceed to Tier 2:
- B1: Deprecate `software-developer` skill (duplicate of `senior-developer`)
- B2: Deprecate `research-pipeline` skill (subsumed by `lit-pm`)
- B3: Resolve `code-review` vs `pr-review-toolkit` plugin overlap
