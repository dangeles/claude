# Sharpen git-strategy-advisor + add session-scope discipline

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

Two related changes to `git-strategy-advisor`:

1. **Triggering**: rewrite the description so it actually fires on the cases it's good at (non-trivial scope decisions, integration crossroads) and stays out of the way on routine commits. Add explicit "NOT for" clauses naming the four sibling skills it currently collides with at the description level.
2. **Behavior**: add session-scope discipline — when generating recommendations, stage ONLY files this session modified and flag every other dirty file for separate review. Prevents inter-agent collisions where one agent's incomplete work gets bundled into another agent's commit.

Surfaced from a session-end audit. Over ~30 commits this session, the skill never fired once — the model defaulted to raw `git commit` via Bash. The description was too broad to compete with that default.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Rewrite the description with sharpened triggers + NOT-for clauses (`/commit-commands:commit`, `/commit-commands:commit-push-pr`, `superpowers:finishing-a-development-branch`, `superpowers:using-git-worktrees`)
- [x] Add a "Session-Scope Discipline" section to the body (between "When NOT to Use" and "Workflow Overview")
- [x] Update Phase 2 metrics extraction to filter to session-set files (out-of-session files counted separately as warnings, not as scope inputs)
- [x] Add `OUT_OF_SESSION_CHANGES` and `STAGED_OUT_OF_SESSION` warning examples to the Output Schema
- [x] Validate frontmatter + sync push + verify live state

## Implementation

**Description** went from a 230-character one-liner to a 798-character paragraph that front-loads the trigger phrases ("how should I structure this?", "branch or commit?", "ready to ship"), explicitly names the session-scope behavior, and routes four sibling skills via NOT-for clauses.

**Session-Scope Discipline** section (60 lines) frames how to operate:

- *Identifying the session set*: three sources in priority order — explicit caller declaration (`session_files: [...]`), conversational scan of recent tool calls (Write/Edit/NotebookEdit/Bash mutations), fall back to asking the user.
- *Partitioning the working tree*: split `git status --porcelain` into "in session set" (commit candidates) and "out of session" (warning candidates).
- *Recommendation discipline*: stage commands MUST list paths explicitly (never `-A` / `.` / `-u`); empty session set with dirty tree → review-only recommendation, no commit/push/PR; staged out-of-session paths are an explicit red flag (`STAGED_OUT_OF_SESSION` warning with `git reset HEAD <path>` remediation).
- *Edge cases*: untracked files, deletes, renames, Read-only files, caller-passed empty session-files list.

**Phase 2.3 (Metrics Extraction)** updated to clarify that scope metrics (files / lines / dirs) are computed over the session set only; out-of-session files don't influence scope classification, they get the warning treatment.

**Output Schema** gained two warning examples (`OUT_OF_SESSION_CHANGES`, `STAGED_OUT_OF_SESSION`) with the `paths: [{path, status}, ...]` shape.

## Expected Outcome

- Triggering: the skill fires when actually useful (non-trivial work, integration crossroads), and routes elsewhere for routine commits / execution / final integration / worktree setup.
- Safety: in multi-agent flows, the skill is structurally incapable of recommending a commit that bundles another agent's dirty files — because the bare `git add -A` / `git add .` patterns are explicitly forbidden, and out-of-session paths surface as warnings rather than scope inputs.
- Verifiable via the live available-skills index reflecting the new description.

## Actual Outcome

All shipped cleanly. The new description is live in `~/.claude/skills/git-strategy-advisor/SKILL.md` and appears in the available-skills index with the session-scope language and all four NOT-for clauses intact. SKILL.md grew from 349 → 405 lines (+56 lines, the Session-Scope Discipline section plus the Phase 2 + Output Schema touch-ups).

## Assessment

**Result**: Success

**Improvements**:
- The description's NOT-for trio for `/commit-commands:commit`, `superpowers:finishing-a-development-branch`, and `superpowers:using-git-worktrees` is the first time these four skills have been cross-routed at the description level. Continues the May 12 tier-3 pattern (`56a667e`) into the git ecosystem.
- The session-scope discipline is a reusable pattern — any future skill that mutates git state could adopt the same `session_files` declaration / partition / explicit-stage-only contract. Worth promoting to a shared reference if more skills join the git lane.
- The skill is still advisory-only (recommendations + warnings, never commits). The discipline added is about WHAT to recommend, not about giving the skill write access.

**Issues**:
- Description length (798 chars) is on the long side. The trade-off is intentional — disambiguating across four siblings inside the always-loaded description field is high-leverage. If it shows up as noisy in future audits, the inline trigger phrases ("how should I structure this?", etc.) are the first candidate to trim.

**Lessons Learned**:
- When a skill "under-fires" relative to its intended use cases, the symptom is usually that the description's trigger surface doesn't compete with the model's default behavior (raw Bash in this case). Sharpening triggers + naming the default-behavior alternative as NOT-for is a stronger fix than adding more body content.
- For session-scope discipline, the "session_files explicit declaration" path is the most reliable; the "conversational scan" path is best-effort. Skills that need this discipline should make the explicit declaration the documented contract, with the scan as a fallback.

## Related Commits

- `d70e57e`: chore(git-strategy-advisor): sharpen description + add session-scope discipline

## Next Steps

- If `weekly-digest`, `meeting-notes-to-action`, or any of the future Notion/Slack-integrated skills mutate git state, adopt the same `session_files` declaration contract documented here.
- Consider promoting the Session-Scope Discipline section to a shared reference (`claude-config/skills/_shared/session-scope.md`) if a second skill needs the same pattern.
