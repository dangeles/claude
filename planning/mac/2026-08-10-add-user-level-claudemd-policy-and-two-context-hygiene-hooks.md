# Add user-level CLAUDE.md policy and two context-hygiene hooks

**Date**: 2026-08-10
**Machine**: mac
**Status**: Success

## Objective

Measured working style, not opinion: across 311 sessions and 34,889 requests, context
averaged **163,573 tokens re-read per request** (worst sessions sustain 400-430K across
700+ requests), against only 852 output tokens per request. Cache reads are 5.71B tokens
— ~88% of cost. Work stays in the main loop, so context grows monotonically and every
later turn pays for it.

Root cause found while investigating why skills never fire: **there was no user-level
`~/.claude/CLAUDE.md`**. Every session in every repo got the harness prompt plus that
repo's project CLAUDE.md, and nothing else. All 58 skills, 27 agents, and
`skills/references/delegation-and-scope.md` are invisible to the main loop unless
something invokes them — and `delegation-and-scope.md` is read only by orchestrator
skills, which fired 0 times in 311 sessions. It was dead code.

The mechanism that demonstrably works was already visible in the data: superpowers skills
account for **29 of 75 total invocations (39%)** from ~9 of ~100 available skills, because
a SessionStart hook injects their entry point into every session. Always-loaded context
changes behaviour; library content awaiting invocation does not.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Add `claude-config/CLAUDE.md` — short, always-loaded working policy
- [x] Add `CLAUDE.md` to `sync.config.yaml` sync_rules so it is version-controlled
- [x] Add `read-burst-nudge.py` (PostToolUse: Read|Grep|Glob)
- [x] Add `context-budget-warn.py` (Stop)
- [x] Register both in `settings.json`
- [x] Tests for both hooks, matching the `test_paper_context_guard.py` pattern

## Expected Outcome

Delegation and targeted reading become the default in the main loop, cutting the context
size that gets re-read on every request. The two hooks enforce at the point of action what
the policy states in prose.

## Actual Outcome

**`claude-config/CLAUDE.md`** (~55 lines, policy only since it is re-read on every
request). Five sections: read targeted not whole (`offset`/`limit`, `Grep -A/-B`, never
re-read, when whole-file reads *are* right); delegate to protect context (with the
threshold list); pass `model` on every dispatch; plan before 3+ file work; budget past
150K. States plainly that project CLAUDE.md and direct instructions override it.

**`read-burst-nudge.py`** — fires once per session at 15 main-thread reads with no
subagent dispatched. Perf contract matters here: PostToolUse on Read fired 3,331 times
historically, so the hot path is an O(1) counter in a JSON state file and the transcript
is walked **at most once per session**, at the threshold. Subagent-exempt (a subagent
reading many files is the recommended pattern). `noqa-read-burst` bypass.

**`context-budget-warn.py`** — Stop hook, bands at 150K/250K/350K, once per band, no
re-fire when context shrinks. Reads the last `cache_read_input_tokens` from the
transcript, which is the prefix the API actually re-read — the real number, not an
estimate. Tail-reads the last 256KB and scans backwards, so it stays cheap on the large
transcripts it is most relevant to. Mute per-session via `{"muted": true}` in state.

`sync.config.yaml` gained a `CLAUDE.md` sync rule; without it the new file would have
been exactly the untracked drift this repo exists to prevent.

### Verification

- 106 tests pass, 1 skipped (was 80/1) — 26 new across the two hooks
- `validate` passes (frontmatter + JSON)
- `context-budget-warn.py` run against a **real 3.6M transcript**: correctly reported
  307K and selected the 250K band. Parsing validated against production data, not fixtures.
- `read-burst-nudge.py` end-to-end: silent for 14 reads, fired on the 15th, state
  `{"reads": 15, "done": true}`, silent after
- `status` reports no repo/live divergence; exec bits preserved on both hooks

## Assessment

**Result**: Success

**Improvements**:
- The main loop now receives delegation, targeted-reading, and model-on-dispatch policy in
  every session — the first time any of this guidance is reachable where decisions happen.
- Two hooks enforce at the point of action, the pattern that already yields 100% frontmatter
  compliance versus 15% skill invocation.
- `context-budget-warn` makes the dominant cost term visible for the first time.

**Issues**:
- The policy file is a nudge, not enforcement. The hooks are what make it stick.
- `read-burst-nudge`'s threshold of 15 is a first guess. It may prove noisy or too lax;
  it is a single constant (`THRESHOLD`) to retune.
- Neither hook can force a dispatch — they inform. If the working style does not change,
  the next step is a harder gate, not a louder message.
- Effect is unmeasured. `stats-cache.json` was 12 days stale, so the 163K/request figure
  from transcripts is the baseline to compare against.

**Lessons Learned**:
- Look for the *absent* artifact. Four turns went into why skills do not fire before
  noticing that the one always-loaded user-level file did not exist. The missing file
  explained more than any of the present ones.
- Guidance in a file nothing loads is not guidance. Check the load path before improving
  the prose.
- Validate transcript parsing against a real transcript. Fixtures agreed with my
  assumptions; only production data proves the format.

## Related Commits

- [pending]: feat(config): add user-level CLAUDE.md policy and context-hygiene hooks

## Next Steps

1. Run a week of normal sessions, then re-run the transcript token report. The metric that
   matters is cache-read-per-request falling from 163,573.
2. Retune `THRESHOLD` in `read-burst-nudge.py` from experience — raise if noisy, lower if
   sessions still reach 300K without it firing.
3. Add `tools/token-report.py` to the repo so the measurement is repeatable rather than
   ad-hoc (the script used for this analysis is not yet committed anywhere).
4. Only then triage the 29 orchestrator-addressed skills. Rewriting them for interactive
   use is the largest remaining item, and it should follow evidence that the working style
   has actually shifted.
5. Still open from the prior entry: `effortLevel: "high"` unchanged pending fresh data
   (output is only ~12% of cost, so this is a minor lever); `skills/README.md` stale counts.
