# Fix sync-config hostname detection bug

**Date**: 2026-05-24
**Machine**: mac
**Status**: Success

## Objective

Fix the hostname-detection bug in `sync-config.py` that caused every `./sync-config.py plan --title "..."` invocation to create planning entries under `planning/192/` instead of the canonical `planning/mac/`. Worked around manually in 4 previous planning entries (plotting-advisor, tier-3 extension, tier-1+2 audit, and this one's predecessor) — due for a real fix.

## Root cause

On this Mac, the OS-level hostname is misconfigured: `hostname`, `hostname -s`, `hostname -f`, and `platform.node()` all return `192.168.1.9` (the LAN IP), not a proper short name. `sync-config.py`'s `detect_hostname()` did `platform.node().split('.')[0]`, which produced `"192"`. That directory name was then used as the canonical machine directory.

## Changes Planned

- [x] Follow CONFIG_MANAGEMENT.md workflow
- [x] Edit `sync-config.py` to add a config-driven hostname resolution with IP-pattern validation
- [x] Add `machine_name: "mac"` to `sync.config.yaml`
- [x] Test all three paths: YAML config, env-var override, validation error
- [x] Verify `./sync-config.py plan` now creates entries under `planning/mac/` on the first try

## Implementation

Rewrote `ClaudeCodeSync.detect_hostname()` with this resolution order:

1. **`SYNC_CONFIG_MACHINE` env var** — per-invocation override; useful for CI or for testing
2. **`machine_name` field in `sync.config.yaml`** — per-machine canonical, the typical fix
3. **`platform.node()` short form** — legacy fallback for correctly-configured hosts

After resolution, validates the result against the regex `\d+(\.\d+){0,3}` — if the candidate is purely numeric or IP-like, raises a clear `ValueError` instructing the user to set `machine_name` or `SYNC_CONFIG_MACHINE`. This makes future OS misconfigurations fail loud rather than silently writing to a non-canonical directory.

Added `machine_name: "mac"` to `sync.config.yaml` (top, near `version` and `source_dir`) with an inline comment explaining the field and the env-var alternative.

## Expected Outcome

- `./sync-config.py plan --title "X"` creates the entry directly under `planning/mac/` (no manual move).
- If a future user has a working `hostname` and no `machine_name` set, the legacy `platform.node()` fallback still works for them.
- If a future user hits the same OS misconfiguration, they get an actionable error message instead of a silent `planning/192/` directory.

## Actual Outcome

All three paths verified end-to-end:

1. **Happy path (YAML config)**: `./sync-config.py plan --title "Fix sync-config hostname detection bug"` created `planning/mac/2026-05-24-fix-sync-config-hostname-detection-bug.md` directly — no manual `mv` needed for the first time.
2. **Env override**: `SYNC_CONFIG_MACHINE=test-env-host ./sync-config.py plan --title "env override test"` correctly created the entry under `planning/test-env-host/`.
3. **Negative path**: With `machine_name` commented out, the tool exits with `Error: Detected machine name '192' looks like an IP address or numeric token. The OS hostname is likely misconfigured. Set 'machine_name' in sync.config.yaml (e.g. 'machine_name: mac') or export SYNC_CONFIG_MACHINE=<short-name>.` — clear, actionable.

## Assessment

**Result**: Success

**Improvements**:
- Planning entries now land in the canonical location on the first invocation.
- Future OS hostname misconfigurations fail loud with an actionable error rather than silently writing to a wrong directory.
- The env-var escape hatch is useful for CI/test contexts where modifying `sync.config.yaml` isn't ideal.
- The config-driven approach also accommodates users whose canonical short name doesn't match `hostname -s` (e.g., a personal naming convention vs. the OS's idea of the hostname).

**Issues**:
- None observed.

**Lessons Learned**:
- "Use `hostname -s`" advice in `CLAUDE.md` doesn't help when the OS hostname is itself broken. A config-driven canonical name is the more robust pattern.
- Pyright surfaced 6 unrelated diagnostics in `sync-config.py` (unused imports, unreachable code, unused variables in unrelated functions) — those are pre-existing and out of scope for this fix; worth a separate cleanup pass if metadata hygiene comes back up.

## Related Commits

- `5f87155`: fix(sync-config): config-driven hostname resolution with IP-pattern validation

## Next Steps

- Optional follow-up: clean up the 6 pre-existing pyright diagnostics in `sync-config.py` (unused `Dict` import, unused for-loop variables, unreachable code branch). Pure hygiene; no behavior change.
- Optional follow-up: skills with no `last_updated` field — set the field consistently across the catalog for metadata hygiene.
