# Phase 0: Archival Guidelines Review — Full Detail

This reference holds the complete session-setup and archival-extraction
walkthrough for Phase 0 of the programming-pm workflow. SKILL.md keeps the phase
summary, quality gate, and a pointer here. Read this file when actively
executing Phase 0.

## Contents

- Process (session directory, stale-session cleanup, state storage)
- Primary Source: .archive-metadata.yaml
- Fallback: CLAUDE.md (Deprecated)
- Output
- Downstream Handoff
- Failure handling, session cleanup, and the Phase 0 → 1 handoff validation script

**Owner**: programming-pm (automatic)
**Checkpoint**: Never (always runs automatically)
**Duration**: 2-5 minutes
**Session Setup**: Creates `~/.claude/programming-pm-sessions/{workflow-id}/`

Initialize workflow session and extract archival guidelines, preferring `.archive-metadata.yaml` over CLAUDE.md, with code-specific extraction focus.

**Process**:
1. **Create session directory**:
```bash
mkdir -p ~/.claude/programming-pm-sessions/
SESSION_DIR="$HOME/.claude/programming-pm-sessions/prog-${PROJECT_SLUG}-$(date +%Y%m%d-%H%M%S)"
mkdir -p "${SESSION_DIR}"
```

2. **Stale session cleanup** (present to user, never auto-delete):
```bash
STALE_THRESHOLD=$((7 * 24 * 60 * 60))  # 7 days
CURRENT_TIME=$(date +%s)
STALE_FOUND=false

for SESSION in "$HOME/.claude/programming-pm-sessions"/*/; do
  [ ! -d "$SESSION" ] && continue
  STATE_FILE="${SESSION}state.yaml"
  if [ -f "$STATE_FILE" ]; then
    LAST_MODIFIED=$(stat -f %m "$STATE_FILE" 2>/dev/null || stat -c %Y "$STATE_FILE" 2>/dev/null || echo 0)
    AGE=$((CURRENT_TIME - LAST_MODIFIED))
    if [ "$AGE" -gt "$STALE_THRESHOLD" ]; then
      echo "STALE SESSION ($(( AGE / 86400 )) days old): $(basename "$SESSION")"
      STALE_FOUND=true
    fi
  fi
done

if [ "$STALE_FOUND" = true ]; then
  echo "Review stale sessions manually: ls $HOME/.claude/programming-pm-sessions/"
  echo "To remove: rm -rf \"$HOME/.claude/programming-pm-sessions/<session-id>\"  (verify non-empty session-id first)"
fi
```

3. **Store session path** in workflow state for downstream agents

## Primary Source: .archive-metadata.yaml
1. Follow the archival compliance check pattern:
   a. Read the reference document: `~/.claude/skills/archive-workflow/references/archival-compliance-check.md`
   b. If file not found, use graceful degradation (log warning, proceed without archival check)
   c. Apply the 5-step pattern to all file creation operations
2. Read `.archive-metadata.yaml` from the repo root
3. Extract code-specific guidelines:
   - `naming_conventions.project_specific_rules` for `*.py`, `*.js`, `*.ts` patterns
   - `structure.summary.source_code` and `structure.summary.tests`
   - `naming_conventions.summary.tests` (test file pattern)
   - `naming_conventions.summary.files` (general file naming)
4. Include the archival_context block in all downstream phase handoffs

## Fallback: CLAUDE.md (Deprecated)
If `.archive-metadata.yaml` is not found:
1. WARN: "Archival guidelines read from CLAUDE.md (fallback). Run archive-workflow
   to generate .archive-metadata.yaml for structured guidelines."
2. Check if `.archive-metadata.yaml` previously existed:
   - Look for `docs/organization/final-organization-report.md`
   - If found: WARN "Archival metadata was previously present but is now missing.
     Re-run archive-workflow."
3. Read CLAUDE.md and extract guidelines using existing prose extraction logic:
   - Code directory structure (`src/`, `modules/`, `experiments/`)
   - Git workflow (commit conventions, no destructive operations, stage specific files)
   - Testing conventions (if present)
   - Documentation conventions (README, inline comments, docstrings)
   - Repository organization for code vs. documentation
4. Produce archival-guidelines-summary.md as before

## Output
Write archival-guidelines-summary.md to the session directory with:
- Source: ".archive-metadata.yaml" or "CLAUDE.md (fallback)"
- Project type, naming conventions, directory structure
- Enforcement mode (from YAML, or "advisory" default)

```yaml
session_setup:
  session_dir: "~/.claude/programming-pm-sessions/{workflow-id}/"
  archival_summary_path: "{session_dir}/archival-guidelines-summary.md"
  guidelines_found: boolean
  guidelines_source: string  # ".archive-metadata.yaml" or "CLAUDE.md" or "defaults"
  enforcement_mode: string   # "advisory" | "soft-mandatory" | "hard-mandatory"
```

## Downstream Handoff
Include archival_context block in all agent dispatches (per the standard
archival context block defined in archival-compliance-check.md):

```yaml
archival_context:
  guidelines_present: true/false
  source: ".archive-metadata.yaml"  # or "CLAUDE.md" or "defaults"
  naming_convention: "snake_case"
  output_directory: "src/"
  enforcement_mode: "advisory"
  user_override: null
```

**Quality Gate**: Session directory created, archival summary written.

**Failure Handling**:
- `.archive-metadata.yaml` malformed: Treat as missing, fall back to CLAUDE.md
- CLAUDE.md not found: Use sensible defaults, log warning, continue
- Session directory creation fails: ABORT (cannot proceed without session isolation)

**Session Cleanup**:
- On successful completion (Phase 6 complete): Delete session directory
- On failure/abort: Retain session directory for debugging (log path to user)

**Timeout**: 5 min (ABORT on timeout - cannot proceed without session)

**Handoff Validation** (Phase 0 → Phase 1):
```bash
# Validate session handoff before proceeding to Phase 1
python3 "${SKILL_DIR}/scripts/validate-handoff.py" \
  "${SESSION_DIR}/handoffs/phase0-session-handoff.yaml" \
  "session_handoff"

if [ $? -ne 0 ]; then
  echo "❌ Phase 0 handoff validation FAILED"
  echo "Options:"
  echo "  (A) Fix issues in session handoff and retry"
  echo "  (B) Override with documented gaps (logged to session state)"
  read -p "Choice [A/b]: " OVERRIDE_CHOICE

  if [ "$OVERRIDE_CHOICE" != "b" ] && [ "$OVERRIDE_CHOICE" != "B" ]; then
    echo "Aborting. Fix session handoff and restart workflow."
    exit 1
  else
    echo "⚠️  Override: Proceeding with gaps (documented in session-state.json)"
    jq '.phase0_handoff_override = true | .phase0_handoff_gaps = "See validation errors above"' \
      "${SESSION_DIR}/session-state.json" > "${SESSION_DIR}/session-state.json.tmp"
    mv "${SESSION_DIR}/session-state.json.tmp" "${SESSION_DIR}/session-state.json"
  fi
else
  echo "✅ Phase 0 handoff validated successfully"
fi
```

**Phase Transition**: Phase 0 complete -> Quality Gate 0 -> PROCEED to Phase 1: Requirements and Scoping
