# Session State Schema

## Session State YAML Schema

The session state is stored in `{session_dir}/session-state.yaml` and tracks all pipeline progress.

```yaml
# Session Identity
pipeline_id: "essay-{uuid}"
started_at: "{ISO8601}"
last_activity: "{ISO8601}"
status: "in_progress"  # in_progress | paused | completed | abandoned

# Configuration
style_profile_path: "{path}"
samples_dir_path: "{path}"
pushback_level: "full"  # full | light | minimal
session_dir: "{path}"

# Pipeline Progress
current_stage: 1  # 1-4
current_stage_name: "thesis_development"  # thesis_development | essay_structuring | argument_development | paragraph_writing

# Stage 1: Thesis Development
stages:
  stage_1:
    status: "not_started"  # not_started | in_progress | completed | needs_review
    thesis_file: null  # "{session_dir}/stage-1-thesis.md" when completed
    development_log: null  # "{session_dir}/stage-1-development-log.md"
    claim_type: null  # argument | explainer | prediction | hybrid
    completed_at: null  # "{ISO8601}"

  # Stage 2: Essay Structuring
  stage_2:
    status: "not_started"
    outline_file: null  # "{session_dir}/stage-2-outline.md"
    total_sections: null
    target_word_count: null
    audience: null  # general | semi-technical | expert
    completed_at: null

  # Stage 3: Argument Development (per-section)
  stage_3:
    status: "not_started"
    total_sections: null
    current_section: null
    sections:
      # Populated dynamically based on outline
      # - section_index: 1
      #   status: "not_started"
      #   argument_map_file: null
      #   verified_claims: 0
      #   deferred_claims: 0
      #   completed_at: null
    birds_eye_approved: false
    completed_at: null

  # Stage 4: Paragraph Writing (per-paragraph)
  stage_4:
    status: "not_started"
    current_section: null
    current_paragraph: null
    draft_file: null  # "{session_dir}/stage-4-draft.md"
    paragraph_log: null  # "{session_dir}/stage-4-paragraph-log.md"
    sections:
      # Populated dynamically
      # - section_index: 1
      #   total_paragraphs: null
      #   completed_paragraphs: 0
      #   status: "not_started"
      #   voice_check_score: null
    completed_at: null

# Cross-Cutting State
voice_calibration:
  first_assessment_user_rating: null  # 1-5
  calibration_notes: null
  voice_reference_mode: "profile"  # profile | sample | hybrid

deferred_verifications: []
  # - claim_text: "[text]"
  #   section: N
  #   paragraph: null
  #   reason: "[timeout | service_unavailable | no_results]"
  #   added_at: "{ISO8601}"
  #   resolved: false

user_overrides: []
  # - claim_text: "[text]"
  #   fact_checker_result: "[status]"
  #   user_action: "override"
  #   user_rationale: "[text]"
  #   logged_at: "{ISO8601}"

revision_history: []
  # - stage_revisited: 1
  #   impact: "major"  # minor | moderate | major
  #   reason: "[user request description]"
  #   downstream_archived: ["stage-2-outline.v1.md", "stage-3-section-1-arguments.v1.md"]
  #   timestamp: "{ISO8601}"

# Output
output:
  final_essay_file: null  # "{session_dir}/final-essay.md"
  sources_verified: null
  voice_consistency_score: null
  total_word_count: null
  completed_at: null
```

## Session Directory Structure

```
{session_dir}/
  session-state.yaml          # Primary state file
  session-state.yaml.bak      # Backup (previous version)
  session-state.yaml.tmp      # Temporary file for atomic writes
  stage-1-thesis.md            # Stage 1 output
  stage-1-development-log.md   # Stage 1 conversation log
  stage-2-outline.md           # Stage 2 output
  stage-3-section-1-arguments.md   # Stage 3 per-section outputs
  stage-3-section-2-arguments.md
  stage-3-section-N-arguments.md
  stage-4-draft.md             # Stage 4 incremental draft
  stage-4-paragraph-log.md     # Stage 4 paragraph journey log
  final-essay.md               # Final assembled essay
  fact-check-log.md            # All verification results and overrides
  voice-check-log.md           # All voice assessments
```

## Atomic Write Protocol

All session state writes must be atomic to prevent corruption:

```bash
# 1. Write to temporary file
echo "$content" > "{session_dir}/session-state.yaml.tmp"

# 2. Backup current state
cp "{session_dir}/session-state.yaml" "{session_dir}/session-state.yaml.bak"

# 3. Atomic rename (replaces old file)
mv "{session_dir}/session-state.yaml.tmp" "{session_dir}/session-state.yaml"
```

The orchestrator performs this via Bash tool. Never write directly to `session-state.yaml`.

## Backup Protocol

- Before every state update, copy the current state to `.bak`
- The `.bak` file always contains the previous valid state
- Stage output files (thesis, outline, argument maps, draft) serve as secondary backups of their respective stage data

## Recovery Protocol

If `session-state.yaml` cannot be parsed (YAML error, corruption):

### Level 1: Backup Recovery
```
1. Read session-state.yaml.bak
2. Parse YAML
3. If valid: Replace session-state.yaml with .bak content
4. Report: "Session state recovered from backup. Last state update may be lost."
```

### Level 2: Artifact Reconstruction
If both primary and backup are corrupt:
```
1. Scan session directory for artifact files
2. Reconstruct state from artifacts:
   - stage-1-thesis.md exists → stage_1.status = "completed"
   - stage-2-outline.md exists → stage_2.status = "completed"
   - stage-3-section-N-arguments.md files → stage_3 section statuses
   - stage-4-draft.md exists → determine progress from content
3. Reconstruct minimal state from artifact evidence
4. Report: "Session state reconstructed from artifact files. Some metadata may be approximate."
```

### Level 3: Manual Recovery
If artifact reconstruction fails:
```
1. Report to user: "Session state could not be automatically recovered."
2. List available artifact files
3. Ask user where to resume
4. Create fresh state based on user direction
```

## State Validation Rules

When loading session state (on resume), validate:

| Condition | Check |
|-----------|-------|
| Stage marked "completed" | Verify artifact file exists at specified path |
| Section marked "completed" | Verify argument map file exists |
| current_section value | Must be within range 1 to total_sections |
| current_paragraph value | Must be within reasonable range |
| File paths | All referenced files must exist on disk |
| Status consistency | No stage can be "completed" if a prerequisite stage is "not_started" |

If validation fails, report the discrepancy and ask the user how to proceed.

## Pause/Resume Protocol

### Pause
1. Save any in-progress work (draft paragraphs, partial argument maps) to files
2. Update session state with current position
3. Atomic write to session-state.yaml
4. Display: "Session paused at [Stage N, Section M, Paragraph P]. To resume, invoke essay-pipeline and select this session."

### Resume
1. Scan for session directories matching the pattern `/tmp/essay-pipeline-*/`
2. Display available sessions with summary:
   ```
   Sessions found:
   1. essay-pipeline-20260214-143022 -- Stage 3, Section 2 -- Last active: 2 hours ago
   2. essay-pipeline-20260212-091500 -- Stage 4, Section 1 -- Last active: 2 days ago
   ```
3. User selects a session
4. Load and validate state
5. Present summary: "Resuming at [Stage N]. Here's where you left off: [brief summary]"
6. Continue from the last checkpoint

### Session Expiration
- **Inactivity timeout**: 24 hours -- Auto-pause
- **Abandonment timeout**: 7 days -- Warn on next access; offer cleanup or resume

## Artifact Versioning Convention

When a stage is revisited and its outputs change, version the old artifacts:

```
stage-1-thesis.md       → stage-1-thesis.v1.md (archived)
stage-1-thesis.md       → (new version, current)
```

Versioning format: `{filename}.v{N}.md` where N increments.

Archived versions are kept for reference but are not used by the pipeline. Only the unversioned file is the current version.

## Change Impact Assessment

When the user revisits an earlier stage, assess the impact on downstream work:

### Minor Impact
- Change does not affect the substance of downstream stages
- Example: Rewording the thesis without changing its meaning
- Action: Flag downstream artifacts for review but do not archive

### Moderate Impact
- Change affects some but not all downstream stages
- Example: Modifying one section in the outline
- Action: Archive affected downstream artifacts; the user reviews and decides what to revise

### Major Impact
- Change fundamentally alters the essay direction
- Example: Changing the thesis from an argument to an explainer
- Action: Archive ALL downstream artifacts; warn the user about scope; restart from the revised stage

### Assessment Process

```
1. Identify what changed (diff between old and new stage output)
2. Classify impact: minor / moderate / major
3. If minor: Flag downstream artifacts with "needs_review" status
4. If moderate: Archive affected artifacts; set their stages to "needs_review"
5. If major: Archive all downstream artifacts; set all subsequent stages to "not_started"
6. Log in revision_history
7. Present impact to user and confirm before proceeding
```
