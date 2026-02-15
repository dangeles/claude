# Monthly Review Checklist

Step-by-step protocol for the orchestrator to follow during a full monthly
review. This file contains the pre-flight validation, session setup, quality
gates, deployment pipeline, and rollback procedures.

---

## Pre-Flight Validation Checklist (Phase 1, Step 0)

### Hard Requirements (abort if missing)

| Check | Command | Failure Action |
|-------|---------|----------------|
| git installed | `which git` | ABORT: "git is required. Install git and retry." |
| git push access per repo | `git ls-remote {repo-url}` for each site in registry | ABORT per-repo: "Cannot access [repo]. Check SSH keys or credentials. Offer: skip this repo, abort." |

**Limitation**: `git ls-remote` verifies read access, not write access. Push access is verified at push time in Phase 5. If a push fails due to insufficient permissions, it will be caught and handled by the Phase 5 error handling protocol.

| Session directory writable | `mkdir -p /tmp/web-presence-session/review-YYYY-MM/` | ABORT: "Cannot create session directory. Check /tmp permissions." |

### Soft Requirements (warn and continue)

| Check | Command | Failure Action |
|-------|---------|----------------|
| Jekyll installed | `which bundle && bundle exec jekyll --version` | WARN: "Jekyll not found. Build validation for Jekyll sites will be skipped. Pushes to Jekyll sites will proceed without build verification (user accepts risk)." |
| LaTeX installed | `which pdflatex` | WARN: "pdflatex not found. CV compilation validation will be skipped." |
| WebSearch available | (Cannot be pre-checked; SEO Manager handles unavailability gracefully per its instruction file) | N/A -- SEO Manager contains explicit WebSearch fallback protocol |
| Disk space | `df -h /tmp` (check available space) | WARN if < 1GB: "Low disk space in /tmp. Session may fail for large repos." |

---

## Session Setup Protocol (Phase 1)

### Step 1: Create Session Directory

```
/tmp/web-presence-session/review-YYYY-MM/
  repos/              -- cloned repositories go here
  outputs/            -- ALL sub-function reports written here:
    design-review.md    -- from Website Designer
    portfolio-review.md -- from Portfolio Manager
    seo-audit.md        -- from SEO Manager
    brand-reference.md  -- from Coherence Manager
    coherence-audit.md  -- from Coherence Manager
    action-items.md     -- from Suggestion Engine
    content-calendar.md -- from Suggestion Engine
    audit-report.md     -- from Suggestion Engine
  session-state.json  -- phase tracking and status
  rollback-info.json  -- pre-push SHAs (created in Phase 5)
  .lock               -- concurrent session prevention
```

### Step 2: Check for Existing Session

If `/tmp/web-presence-session/review-YYYY-MM/` already exists:

1. Read `session-state.json` to determine where the previous session stopped.
2. Present to user: "Found existing session from [date]. Status: Phase [N],
   [status]. Options: (a) Resume from Phase [N], (b) Start fresh (deletes
   existing session)."
3. If user chooses resume: load session state and continue from last completed
   phase.
4. If user chooses fresh: delete existing directory and start over.

### Step 3: Check for Lock File

If `.lock` file exists:

1. Read the lock file (JSON format: `{"locked_at": "ISO-8601", "session_id": "review-YYYY-MM"}`).
2. Check timestamp age against 120-minute global ceiling.
3. If lock is younger than 120 minutes: "Another review session is active (started [timestamp], session [ID]). Cannot run concurrently. Wait for it to finish or manually remove the lock file."
4. If lock is older than 120 minutes (stale): "Found stale lock file from [timestamp]. Previous session likely crashed or timed out. Removing lock and proceeding."

Create lock file:
```bash
echo '{"locked_at": "'$(date -u +%Y-%m-%dT%H:%M:%SZ)'", "session_id": "review-YYYY-MM"}' > .lock
```

### Step 4: Clone Repositories

For each site in the registry:

```bash
git clone --depth 25 --single-branch --branch {branch} \
  git@github.com:{repo}.git \
  /tmp/web-presence-session/review-YYYY-MM/repos/{repo-name}/
```

Use shallow clones (`--depth 25`) to save disk space while providing enough
history for Portfolio Manager's `git log --oneline -20` and rollback operations.
Record clone status in session-state.json per repo.

If a clone fails: record the error, continue with remaining repos. The
orchestrator will offer partial review or abort after all clones are attempted.

### Step 5: Load Previous Audit

Check for previous audit report in persistent location:
`~/.web-presence-audits/audit-YYYY-MM.md`

If found: copy to session directory as `previous-audit.md` for trend
comparison by the Suggestion Engine.

If not found: set `first_run: true` in session state. This is normal for the
first review.

### Step 6: Initialize Session State

Write `session-state.json` with the schema defined below.

---

## Session State Schema

```json
{
  "session_id": "review-YYYY-MM",
  "mode": "full-review",
  "started_at": "ISO-8601 timestamp",
  "current_phase": 1,
  "status": "in_progress",
  "environment": {
    "git": "version string",
    "ruby": "version string or null",
    "bundler": "version string or null",
    "jekyll_build_available": true,
    "latex_available": true
  },
  "sites": {
    "site-key": {
      "clone_status": "cloned | failed | skipped",
      "branch": "main",
      "pre_push_sha": "sha or null",
      "error": "error message or null"
    }
  },
  "phase_2": {
    "website_designer": {
      "status": "pending | running | complete | timeout | error",
      "started_at": "ISO-8601 or null",
      "completed_at": "ISO-8601 or null"
    },
    "portfolio_manager": { "status": "pending" },
    "seo_manager": { "status": "pending" }
  },
  "phase_3": {
    "coherence_manager": { "status": "pending" }
  },
  "phase_4": {
    "suggestion_engine": { "status": "pending" }
  },
  "phase_5": {
    "changes_approved": [],
    "push_progress": {
      "site-key": "pending | committed | pushed | failed | reverted"
    },
    "audit_saved": false
  },
  "previous_audit_path": "path or null",
  "first_run": true,
  "errors": []
}
```

Update session state after each phase transition and after each significant
event (clone success/failure, sub-function completion, push success/failure).

---

## Sub-Function Responsibility Boundaries

| Concern | Responsible | Consulted |
|---------|-------------|-----------|
| Visual design quality per site | Website Designer | -- |
| Accessibility per site | Website Designer | -- |
| Portfolio currency per site | Portfolio Manager | -- |
| Meta tags, structured data, SEO | SEO Manager | -- |
| Project description consistency | Portfolio Manager (detects) | Coherence Manager (resolves) |
| Color and font consistency | Coherence Manager | Website Designer |
| Cross-site bio alignment | Coherence Manager | -- |
| Cross-site timeline alignment | Coherence Manager | -- |
| Prioritization of all actions | Suggestion Engine | -- |
| Deployment and push operations | Orchestrator (SKILL.md) | -- |

---

## Inter-Phase Quality Gates

### Gate 1: Phase 1 -> Phase 2

Pass criteria:
- [ ] At least 1 repo cloned successfully (for partial review) OR all repos
      cloned (for full review)
- [ ] Site registry validated (all entries have required fields)
- [ ] Session state initialized and writable

Failure action: If zero repos cloned, ABORT. If partial: present list to user,
ask whether to proceed with partial review or abort.

### Gate 2: Phase 2 -> Phase 3

Pass criteria:
- [ ] At least 2 of 3 Phase 2 reports exist in `outputs/`
- [ ] Each existing report contains the required section headers:
  - `design-review.md`: "Current Assessment" AND "Accessibility Check"
  - `portfolio-review.md`: "Current State" AND "Updates Needed"
  - `seo-audit.md`: "Plugin Stack Status" AND "Critical Issues"
- [ ] No report consists solely of error messages or placeholder text

Failure action: If 1 report missing, proceed with warning passed to Coherence
Manager. If 2+ reports missing: pause and ask user whether to proceed (limited
coherence analysis) or re-run failed sub-functions.

### Gate 3: Phase 3 -> Phase 4

Pass criteria:
- [ ] `coherence-audit.md` exists in session directory
- [ ] Contains both "Narrative Coherence" and "Visual Coherence" sections
- [ ] Contains at least one score

Failure action: If missing or incomplete, skip coherence input to Suggestion
Engine. Flag in session state. Suggestion Engine will produce output with
reduced confidence for coherence-related items.

### Gate 4: Phase 4 -> Phase 5

Pass criteria:
- [ ] `action-items.md` exists with non-trivial content (> 100 words)
- [ ] `content-calendar.md` exists with non-trivial content (> 100 words)
- [ ] `audit-report.md` exists with non-trivial content (> 200 words)

Failure action: If any deliverable missing, present raw Phase 2/3 outputs to
user directly. Flag the missing synthesis.

### Gate 5: Before Each Push (Phase 5)

Pass criteria:
- [ ] Build validation passed for this repo (or build_validation is `none`)
- [ ] Rollback info recorded in `rollback-info.json` (pre-push SHA saved)
- [ ] User explicitly approved changes for this specific repo
- [ ] Commit uses conventional format with HEREDOC message
- [ ] Changes staged specifically (NOT `git add -A` or `git add .`)

Failure action: Block push for this repo. Report the failing check. Offer:
fix the issue, skip this repo, or abort remaining pushes.

### Output Contract Map

Traces data flow between producers and consumers. Gate checks validate required sections.

| Producer | Output File | Required Sections | Consumer(s) |
|----------|------------|-------------------|-------------|
| Website Designer | `outputs/design-review.md` | Current Assessment, Accessibility Check, Quick Wins, Recommendations | Coherence Manager, Suggestion Engine |
| Portfolio Manager | `outputs/portfolio-review.md` | Current State, Updates Needed, Consistency Check | Coherence Manager, Suggestion Engine |
| SEO Manager | `outputs/seo-audit.md` | Plugin Stack Status, Critical Issues, Structured Data Check, Canonical URL Check | Suggestion Engine |
| Coherence Manager | `outputs/coherence-audit.md` | Narrative Coherence Score, Visual Coherence Score | Suggestion Engine |
| Coherence Manager | `outputs/brand-reference.md` | Visual Identity, Narrative Identity | (reference only) |
| Suggestion Engine | `outputs/action-items.md` | Must Do (this month), Should Do, Backlog | Orchestrator (Phase 5) |
| Suggestion Engine | `outputs/audit-report.md` | Executive Summary, Scores, Priority Matrix | Orchestrator (Phase 5) |
| Suggestion Engine | `outputs/content-calendar.md` | Blog Post Ideas, Portfolio Updates | Orchestrator (Phase 5) |

---

## Phase 5 Deployment Pipeline

### Step 5a: Present Audit Summary

Read `audit-report.md` Executive Summary and Scores table. Present to user
as the opening of Phase 5.

### Step 5b: Present Action Items for Approval

Read `action-items.md` "Must Do" list. Ask user which items to implement now.
Present as a numbered list for easy selection.

### Step 5c: Group Changes by Repository

For each approved action item, determine which repository it affects. Group
changes by repo for efficient per-repo deployment.

### Step 5d: Per-Repo Deployment

For each repo with approved changes:

1. **Make edits** in the cloned repo working copy.
2. **Stage specific files** -- NEVER use `git add -A` or `git add .`. Stage
   only the files that were modified for approved changes.
3. **Commit with HEREDOC message**:
   ```bash
   git commit -m "$(cat <<'EOF'
   <type>(<scope>): <description>

   Reviewed-by: web-presence-manager monthly audit (YYYY-MM)
   EOF
   )"
   ```
4. **Run build validation** using the `build_command` from site registry for each site:
   - If `build_validation` is `required`: build MUST pass. Block push on failure.
   - If `build_validation` is `optional`: build failure produces a WARNING but push proceeds.
   - If `build_validation` is `none` or `build_command` is `none`: skip validation.
   - For LaTeX: consider using `latexmk -pdf {file}` instead of raw `pdflatex` for proper multi-pass compilation.
5. **Show diff to user**: `git diff HEAD~1` for this repo.
6. **Get per-repo push confirmation**: "Push these changes to [repo]? (y/n)"
7. **Push**: `git push origin {branch}`
8. **Record status** in session state.

### Step 5e: Post-Push Verification

After all pushes complete:
- Report which repos were pushed successfully
- Report any failures
- Verify cross-repo consistency (if coherence changes were made)

### Step 5f: Save Audit Report

Save audit to persistent location:

```bash
mkdir -p ~/.web-presence-audits/
cp /tmp/web-presence-session/review-YYYY-MM/outputs/audit-report.md \
   ~/.web-presence-audits/audit-YYYY-MM.md
```

Note: The `mkdir -p` is essential on first run -- without it, the copy fails silently and trend tracking is permanently broken.

---

## Two-Phase Commit Protocol for Multi-Repo Push

When pushing to multiple repos, use this protocol to maintain consistency:

### Phase A: Commit All Locally

1. For each repo with approved changes, create the commit locally.
2. Do NOT push yet.
3. Verify all commits are valid (no empty commits, correct files staged).
4. Record pre-push SHA for each repo in `rollback-info.json`.

### Phase B: Push Sequentially

1. Push repos one at a time in this order: primary site first, then secondary
   sites in registry order.
2. After each push, record success/failure in session state.
3. If a push fails: STOP. Report which repos succeeded and which failed.

### On Partial Push Failure

Present user with options:
1. **Retry**: Attempt the failed push again.
2. **Revert successful pushes**: Roll back repos that were already pushed
   (see Rollback Protocol).
3. **Accept inconsistency**: Leave successful pushes in place, report the
   failed repos for manual resolution.

---

## Rollback Protocol

### Before Push: Record State

For each repo about to be pushed, record in `rollback-info.json`:

```json
{
  "repos": {
    "dangeles.github.io": {
      "pre_push_sha": "abc123...",
      "branch": "main",
      "push_status": "pending"
    }
  },
  "created_at": "ISO-8601"
}
```

### Standard Rollback (safe)

```bash
git revert HEAD --no-edit
git push origin main
```

This creates a new revert commit, preserving history. Preferred method.

### Force Rollback (last resort, user must explicitly approve)

```bash
git reset --hard {pre_push_sha}
git push --force-with-lease origin main
```

Only use if the revert approach fails or if the commit introduced a build
break that cannot be reverted cleanly. Requires explicit user confirmation:
"This will force-push to [repo], rewriting history. Confirm? (y/n)"

---

## Commit Message Format

All commits created by the web-presence-manager use conventional commit format:

```
<type>(<scope>): <description>

Reviewed-by: web-presence-manager monthly audit (YYYY-MM)
```

### Types

| Type | Use When |
|------|----------|
| `fix` | Fixing broken links, correcting dates, fixing HTML errors |
| `feat` | Adding new portfolio entries, new blog posts, new sections |
| `style` | CSS changes, design improvements, formatting |
| `content` | Updating bio, descriptions, text content |
| `seo` | Meta tags, structured data, sitemap, robots.txt |
| `docs` | README updates, documentation improvements |

**Note**: The `content` and `seo` types are custom additions for this workflow. If the target repository uses commitlint or similar tools that enforce standard conventional commit types, use `chore(content)` or `chore(seo)` as fallback types.

### Scopes

| Scope | Use When |
|-------|----------|
| `blog` | Changes to the Jekyll blog/personal site |
| `profile` | Changes to the GitHub profile README |
| `cv` | Changes to the LaTeX CV |
| `all` | Changes that affect multiple sites |

---

## Session Cleanup

### At Session Completion

Offer the user two options:
1. **Keep audit reports only**: Delete cloned repos from `/tmp/`, keep
   `audit-report.md`, `action-items.md`, `content-calendar.md` in session
   directory. Copy audit to persistent location.
2. **Keep everything**: Leave session directory intact for reference.

### At Session Start (housekeeping)

Check for old sessions in `/tmp/web-presence-session/`:
- Sessions > 60 days old: "Found old review session from [date]. Delete? (y/n)"
- Remove stale lock files (process not running).
