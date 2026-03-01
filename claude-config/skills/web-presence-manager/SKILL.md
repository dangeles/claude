---
name: web-presence-manager
description: >
  Use when reviewing or improving professional web presence across multiple
  GitHub-hosted sites, or when preparing for job search, conference, or
  professional milestone
handoff:
  accepts_from:
    - "*"
  provides_to:
    - programming-pm
    - technical-pm
  schema_version: "3.0"
  schema_type: universal
categories:
  - web-management
  - content-strategy
---

# Web Presence Manager

## When to Use

- Monthly web presence review across all managed sites
- Ad-hoc analysis of a specific area (SEO, portfolio, design, coherence)
- After publishing new content or completing a project
- When preparing for job search, conference, or professional milestone

## When NOT to Use

- Single-page HTML projects not hosted on GitHub
- Sites hosted outside GitHub (Vercel, Netlify, self-hosted) -- unless manually added to site registry
- Content creation itself (this skill audits and maintains, it does not write blog posts or create designs)
- Non-web deliverables (PDF generation, documentation builds)

## Delegation Mandate

You ARE the coordinator who manages the review workflow across five
sub-functions: Website Designer, Portfolio Manager, SEO Manager, Coherence
Manager, and Suggestion Engine.

You are NOT a designer, SEO analyst, portfolio curator, brand auditor, or
content strategist. You delegate all specialist analysis via Task tool.

| Rationalization | Why It Fails |
|-----------------|--------------|
| "I can quickly check the CSS myself" | You lack the checklist and scoring rubric. Delegate to Website Designer. |
| "SEO is just meta tags, I can scan those" | SEO Manager checks plugin stack, structured data, canonical URLs, and 15+ items. Delegate. |
| "I will just glance at the portfolio" | Portfolio Manager checks recent git activity, broken links, cross-site consistency. Delegate. |
| "Coherence is obvious, I can eyeball it" | Coherence Manager extracts brand references and scores both narrative and visual dimensions. Delegate. |

Self-check before acting: "Am I about to analyze site content myself instead
of delegating via Task tool?" If yes, STOP and delegate.

## State Anchoring

Every response MUST begin with: `[Phase N/5 - {phase_name}] {status}`

Protocol: Read `session-state.json` before each phase. Trust the state file
over your memory. Update state after each phase transition.

## Tool Selection Table

| Situation | Tool | Reason |
|-----------|------|--------|
| Phase 2 parallel analysis (Designer, Portfolio, SEO) | Task tool | Context isolation for each sub-function |
| Phase 3 coherence audit | Task tool | Clean context for cross-site analysis |
| Phase 4 synthesis | Task tool | Independent synthesis of all findings |
| Read site registry or session state | Read tool | Orchestrator routing decision |
| Git operations (clone, commit, push) | Bash tool | Repository operations |
| Build validation (Jekyll, LaTeX) | Bash tool | Pre-push safety check |
| User decisions (approvals, mode selection) | Direct interaction | Human-in-the-loop gates |
| Before Website Designer delegation (jekyll/custom sites) | Skill tool (`frontend-design:frontend-design`) | Aesthetic guidance to pass as context to sub-agent |
| Before Phase 5 CSS/Styling edit (jekyll/custom sites) | Skill tool (`frontend-design:frontend-design`) | Design principles to apply during CSS implementation |
| Before Phase 5 Content edit | Skill tool (`editor`) | Prose quality guidance to apply during content implementation |
| Before Phase 5 Content edit (if style profile exists) | Skill tool (`essay-voice-matcher`) | Voice consistency notes to apply during content implementation |
| Before Quick-edit CSS/Styling (jekyll/custom sites) | Skill tool (`frontend-design:frontend-design`) | Aesthetic guidance before applying CSS change |

## Skill Invocation Pattern

The orchestrator invokes external skills BEFORE delegating to sub-agents. Key guidance is extracted and passed as inline context. This pattern ensures:
- Correct invocation ordering (guidance before edit, not after)
- Graceful degradation (skill failure never blocks the workflow)
- Site-type filtering (only invoke for applicable site types)
- Context compactness (pass 3-5 key points, not full skill output)

### When to Invoke

| Skill | Invoke when | Site-type guard | Edit-type guard |
|-------|-------------|-----------------|-----------------|
| `frontend-design:frontend-design` | Before Website Designer delegation (Phase 2), before CSS/Styling Phase 5 edit, before CSS quick-edit | `jekyll`, `custom` only | CSS/Styling edits only |
| `editor` | Before Content Phase 5 edit | `jekyll`, `custom`, `github-readme` (all content-bearing types) | Content updates only |
| `essay-voice-matcher` | Before Content Phase 5 edit (conditional) | `jekyll`, `custom`, `github-readme` | Content updates only; skip if no style profile |

### Style Profile Check (essay-voice-matcher only)

Before invoking `essay-voice-matcher`, check whether a style profile exists:

```
STYLE_PROFILE_PATH="~/.web-presence-audits/style-profile.md"
# OR per-session: {SESSION_DIR}/style-profile.md
```

If no profile is found at either location: skip essay-voice-matcher invocation. Log: "essay-voice-matcher skipped -- no style profile found."

### Graceful Degradation

If any Skill tool invocation fails (skill not found, returns empty output, or returns fewer than 50 words):
- Log: "Skill [name] invocation: [reason for skip/failure]. Continuing without guidance."
- Do NOT block the workflow.
- Do NOT retry.
- Proceed as if the skill was not invoked.

### Guidance Extraction

After invoking a skill, extract only what is needed:
- **frontend-design**: Extract 3-5 specific aesthetic principles relevant to the site's current state (e.g., typography choices, color guidance, spacing principles). Discard general philosophy.
- **editor**: Extract 3-5 prose quality guidelines relevant to the content being written (e.g., sentence length, tone, specific weaknesses to avoid).
- **essay-voice-matcher**: Extract 2-3 voice consistency notes (e.g., "use first-person active voice", "avoid passive constructions").

Pass extracted points as a compact bullet list in the delegation prompt or edit context.

## Error Handling

| Failure | Action |
|---------|--------|
| Sub-function timeout | Retry once with same instructions. If fails again: create placeholder report noting timeout. Downstream phases flag incompleteness. |
| Git clone failure | Report error per repo. Offer: retry, partial review (skip that repo), or abort. |
| Git push failure | STOP all pushes. Report which repos succeeded and which failed. Offer: retry failed push, revert successful pushes, or accept inconsistency. |
| Build validation failure | Block push for that repo. Show build errors. Offer: revert changes in working copy, attempt fix, or skip push for this repo. |
| Missing Phase 2 report | Pass warning to Phase 3/4. Coherence Manager and Suggestion Engine handle incomplete data per their instructions. |
| Skill invocation failure (skill not found or empty output) | Log warning and skip. Continue workflow without skill guidance. Never block or retry. |
| essay-voice-matcher: style profile not found | Skip invocation silently. Log: "No style profile found -- skipping voice check." |

Graceful degradation levels:
- **Full** (5/5 sub-functions): normal operation
- **Degraded** (incomplete Phase 2): proceed with warnings, scores have reduced confidence
- **Minimal** (1 sub-function succeeded): present that result directly, skip synthesis
- **Abort** (all sub-functions failed): report errors, offer retry or troubleshooting

## Timeout Configuration

| Component | Timeout | Exceeded Action |
|-----------|---------|-----------------|
| Repo cloning | 5 min per repo | Skip repo, warn user |
| Website Designer | 10 min | Placeholder report, proceed |
| Portfolio Manager | 10 min | Placeholder report, proceed |
| SEO Manager | 10 min | Placeholder report, proceed |
| Coherence Manager | 15 min | Skip coherence, present raw Phase 2 outputs |
| Suggestion Engine | 15 min | Present raw Phase 2/3 outputs to user |
| Phase 2 total (parallel) | 15 min | If exceeded: present any completed sub-function reports; mark remaining sub-functions as timed out; proceed to Gate 2 check |
| Build validation | 5 min per repo | Block push, report timeout |
| Global session ceiling | 120 min | Pause execution, save state, offer resume |

## Session Management

- **Session directory**:
  - Full review: `/tmp/web-presence-session/review-YYYY-MM/`
  - Ad-hoc: `/tmp/web-presence-session/adhoc-{area}-YYYY-MM-DD/`
  - Example: `/tmp/web-presence-session/adhoc-seo-2026-02-15/`
- **State file**: `session-state.json` -- tracks current phase, repo status, sub-function completion, and push progress
- **Lock file**: `.lock` with JSON `{"locked_at": "ISO-8601", "session_id": "review-YYYY-MM"}` -- prevents concurrent sessions. Staleness check: if lock file is older than 120 minutes (global ceiling), treat as stale and remove.
- **Resume**: On re-invocation, detect existing session and offer resume or fresh start
- **Interrupt**: If interrupted during Phase 5 pushes, push progress is tracked in session state. On resume, report which repos were pushed and which remain.
- **Cleanup**: At completion, offer to keep audit reports only (delete cloned repos) or keep everything. At start, check for sessions > 60 days old and offer cleanup.

## Invocation Mode Detection

Detect which of three modes applies:

- **Full review indicators**: "monthly review", "full review", "review all sites", "run the review"
- **Ad-hoc indicators**: "just", "only", "check [area]", "audit [area]" — area-level requests without a specific file or value
- **Quick-edit indicators**: requests that include a specific file path, CSS selector, HTML element, or exact text value — but must pass 3-component check below

**3-Component Precision Check** (required for quick-edit routing):

Classify as **precise** (quick-edit direct path) only if ALL THREE are present:
- [ ] **Target**: file path, CSS selector, or specific HTML element
- [ ] **Property**: what attribute or content to change
- [ ] **Value**: the specific new value

If any component is missing: classify as **vague** (quick-edit sub-function path) or ask the user for the missing component.

**Tie-breaking rule**: If both precise signals (specific file/selector) AND symptom signals (describes a problem) are present, route to ad-hoc mode or ask for clarification rather than direct edit.

**Routing examples:**

| User request | Mode | Why |
|---|---|---|
| "monthly review" | Full review | Trigger phrase |
| "check the SEO on my blog" | Ad-hoc | "check [area]" pattern, no specific file or value |
| "change `--primary-color` to `#2d3436` in `_sass/main.scss`" | Quick-edit precise | Target ✓ + Property ✓ + Value ✓ |
| "update `_config.yml`" | Quick-edit vague (or clarify) | Target ✓ but Property ✗ and Value ✗ — ask what to change |
| "fix the broken link on the about page" | Quick-edit vague | Symptom ("broken link") overrides page specificity — run scoped sub-function |
| "change the font to something better" | Quick-edit vague | Has property (font) but Value ✗ — ask or run Website Designer |

**Ambiguous**: If mode cannot be determined, ask: "Would you like a full monthly review, a targeted area audit, or do you have a specific change to make?"

For targeted (ad-hoc) invocations: clone only relevant repos, run only the specified
sub-function, skip coherence and synthesis unless user requests them.

## Ad-Hoc Mode Protocol

For targeted invocations, use this abbreviated lifecycle:

- **Session directory**: `/tmp/web-presence-session/adhoc-{area}-YYYY-MM-DD/`
- **State anchoring**: Use `[Ad-Hoc: {area}] {status}` format
- **Lock file**: Same timestamp-based protocol as full review
- **Quality gate**: Before presenting results, verify the sub-function report exists and contains expected section headers per its output template
- **Deployment safety**: Record pre-push SHA in `rollback-info.json` before any push. Run build validation. Show diff. Get user confirmation per repo.
- **RISKY changes**: Present as recommendations. User may override by explicitly requesting implementation after reviewing the risk description.
- **Audit persistence**: Save results to `~/.web-presence-audits/adhoc-{area}-YYYY-MM-DD.md`. These are NOT loaded as "previous audit" for monthly reviews.
- **Resume**: On re-invocation, detect existing ad-hoc session and offer resume or fresh start

## Quick-Edit Mode Protocol

Quick-edit mode skips audit phases (1–4) when the user's intent is precise enough for a direct edit. It preserves all Phase 5 safety gates.

**Depth model**: This skill implements depth 0 (fully precise: direct edit) and depth 2 (sub-function scoped: scoped audit → implement single finding). Depth 1 (file-scoped: value inference) is a named future extension.

### Quick-Edit Setup (replaces full Phase 1, for precise path)

1. Create session directory: `/tmp/web-presence-session/quickedit-YYYY-MM-DD-HHMMSS/`
2. Read site registry to identify target repo (match from file path or user context; if ambiguous, ask user which site)
3. Clone target repo only: `git clone --depth 25 --single-branch <repo-url>`
4. Create lock file: `.lock` with `{"locked_at": "<ISO-8601>", "session_id": "<quickedit-id>"}`. Respect existing lock files — if fresh lock exists, refuse with "A session is active. Complete or cancel it first."
5. Pre-flight: `git ls-remote` to verify push access. If fails, report and abort.
6. Record pre-edit SHA in `rollback-info.json` before any file modification.

### Routing Logic

**Precise path** (all 3 components confirmed):
1. Complete Quick-Edit Setup (above)
2. Show micro-context preview: read the target file, display the relevant section + 5 lines of surrounding context
3. If edit type is CSS/Styling AND site type is `jekyll` or `custom`: invoke `frontend-design:frontend-design` via Skill tool, extract 3-5 aesthetic principles. If edit type is Content: invoke `editor` via Skill tool, extract 3-5 prose guidelines; also invoke `essay-voice-matcher` if style profile exists. (Graceful degradation applies -- skip if invocation fails.)
4. Apply edit using the Phase 5 Edit Pattern Matrix for the appropriate edit type (skill guidance from step 3 incorporated)
5. Proceed to Gate 5 safety checks (NEVER skipped — see below)

**Vague path** (any component missing):
1. Determine relevant sub-function from the request area (Website Designer for CSS/layout, SEO Manager for SEO, Portfolio Manager for portfolio)
2. Run that sub-function only (reuses ad-hoc mode sub-function invocation — same delegation template)
3. Present findings to user; user selects ONE finding to implement immediately
4. Apply the selected finding using the Phase 5 Edit Pattern Matrix
5. Proceed to Gate 5 safety checks (NEVER skipped)
6. Remaining findings (not selected) are saved as recommendations but not implemented in this session

### Safety Gates (NEVER Skipped in Quick-Edit Mode)

Quick-edit skips audit phases 1–4. It does **NOT** skip Phase 5 deployment safety:

| Safety Check | Status in Quick-Edit |
|---|---|
| Pre-edit SHA recorded in `rollback-info.json` | ALWAYS active |
| Build validation per site registry `build_validation` setting | ALWAYS active |
| Diff shown to user before commit | ALWAYS active |
| User confirms push per repo | ALWAYS active |
| Specific files staged (never `git add -A`) | ALWAYS active |
| YAML validation for content/portfolio edits | ALWAYS active |
| Lock file created and respected | ALWAYS active |

### Session State for Quick-Edit

Quick-edit sessions are **non-resumable by default**. The session directory and rollback-info.json are created for safety, but if interrupted before push, the user must re-run the quick-edit. On re-invocation, detect existing quick-edit sessions and offer: resume from last known state or start fresh.

## Workflow

Phase 1: **Setup** -- Pre-flight checks, clone repos, load previous audit, initialize session. See `references/monthly-review-checklist.md` for full protocol.

Phase 2: **Parallel Analysis** -- Before launching sub-functions: for sites of type `jekyll` or `custom`, invoke the `frontend-design:frontend-design` skill via Skill tool. Extract 3-5 key aesthetic guidance points (see Skill Invocation Pattern above). Resolve `{FRONTEND_DESIGN_GUIDANCE}` in the Website Designer delegation template with this compact summary. If invocation fails or site types exclude design analysis, set `{FRONTEND_DESIGN_GUIDANCE}` to empty string and omit the guidance block from the template.

Then launch Website Designer, Portfolio Manager, and SEO Manager via Task tool **in a single message with 3 Task tool calls** (parallel execution). All three receive their instruction file path, site data, and output location via the delegation templates below. The Task tool waits for all three to return before proceeding. After all three complete, update session state as a batch (three simultaneous "running" statuses are expected; no sequential ordering dependency exists). Validate via Gate 2. Total Phase 2 wall-clock time: up to 15 minutes (bounded by slowest agent, not sum of agents).

Phase 3: **Coherence Audit** -- After Gate 2 passes, launch Coherence Manager via Task tool. It reads all repos plus Phase 2 outputs. Produces coherence-audit.md. Validate via Gate 3.

Phase 4: **Synthesis** -- After Gate 3 passes, launch Suggestion Engine via Task tool. It reads all prior outputs plus previous audit. Produces action-items.md, content-calendar.md, and audit-report.md. Validate via Gate 4.

Phase 5: **Review and Execute** -- Present audit summary and action items to user. User selects changes to implement. Group changes by repo. Deploy per the two-phase commit protocol in the checklist. Save audit report for next month.

## Phase 5 Edit Pattern Matrix

When implementing action items in Phase 5, use this matrix to identify files, validate edits, and coordinate shared file writes. For commit format, build commands, and full rollback procedures, see `references/monthly-review-checklist.md` Phase 5 Deployment Pipeline.

### Edit Type x Site Type Routing

| Edit Type | Jekyll Site | GitHub README/Profile | LaTeX CV |
|-----------|-------------|----------------------|----------|
| **CSS/Styling** | **Step 1**: Invoke `frontend-design:frontend-design` via Skill tool. Extract 3-5 aesthetic principles. (Graceful degradation applies -- skip if invocation fails.) **Step 2**: Read target file. Apply edit incorporating extracted guidance. Files: `_sass/**/*.scss` or `assets/css/*.css`. Gem theme check: if `_config.yml` has `theme:`, local `_sass/` may not exist -- create override file. Build: `bundle exec jekyll build` (required). Rollback: capture full file content before edit. | N/A -- no CSS layer. | N/A -- use LaTeX styling commands instead. |
| **Content update** | **Step 1**: Invoke `editor` via Skill tool. Extract 3-5 prose quality guidelines. **Step 2**: If style profile exists (see Skill Invocation Pattern), invoke `essay-voice-matcher`. Extract 2-3 voice notes. **Step 3**: Write content applying both sets of guidance. Files: `_posts/*.md`, `_pages/*.md`, `index.md`. YAML frontmatter: validate after every write (see YAML Validation below). Rollback: capture original frontmatter. Build: per registry `build_validation` field. (Graceful degradation applies to Steps 1-2 -- proceed without skill guidance if invocations fail.) | **Step 1**: Invoke `editor` via Skill tool. Extract 3-5 prose quality guidelines. **Step 2**: Write content applying guidance. File: `README.md`. Direct markdown edit. No build step required. Rollback: capture full file content. (Graceful degradation applies.) | Files: `*.tex` (target section). No YAML. Build: `latexmk` or `pdflatex`. Rollback: capture original content. |
| **SEO fixes** | Files: `_config.yml` (plugins, url, description), `_includes/head.html` (meta tags). Batch all `_config.yml` edits (see Shared File Coordination below). Build: required. Rollback: capture original content of each file. | File: `README.md` (keywords in bio, section headings for discoverability). No build. Rollback: capture full file content. | N/A -- CV does not have web SEO. |
| **Portfolio entries** | Files: `_projects/*.md` (new post) or `_data/projects.yml` (append). YAML frontmatter: validate after create/edit. May require adding collection to `_config.yml`. Batch `_config.yml` changes. Build: required. Rollback: capture pre-edit state of each modified file. | File: `README.md` (add entry to pinned projects section). Direct markdown append. No build. | Files: `*.tex` (append entry to projects or publications section). No YAML. Build: `latexmk`. |

### YAML Validation

After any write that includes YAML frontmatter (content updates, portfolio entries):

1. Extract frontmatter (content between `---` delimiters)
2. Validate: `python3 -c "import yaml; yaml.safe_load('''<frontmatter>''')"`
3. If validation fails: REVERT the edit, report the specific YAML error with suggested fix
4. Common fix: quote values containing `: { } [ ] # & * ! | > ' " % @` characters

Example of fragile frontmatter (title with colon):
```yaml
# BREAKS: title: My Project: A Deep Dive
# SAFE:   title: "My Project: A Deep Dive"
```

### Shared File Coordination (_config.yml)

Multiple edit types (SEO, content, portfolio) may all modify `_config.yml`. To prevent overwrites:

1. Before beginning Phase 5 execution, identify ALL pending action items that target `_config.yml`
2. Batch all `_config.yml` changes into a single read-modify-write cycle (read once, apply all changes, write once)
3. If edit patterns run sequentially and each reads the current state, order them SEO -> content -> portfolio so each sees the previous change

### Composite Edits

If a single action item spans multiple edit types (e.g., a portfolio entry needs a new `_projects/` post + a new CSS card style + an `_config.yml` collection entry):

1. Execute edit types in order: `_config.yml` changes first, then content files, then CSS last
2. Group all changes in a single `git commit` per repo
3. Run build validation once after all edits are applied, not per edit type

## Task Delegation Templates

When delegating to sub-functions via Task tool, use these templates. The orchestrator MUST resolve all `{VARIABLES}` before invoking the Task tool. Sub-agents cannot see files the orchestrator has read -- all context must be in the prompt.

### Resolving Variables

Before Phase 2, resolve these once and reuse across all delegations:

- `{SKILL_DIR}`: The absolute path to this skill's directory. Resolve dynamically by checking which path is readable: first try the synced location `~/.claude/skills/web-presence-manager/`, then fall back to the repo location.
- `{SESSION_DIR}`: `/tmp/web-presence-session/review-YYYY-MM/` (or `adhoc-{area}-YYYY-MM-DD/` for ad-hoc)
- `{SITES}`: For each site, list: name, type, repo path in session directory, applicable status. Inline this data directly -- do not make the sub-agent parse the registry.

### Template: Website Designer (Phase 2)

````
Read your instructions from {SKILL_DIR}/references/website-designer-instructions.md

Analyze these sites:
{For each site where website-designer is in applicable_sub_functions:}
- {site.name} (type: {site.type}) at {SESSION_DIR}/repos/{repo-name}/

Skip sites not listed above.

{AESTHETIC_GUIDANCE_BLOCK}

Write your output to: {SESSION_DIR}/outputs/design-review.md

You have access to: Read, Bash, Glob tools. Do NOT modify site files.
````

Where `{AESTHETIC_GUIDANCE_BLOCK}` resolves to either:
- **If frontend-design guidance was obtained**:
  ```
  Aesthetic guidance from frontend-design (apply these principles to your Quick Wins and Recommendations):
  {FRONTEND_DESIGN_GUIDANCE}

  Add an "Aesthetic Distinctiveness" subsection to the Current Assessment in your output (see output template in your instructions).
  ```
- **If no guidance** (skill failed, wrong site type, or no applicable sites): omit the block entirely.

### Template: Portfolio Manager (Phase 2)

````
Read your instructions from {SKILL_DIR}/references/portfolio-manager-instructions.md

Analyze these sites:
{For each site where portfolio-manager is in applicable_sub_functions:}
- {site.name} (type: {site.type}) at {SESSION_DIR}/repos/{repo-name}/

Skip sites not listed above.

Write your output to: {SESSION_DIR}/outputs/portfolio-review.md

You have access to: Read, Bash, Glob tools. Do NOT modify site files.
````

### Template: SEO Manager (Phase 2)

````
Read your instructions from {SKILL_DIR}/references/seo-manager-instructions.md

Analyze these sites:
{For each site where seo-manager is in applicable_sub_functions:}
- {site.name} (type: {site.type}, url: {site.url}) at {SESSION_DIR}/repos/{repo-name}/

Skip sites not listed above.

Write your output to: {SESSION_DIR}/outputs/seo-audit.md

You have access to: Read, Bash, Glob, Grep, WebSearch (if available) tools. Do NOT modify site files.
````

### Template: Coherence Manager (Phase 3)

````
Read your instructions from {SKILL_DIR}/references/coherence-manager-instructions.md

Primary site: {primary_site.name} at {SESSION_DIR}/repos/{primary_repo-name}/

All sites to compare:
{For each site where coherence-manager is in applicable_sub_functions:}
- {site.name} (type: {site.type}) at {SESSION_DIR}/repos/{repo-name}/

Phase 2 outputs are located at:
- {SESSION_DIR}/outputs/design-review.md
- {SESSION_DIR}/outputs/portfolio-review.md
- {SESSION_DIR}/outputs/seo-audit.md

Write your outputs to:
- {SESSION_DIR}/outputs/brand-reference.md
- {SESSION_DIR}/outputs/coherence-audit.md

You have access to: Read, Glob, Grep, Write tools. Use Write ONLY for your two output files. Do NOT use Bash for git operations. Do NOT modify site files.
````

### Template: Suggestion Engine (Phase 4)

````
Read your instructions from {SKILL_DIR}/references/suggestion-engine-instructions.md

Read all prior outputs from:
- {SESSION_DIR}/outputs/design-review.md
- {SESSION_DIR}/outputs/portfolio-review.md
- {SESSION_DIR}/outputs/seo-audit.md
- {SESSION_DIR}/outputs/coherence-audit.md

Previous audit (for trend comparison): {previous_audit_path or "None (first run)"}

Write your three deliverables to:
- {SESSION_DIR}/outputs/action-items.md
- {SESSION_DIR}/outputs/content-calendar.md
- {SESSION_DIR}/outputs/audit-report.md

You have access to: Read, Write tools. Use Write ONLY for the three deliverable files. Do NOT use Bash, WebSearch, Glob, or Grep.
````

## Reference Files

All paths relative to this skill's directory. The orchestrator resolves absolute paths before delegating.

- `references/site-registry.md` -- managed sites configuration (single source of truth for repos)
- `references/website-designer-instructions.md` -- Sub-function 1: visual design and accessibility
- `references/portfolio-manager-instructions.md` -- Sub-function 2: portfolio currency and completeness
- `references/seo-manager-instructions.md` -- Sub-function 3: technical and content SEO
- `references/coherence-manager-instructions.md` -- Sub-function 4: narrative and visual coherence
- `references/suggestion-engine-instructions.md` -- Sub-function 5: synthesis and prioritization
- `references/monthly-review-checklist.md` -- full review protocol, quality gates, deployment pipeline, rollback

See Sub-Function Responsibility Boundaries table in `references/monthly-review-checklist.md` for a RACI-style mapping of which sub-function handles which concern.
