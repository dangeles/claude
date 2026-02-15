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

## Error Handling

| Failure | Action |
|---------|--------|
| Sub-function timeout | Retry once with same instructions. If fails again: create placeholder report noting timeout. Downstream phases flag incompleteness. |
| Git clone failure | Report error per repo. Offer: retry, partial review (skip that repo), or abort. |
| Git push failure | STOP all pushes. Report which repos succeeded and which failed. Offer: retry failed push, revert successful pushes, or accept inconsistency. |
| Build validation failure | Block push for that repo. Show build errors. Offer: revert changes in working copy, attempt fix, or skip push for this repo. |
| Missing Phase 2 report | Pass warning to Phase 3/4. Coherence Manager and Suggestion Engine handle incomplete data per their instructions. |

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
| Phase 2 total (sequential) | 30 min | If exceeded: present completed reports, skip remaining sub-functions |
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

Detect whether this is a full review or targeted ad-hoc invocation:

- **Full review indicators**: "monthly review", "full review", "review all sites", "run the review"
- **Targeted indicators**: "just", "only", "check [area]", "update [thing]", "fix [specific item]"
- **Ambiguous**: Ask user to clarify: "Would you like a full monthly review across all sites, or a targeted [area] check?"

For targeted invocations: clone only relevant repos, run only the specified
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

## Workflow

Phase 1: **Setup** -- Pre-flight checks, clone repos, load previous audit, initialize session. See `references/monthly-review-checklist.md` for full protocol.

Phase 2: **Sequential Analysis** -- Launch Website Designer, Portfolio Manager, and SEO Manager via Task tool sequentially (one at a time; Task tool blocks until each completes). Each sub-function receives its instruction file path, site data, and output location via the delegation template below. Wait for all three. Validate via Gate 2. Total Phase 2 time: up to 30 minutes.

Phase 3: **Coherence Audit** -- After Gate 2 passes, launch Coherence Manager via Task tool. It reads all repos plus Phase 2 outputs. Produces coherence-audit.md. Validate via Gate 3.

Phase 4: **Synthesis** -- After Gate 3 passes, launch Suggestion Engine via Task tool. It reads all prior outputs plus previous audit. Produces action-items.md, content-calendar.md, and audit-report.md. Validate via Gate 4.

Phase 5: **Review and Execute** -- Present audit summary and action items to user. User selects changes to implement. Group changes by repo. Deploy per the two-phase commit protocol in the checklist. Save audit report for next month.

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

Write your output to: {SESSION_DIR}/outputs/design-review.md

You have access to: Read, Bash, Glob tools. Do NOT modify site files.
````

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
