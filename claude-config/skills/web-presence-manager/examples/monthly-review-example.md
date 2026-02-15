# Example: Full Monthly Review

## User Invocation

> Run my monthly web presence review

## Phase 1: Setup and Validation

```
[Phase 1/5 - Setup] Starting monthly web presence review for February 2026

Pre-flight checks:
  [PASS] git 2.42.0 installed
  [PASS] git push access to dangeles/dangeles.github.io
  [PASS] git push access to dangeles/dangeles
  [PASS] git push access to dangeles/cv
  [WARN] Jekyll not found -- build validation will be skipped for Jekyll sites
  [PASS] pdflatex available
  [PASS] Session directory created: /tmp/web-presence-session/review-2026-02/

Cloning repositories (shallow):
  [OK] dangeles.github.io -> repos/dangeles.github.io/
  [OK] dangeles -> repos/dangeles/
  [OK] cv -> repos/cv/

Previous audit: Loaded from ~/.web-presence-audits/audit-2026-01.md
Session state initialized. Proceeding to Phase 2.
```

## Phase 2: Parallel Analysis

```
[Phase 2/5 - Analysis] Launching 3 sub-functions (wave-based, 15s stagger)

  T=0s:  Website Designer  -> analyzing dangeles.github.io...
  T=15s: Portfolio Manager -> analyzing all 3 sites...
  T=30s: SEO Manager       -> analyzing dangeles.github.io...

Sub-function results:
  [DONE] Website Designer:  design-review.md (Design Score: 7/10)
  [DONE] Portfolio Manager: portfolio-review.md (Portfolio Score: 6/10)
  [DONE] SEO Manager:       seo-audit.md (SEO Score: 62/100)

Gate 2 validation: 3/3 reports present, all > 200 words. PASS.
```

## Phase 3: Coherence Audit

```
[Phase 3/5 - Coherence] Launching Coherence Manager

  Extracting brand reference from primary site (dangeles.github.io)...
  Comparing narrative consistency across 3 sites...
  Comparing visual identity across applicable sites...

  [DONE] Coherence Manager: coherence-audit.md
    Narrative Coherence: 7/10
    Visual Coherence: 8/10

Gate 3 validation: Report complete with both sections. PASS.
```

## Phase 4: Synthesis

```
[Phase 4/5 - Synthesis] Launching Suggestion Engine

  Reading all Phase 2 and Phase 3 outputs...
  Loading previous audit (January 2026) for trend comparison...
  Generating deliverables...

  [DONE] action-items.md (8 items: 3 Must Do, 3 Should Do, 2 Backlog)
  [DONE] content-calendar.md (2 blog ideas, 1 portfolio update, 1 profile update)
  [DONE] audit-report.md (executive summary + scores + priority matrix)

Gate 4 validation: All 3 deliverables present with content. PASS.
```

## Phase 5: Review and Execute

```
[Phase 5/5 - Review] Presenting results

=== Web Presence Audit: February 2026 ===

Executive Summary:
  Web presence is in good shape with a few areas needing attention. Portfolio
  is missing two recent projects, and SEO score dropped due to missing
  structured data. Design improved since January after CSS cleanup.

Scores:
  | Area                | Score  | Trend  |
  |---------------------|--------|--------|
  | Design              | 7/10   | up     |
  | Portfolio            | 6/10   | stable |
  | SEO                 | 62/100 | down   |
  | Coherence (Narrative)| 7/10   | stable |
  | Coherence (Visual)  | 8/10   | up     |

Must Do items:
  1. Add Person schema JSON-LD to personal site (SEO, SAFE edit)
  2. Add missing alt text to 4 blog images (SEO/Accessibility, SAFE edit)
  3. Update CV with recent publication from December 2025
```

### User Interaction

> Implement items 1 and 2. I will handle item 3 manually since I need to
> update the BibTeX file too.

```
[Phase 5/5 - Execute] Implementing approved changes

Repository: dangeles.github.io

  Staging: _includes/head.html (Person schema JSON-LD added)
  Staging: _posts/2025-11-15-ml-paper.md (alt text added to 2 images)
  Staging: _posts/2025-12-01-year-review.md (alt text added to 2 images)

  Commit message:
    seo(blog): add Person schema JSON-LD and missing image alt text

    - Added schema.org/Person JSON-LD with sameAs links to GitHub and CV
    - Added descriptive alt text to 4 blog post images

    Reviewed-by: web-presence-manager monthly audit (2026-02)

  Build validation: SKIPPED (Jekyll not available)
  Diff preview: [showing 4 files changed, +42 -2]

  Push dangeles.github.io to origin/main? (y/n)
```

> y

```
  [PUSHED] dangeles.github.io -> origin/main (sha: a1b2c3d)

No changes for other repositories.

=== Review Complete ===

  Pushed: dangeles.github.io (1 commit)
  Skipped: dangeles (no changes), cv (user will handle manually)

  Audit report saved to: ~/.web-presence-audits/audit-2026-02.md

  Next review recommended: March 2026
  Key items to track:
    - SEO score (currently 62, target 80+ after schema addition)
    - Portfolio currency (2 projects still need adding)
    - CV update (user handling manually)

Session cleanup: keeping audit reports, deleting cloned repos.
```
