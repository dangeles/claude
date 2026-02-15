# Example: Ad-Hoc Targeted Invocation

## User Invocation

> Check the SEO on my blog

## Mode Detection

```
[Mode Detection] Detected targeted invocation.

  Trigger phrase: "Check the SEO on my blog"
  Matched pattern: "check [area]" -> ad-hoc targeted mode
  Target sub-function: SEO Manager
  Target site: dangeles/dangeles.github.io (Jekyll blog)

  Running SEO Manager only. Skipping: Website Designer, Portfolio Manager,
  Coherence Manager, Suggestion Engine.

  Proceed with targeted SEO audit of dangeles.github.io? (y/n)
```

> y

## Setup (abbreviated)

```
[Setup] Targeted mode -- cloning only relevant repository

  Cloning dangeles/dangeles.github.io (shallow)...
  [OK] dangeles.github.io -> /tmp/web-presence-session/adhoc-seo-2026-02/repos/

  Pre-flight:
    [PASS] git push access verified
    [WARN] Jekyll not found -- build validation skipped
```

## Targeted Execution

```
[Analysis] Running SEO Manager on dangeles.github.io

  Checking Jekyll SEO plugin stack...
  Analyzing meta tags, heading structure, structured data...
  Checking sitemap, robots.txt, canonical URLs...

  [DONE] seo-audit.md

  SEO Score: 62/100
```

## Results

```
=== SEO Audit: dangeles.github.io ===

Score: 62/100

Plugin Stack:
  | Plugin          | Status    | Action                          |
  |-----------------|-----------|----------------------------------|
  | jekyll-seo-tag  | installed | none                            |
  | jekyll-sitemap  | missing   | Add to Gemfile and _config.yml  |
  | jekyll-feed     | installed | none                            |

Critical Issues:
  1. Missing Person schema JSON-LD (SAFE edit)
     - No structured data for author identity
     - Fix: Add JSON-LD to _includes/head.html

  2. Missing jekyll-sitemap plugin (SAFE edit)
     - No sitemap.xml being generated
     - Fix: Add to Gemfile and plugins list in _config.yml

  3. Canonical URL mismatch (RISKY -- recommendation only)
     - _config.yml url: "" (empty)
     - Should be: "https://dangeles.github.io"
     - This affects ALL canonical URLs site-wide

Improvements:
  4. 4 images missing alt text across blog posts (SAFE edit)
  5. Blog post descriptions missing in front matter (MODERATE edit)

Would you like to implement any of these fixes?
```

## User Selects Changes

> Fix items 1, 2, and 4. For item 3 I want to review the config change first
> before committing.

```
[Execute] Implementing approved changes

  Adding Person schema JSON-LD to _includes/head.html... done
  Adding jekyll-sitemap to Gemfile... done
  Adding jekyll-sitemap to _config.yml plugins list... done
  Adding alt text to 4 blog images... done

  For item 3 (canonical URL), here is the proposed change:

    File: _config.yml
    Line 1:  url: ""
    Change:  url: "https://dangeles.github.io"

    This is categorized as RISKY because changing the url field affects
    all canonical URLs, sitemap URLs, and feed URLs site-wide. It is the
    correct value, but verify this matches your actual serving domain.

    Apply this change? (y/n)
```

> y, that looks right

```
  Updating _config.yml url field... done

  Staging files:
    _includes/head.html
    Gemfile
    _config.yml
    _posts/2025-11-15-ml-paper.md
    _posts/2025-12-01-year-review.md

  Commit:
    seo(blog): add sitemap plugin, Person schema, alt text, and canonical URL

    - Added jekyll-sitemap plugin for automatic sitemap.xml generation
    - Added schema.org/Person JSON-LD with sameAs links
    - Added descriptive alt text to 4 blog images
    - Set canonical URL to https://dangeles.github.io

    Reviewed-by: web-presence-manager ad-hoc audit (2026-02)

  Build validation: SKIPPED (Jekyll not available)
  Diff: 5 files changed, +38 -3

  Push to dangeles.github.io origin/main? (y/n)
```

> y

```
  [PUSHED] dangeles.github.io -> origin/main (sha: d4e5f6a)

  Done. SEO score estimated improvement: 62 -> ~78 after changes are indexed.

  Remaining items for future:
    - Blog post descriptions in front matter (5 posts missing)
    - Consider adding google_site_verification to _config.yml

  Tip: Run a full monthly review for a complete audit across all sites.
```
