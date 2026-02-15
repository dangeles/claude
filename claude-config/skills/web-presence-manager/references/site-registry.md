# Site Registry

Configuration for all managed web properties. This file is the single source
of truth for which repositories are managed, what type they are, and which
sub-functions apply to each. To add a new site, copy an existing entry and
fill in all required fields.

## Schema

Each site entry MUST include:

- **name**: Display name for reports
- **repo**: GitHub repository path (owner/repo)
- **type**: `jekyll` | `latex` | `github-readme` | `custom`
- **branch**: Branch to operate on (default: `main`)
- **content_areas**: List from `[bio, portfolio, blog, cv, services, team]`
- **applicable_sub_functions**: Which sub-functions analyze this site
- **build_command**: Command to validate before push (or `none`)
- **build_validation**: `required` (block push on failure) | `optional` (warn on failure, allow push) | `none` (no build step)

Optional fields:

- **url**: Published URL (for SEO analysis of live site)
- **primary**: `true` if this is the authoritative source for coherence resolution
- **notes**: Free-text notes for sub-functions

---

## Sites

### Site: Personal Website

- **name**: Personal Website
- **repo**: dangeles/dangeles.github.io
- **type**: jekyll
- **branch**: main
- **url**: https://dangeles.github.io
- **primary**: true
- **content_areas**: [bio, portfolio, blog]
- **applicable_sub_functions**: [website-designer, portfolio-manager, seo-manager, coherence-manager]
- **build_command**: `bundle exec jekyll build`
- **build_validation**: required
- **notes**: Primary site. Jekyll theme. Source of truth for coherence resolution. All sub-functions except suggestion-engine apply directly.

### Site: GitHub Profile

- **name**: GitHub Profile README
- **repo**: dangeles/dangeles
- **type**: github-readme
- **branch**: main
- **url**: https://github.com/dangeles
- **content_areas**: [bio, portfolio]
- **applicable_sub_functions**: [portfolio-manager, coherence-manager]
- **build_command**: none
- **build_validation**: none
- **notes**: Profile README. Markdown only. No build step. Website Designer skips this site. SEO Manager skips this site.

### Site: CV

- **name**: Academic CV
- **repo**: dangeles/cv
- **type**: latex
- **branch**: main
- **content_areas**: [cv, bio, portfolio]
- **applicable_sub_functions**: [portfolio-manager, coherence-manager]
- **build_command**: `pdflatex main.tex`
- **build_validation**: optional
- **notes**: LaTeX CV. Portfolio Manager extracts structured data (publications, experience, education). Parsing may require fallback to raw text search if macros are non-standard. Website Designer and SEO Manager skip this site.

---

## Adding a New Site

To add a site, copy one of the entries above and fill in all required fields.
The orchestrator reads this registry at Phase 1 startup and routes sub-functions
accordingly. No changes to SKILL.md or instruction files are needed.

<!-- Uncomment and fill in when company website is ready:
### Site: Company Website

- **name**: Company Website
- **repo**: company/website
- **type**: custom
- **branch**: main
- **url**: https://company.com
- **content_areas**: [services, team, bio]
- **applicable_sub_functions**: [website-designer, seo-manager, coherence-manager]
- **build_command**: `npm run build`
- **build_validation**: required
- **notes**: Custom stack. Verify build_command is correct before first review. Check for framework-specific SSR/SSG patterns.
-->
