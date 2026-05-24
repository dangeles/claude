# Plotting Advisor — Design Spec

**Status**: Draft for review
**Author**: dangeles (via brainstorming session, 2026-05-24)
**Repo**: `claude-config/` (skills live in `claude-config/skills/plotting-advisor/`)

## Objective

Add a new skill, `plotting-advisor`, that acts as a **rules engine** for Python plotting. When a user or another agent is about to create a plot, the advisor returns a structured checklist of choices (chart type, palette, axes, annotation, accessibility) plus a paste-ready YAML decision card. Library mechanics are delegated to the existing `scientific-skills:matplotlib`, `scientific-skills:seaborn`, and `scientific-skills:plotly` skills. A passive lint mode also exists for checking already-rendered figures.

## Scope

**In scope:**

- Style guidance (palette, axes, annotation, accessibility) drawn from Tufte, Cleveland, Wong, Wilke, Healy, Bergstrom & West
- Chart-type selection by data shape × intent
- Static-first coverage (matplotlib, seaborn) with a short interactive-adaptation reference for plotly/bokeh/altair
- Lint flow for existing figures (PNG, matplotlib `Figure` via spec extraction, or textual description)
- A separate `scientific-conventions.md` reference covering volcano, UMAP/t-SNE/PCA, Manhattan, forest, hierarchical-clustered heatmap, Kaplan-Meier, ROC, Bland-Altman, raincloud

**Out of scope:**

- Writing plot code itself (advisor pattern — delegates syntax)
- Non-Python plotting (R/ggplot, D3, etc.)
- Domain-specific scientific conventions inside `chart-selection.md` (those live in `scientific-conventions.md`)
- Replacing or duplicating `scientific-skills:matplotlib` / `seaborn` / `plotly` content

## Success criteria

- Chart-type recommendation justified by data shape + intent with an inline citation
- Palette selected is colorblind-safe OR exception flagged with reason
- Axis treatment (scale, zero, breaks) explicitly addressed
- Decision card produced as a paste-ready YAML block including a `delegate_to` field pointing to the right software skill
- Accessibility floor (palette safety, font ≥8pt print / 14pt slide, redundant encoding) met or violations flagged
- Lint flow exits 0 by default; `--strict` flag exits 2 on critical violations for CI use

## Architecture

### Skill identity

```yaml
---
name: plotting-advisor
version: 1.0
last_updated: 2026-05-24
description: "Use BEFORE writing any Python plotting code (matplotlib, seaborn, plotly) — recommends chart type, color palette, axis treatment, and accessibility choices using Tufte, Cleveland, Wong, and Wilke principles. Returns a structured checklist + a paste-ready decision card. Also lints existing figures on request. Triggers on 'plot', 'chart', 'visualize', 'graph', 'figure', 'heatmap', 'histogram', 'scatter', 'bar plot', or any mention of matplotlib/seaborn/plotly."
prerequisites:
  - Data shape known (variable types, n)
  - Plotting intent stated (compare, distribution, trend, composition, relation)
constraints:
  - Does NOT write plotting code — delegates syntax to scientific-skills:matplotlib / seaborn / plotly
  - Lint is passive — caller provides Figure object or PNG path; advisor never imports matplotlib to render
success_criteria:
  - Chart-type recommendation justified by data shape + intent
  - Palette selected is colorblind-safe OR exception flagged with reason
  - Axis treatment (scale, zero, breaks) explicitly addressed
  - Decision card produced as paste-ready block
  - Accessibility minimums met (font size, redundant encoding when color encodes a meaningful variable)
extended_thinking_budget: 2048
---
```

### File layout

```
claude-config/skills/plotting-advisor/
├── SKILL.md
├── references/
│   ├── chart-selection.md          (~800 words; general only — no domain entries)
│   ├── palettes.md                 (~600 words; Okabe-Ito, viridis family, ColorBrewer, Tol)
│   ├── anti-patterns.md            (~900 words; 12-15 entries with citations)
│   ├── accessibility.md            (~500 words; WCAG, colorblind sim, font, B/W fallback)
│   ├── interactive-adaptation.md   (~300 words; how rules shift for plotly/bokeh/altair)
│   ├── scientific-conventions.md   (~600 words; volcano, UMAP, Manhattan, forest, etc.)
│   └── SOURCES.md                  (single bibliography; inline cites use short tags)
└── scripts/
    ├── style_lint.py               (CLI entry point for lint flow)
    ├── figure_spec.py              (~50 lines; matplotlib Figure → JSON spec — usable as both import and CLI)
    └── palettes.py                 (canonical hex lists — used by lint and refs)
```

### Two flows, shared rule base

**Advisor flow** (default trigger):

1. Intake: data shape, intent, audience, constraints. If intent missing, ask one question.
2. Chart selection — consult `references/chart-selection.md`, walk the decision tree.
3. Palette — consult `references/palettes.md`.
4. Axis & scale, annotation, accessibility floor — inline rules in SKILL.md.
5. Anti-pattern veto check — inline (top 5: 3D, dual y-axis, rainbow on continuous, truncated baseline, pie >5).
6. Produce checklist + decision card (formats below).

**Lint flow** (when caller hands over a figure):

1. Input mode (mutually exclusive):
   - `--image figure.png` (Pillow + optional Tesseract OCR)
   - `--figure-spec figure.json` (caller's `extract_spec(fig)` output)
   - `--describe "log-y bar chart of 7 categories using rainbow palette, no axis labels"`
2. Inspection → rule check → grouped violations.
3. Exit 0 by default; `--strict` exits 2 on critical violations.

### What lives inline in SKILL.md vs references/

Split is by **size + universality**:

| Lives | Why |
|---|---|
| **Inline** — axis & scale, annotation, accessibility floor, top-5 anti-pattern vetoes | Short (4-6 bullets each), applies to every plot. Moving to references would force an extra file-read on every invocation. |
| **References** — chart-selection tree, palette catalog, long anti-pattern list, accessibility deep-dive, interactive adaptation, scientific conventions, bibliography | Large or branchy. Caller needs only one branch / one subsection at a time. Progressive disclosure earns its keep. |

The `accessibility` topic appears in both because of a floor-vs-deep split: 5 universal minimums are always checked inline; WCAG specifics, colorblind simulation, journal-specific overrides live in the reference.

## Output formats

### Advisor — Checklist + Decision card

```markdown
## Plotting Advisor: [chart type recommendation in one phrase]

### Intent & data
- Intent: ...
- Data: ...
- Audience: ...

### Recommended chart: [type]
- Rule: [citation tag, e.g., Cleveland 1985 — position encodes more accurately than angle]
- Alternatives considered: ...

### Palette: [name]
- Hex: ...
- Why: ...
- Reference: references/palettes.md#anchor

### Axes & scale
- ...

### Annotation
- ...

### Accessibility floor
- Colorblind-safe / redundant encoding / font / contrast — checked

### Anti-pattern check
- No 3D / no dual axes / no rainbow / no truncated baseline / no pie misuse — pass/flag

### Decision card

​```yaml
chart: <type>
library: <matplotlib | seaborn | plotly>
palette:
  type: <categorical | sequential | diverging>
  name: <okabe_ito | viridis | rdbu | ...>
  hex: [...]
axes:
  x: {type: ..., order: ...}
  y: {type: ..., include_zero: ..., label: "..."}
encoding:
  position: <variable>
  color: <variable>
annotation:
  direct_labels: <bool>
  sample_size_per_group: <bool>
  title: "..."
  subtitle: "..."
accessibility:
  colorblind_safe: <bool>
  min_font_pt: <int>
output:
  dpi: <int>
  format: <pdf | svg | png>
delegate_to: scientific-skills:<matplotlib | seaborn | plotly>
​```
```

### Lint — grouped violations

```markdown
## Plotting Advisor: figure check — N issues found

### Critical (n)
- <issue>. <impact>.
  - Fix: <one-line fix>
  - Rule: references/<file>.md#<anchor>

### Major (n)
- ...

### Minor (n)
- ...
```

Plus JSON to stdout for machine consumption.

## Lint script (`scripts/style_lint.py`)

### Modes

1. `--image PATH` — Pillow + optional OCR. Quantizes anti-aliased PNGs to ~16 representative colors before checking palette.
2. `--figure-spec PATH` — reads JSON produced by `figure_spec.extract_spec(fig)`. Cleanest input mode. The helper is usable two ways: (a) as an import via `sys.path.insert(0, "<skill>/scripts"); from figure_spec import extract_spec`, or (b) as a CLI: `python3 scripts/figure_spec.py --pickle fig.pkl --out figure.json`. Skill docs include a copy-paste snippet.
3. `--describe TEXT` — parses textual description; for ad-hoc checks without a saved figure.

### Checks (matching the rule base)

| Category | Check | Severity |
|---|---|---|
| Axis | Both axes labeled | Critical |
| Axis | Units in label (parenthesized) | Major |
| Axis | Linear axis includes zero OR annotated truncation | Critical |
| Axis | Log only when range ≥2 decades | Major |
| Scale | ≤~7 major ticks per axis | Minor |
| Color | Palette colorblind-safe (whitelist OR ΔE in CIELAB) | Major |
| Color | ≤8 distinct categorical colors | Major |
| Color | No `jet`/`gist_rainbow`/`hsv` on continuous | Critical |
| Encoding | Redundant encoding when color encodes categorical | Minor |
| Font | ≥8pt print / ≥14pt slide (caller declares context) | Major |
| Anti-pattern | No 3D projection | Critical |
| Anti-pattern | No dual y-axis | Critical |
| Anti-pattern | Pie with >5 slices | Major |
| Annotation | Title present | Major |
| Annotation | Legend when ≤5 series — suggest direct labels | Minor |

### What lint does NOT do

- Does not import matplotlib to render
- Does not modify the figure
- Does not download external data
- Does not call out to an LLM — pure rule matcher

### Failure modes

- Corrupt image → `severity: error` record, exit 0 (lint is advisory)
- OCR unavailable → skip label-presence in image mode, note the gap
- Malformed spec JSON → one-line diagnostic, exit 0
- Anti-aliased palette ambiguity → quantize to 16, check whitelist + ΔE

### Exit behavior

- Default: exit 0 always (advisory)
- `--strict`: exit 2 if any critical violation present (for CI use)

## Reference doc contents (summary)

### `chart-selection.md`

Decision tree organized by **intent × data shape**. Intents: compare, distribution, trend, composition, relation, ranking, geographic, network. Each branch cites the rule that selected it (Cleveland, Tufte, Wilke, Healy). **Excludes** domain-specific charts.

### `palettes.md`

Palette catalog: Okabe-Ito (categorical ≤8), viridis/cividis (sequential), ColorBrewer RdBu/PuOr (diverging), Tol bright (alt categorical ≤7). Hex codes, when-to-use, when-not, colorblind-safety guarantee. Citations: Okabe & Ito 2008, Wong 2011, Tol 2018.

### `anti-patterns.md`

12-15 entries with practice → why-wrong (perceptual or data-integrity rule) → fix → citation. Includes 3D, dual axes, rainbow, truncated baseline, pie misuse, dynamite plots, aspect ratio gaming, legend overuse, gridline competition, chartjunk, alphabet-order, count-vs-rate, overplotting, uninformative heatmap order.

### `accessibility.md`

Colorblind simulation workflow, WCAG 2.2 contrast (4.5:1 text, 3:1 graphical), redundant encoding, font sizing by medium, journal-specific overrides (Nature, Cell, eLife), screen-reader caption guidance, B/W fallback test.

### `interactive-adaptation.md`

How rules shift for plotly/bokeh/altair: hover relaxes label density, zoom invalidates "always include zero" (annotate default view), color budget can expand, animation respects reduced-motion, static export fallback required.

### `scientific-conventions.md`

Domain-specific charts: volcano, UMAP/t-SNE/PCA, Manhattan, forest, hierarchical-clustered heatmap, Kaplan-Meier, ROC, Bland-Altman, raincloud. Each: the convention, canonical example, common failure modes.

### `SOURCES.md`

Single bibliography. Inline cites use short tags (`[Cleveland 1985]`); SOURCES.md has full entries.

## Triggering

The skill is **auto-triggered** via its description field — front-loaded verbs (`plot`, `chart`, `visualize`, `graph`, `figure`, `heatmap`, `histogram`, `scatter`, `bar plot`) and library mentions (`matplotlib`, `seaborn`, `plotly`). The "Use BEFORE writing any Python plotting code" phrase establishes the discipline: consult before, not after.

When a figure already exists and the caller asks "check this" / "review this figure" / "make this better", the skill enters lint mode.

## Integration with existing skills

- **`scientific-skills:matplotlib` / `seaborn` / `plotly`**: explicit `delegate_to` field in the decision card names the right software skill. The advisor never duplicates their syntax content.
- **`bioinformatician`**: its current reference to a missing `visualization_best_practices.md` is satisfied by pointing at `plotting-advisor/references/scientific-conventions.md` (follow-up; not part of this skill's PR).
- **`notebook-writer`**: notebooks that include plotting cells should invoke `plotting-advisor` for each figure cell.
- **`editor`** / `consistency-auditor`: model precedent for the rules-engine-as-skill pattern.

## Open questions / deferred

- Whether to ship a small Jupyter cell magic (`%%plotting-advisor`) that wraps the lint call. Deferred — adds dependency surface; can be a follow-up.
- Whether `bioinformatician/SKILL.md` should be updated in the same PR to point at the new scientific-conventions reference. Recommend: separate PR after this one merges.
- Whether to add `pytest` smoke tests for the lint script. Recommend: yes, lightweight tests for each input mode + each severity-bearing rule.

## Risks

- **Reference doc rot** — citations and palette tables go stale. Mitigation: `SOURCES.md` as single source of truth; date-stamped at top.
- **Over-triggering** — description verbs may fire on debug/exploratory plot requests. Mitigation: explicit "Do NOT use this skill when… exploratory glance during debugging" in SKILL.md.
- **Lint false positives** — image-mode color extraction is approximate. Mitigation: quantize before checking; never gate behavior on image-mode alone for borderline cases; figure-spec mode is the recommended path.
- **Scope creep into "author" territory** — lint script could grow into a fix-it. Mitigation: spec explicitly forbids modifying figures; reinforced in CONTRIBUTING note in SKILL.md.

## Build sequence

1. Create skill directory + SKILL.md with frontmatter
2. Write inline rule sections (axis, annotation, accessibility floor, top-5 anti-patterns)
3. Write reference docs in order: SOURCES → chart-selection → palettes → anti-patterns → accessibility → interactive-adaptation → scientific-conventions
4. Write `scripts/palettes.py`, `scripts/figure_spec.py`, `scripts/style_lint.py`
5. Add smoke tests for the lint script
6. Run `./sync-config.py push --dry-run`, then `push`
7. Smoke-test invocation with a fake "plot 4 conditions, n=30" request
8. Commit + planning entry

## Approval

This spec was developed interactively with the user across 6 design sections. Awaiting written-spec review before transition to `writing-plans`.
