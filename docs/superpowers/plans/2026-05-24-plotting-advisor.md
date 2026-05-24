# Plotting Advisor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `plotting-advisor`, a rules-engine skill for Python plotting that recommends chart type, palette, axes, annotation, and accessibility choices, plus a passive lint flow for existing figures.

**Architecture:** A `claude-config/skills/plotting-advisor/` directory containing `SKILL.md` (frontmatter + workflow + inline rules), seven reference docs in `references/`, and three executable scripts in `scripts/` (`palettes.py` for canonical hex constants, `figure_spec.py` for matplotlib `Figure` → JSON extraction, `style_lint.py` as the lint CLI). Tests use stdlib `unittest` (zero new dependencies) and live at `scripts/tests/`. Skill is synced to `~/.claude/` via `sync-config.py push`.

**Tech Stack:** Markdown (skill content), Python 3.9+ (scripts), stdlib `unittest` (tests), Pillow (optional for image-mode lint), `colorspacious` (optional for rigorous colorblind ΔE — falls back to whitelist + RGB distance).

**Spec:** `docs/superpowers/specs/2026-05-24-plotting-advisor-design.md` (committed at `69af12c`).

---

## File Structure

```
claude-config/skills/plotting-advisor/
├── SKILL.md
├── references/
│   ├── SOURCES.md
│   ├── chart-selection.md
│   ├── palettes.md
│   ├── anti-patterns.md
│   ├── accessibility.md
│   ├── interactive-adaptation.md
│   └── scientific-conventions.md
└── scripts/
    ├── palettes.py
    ├── figure_spec.py
    ├── style_lint.py
    └── tests/
        ├── test_palettes.py
        ├── test_figure_spec.py
        └── test_style_lint.py
```

Tests sync along with the skill (sync.config.yaml syncs `skills/` wholesale). Total payload ~30 KB across ~15 files.

---

## Task 1: Create planning entry

**Files:**
- Create: `planning/mac/2026-05-24-plotting-advisor-skill.md` (created by `sync-config.py plan`)

- [ ] **Step 1: Verify clean working tree**

```bash
git status
```

Expected: `nothing to commit, working tree clean` (or only the spec/plan files present).

- [ ] **Step 2: Create planning entry**

```bash
./sync-config.py plan --title "Add plotting-advisor skill"
```

Expected: prints path to new file under `planning/mac/`.

- [ ] **Step 3: Populate the planning entry**

Edit the created file. Replace the template body with:

```markdown
# Add plotting-advisor skill

**Date**: 2026-05-24
**Machine**: mac
**Status**: In Progress

## Objective

Add a new `plotting-advisor` skill (rules engine, advisor pattern) that recommends chart type, palette, axes, annotation, and accessibility choices for Python plotting. Delegates library mechanics to existing `scientific-skills:matplotlib`/`seaborn`/`plotly` skills. Includes passive lint flow (`scripts/style_lint.py`) for checking existing figures.

## Changes Planned

- [ ] Follow CONFIG_MANAGEMENT.md workflow
- [ ] Create `claude-config/skills/plotting-advisor/` skill scaffolding + SKILL.md
- [ ] Write seven reference docs (`SOURCES`, `chart-selection`, `palettes`, `anti-patterns`, `accessibility`, `interactive-adaptation`, `scientific-conventions`)
- [ ] Write three scripts (`palettes.py`, `figure_spec.py`, `style_lint.py`) with unittest coverage
- [ ] Run `./sync-config.py push --dry-run` and `./sync-config.py push`
- [ ] Verify skill appears at `~/.claude/skills/plotting-advisor/` and YAML parses

## Expected Outcome

`plotting-advisor` is invocable, fires on plotting verbs, returns structured checklist + YAML decision card, and provides lint via `python3 scripts/style_lint.py`.

## Actual Outcome

[Filled in after implementation.]

## Assessment

[Filled in after implementation.]

## Related Commits

- [pending]

## Next Steps

- Update `bioinformatician/SKILL.md` to point its missing visualization reference at this skill's `references/scientific-conventions.md` (separate PR).
```

- [ ] **Step 4: Commit the planning entry**

```bash
git add planning/mac/2026-05-24-plotting-advisor-skill.md
git commit -m "$(cat <<'EOF'
chore(planning): plan plotting-advisor skill

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: commit succeeds. Note the SHA for later.

---

## Task 2: Scaffold skill directory + write SKILL.md

**Files:**
- Create: `claude-config/skills/plotting-advisor/SKILL.md`
- Create: `claude-config/skills/plotting-advisor/references/` (empty for now)
- Create: `claude-config/skills/plotting-advisor/scripts/` (empty for now)
- Create: `claude-config/skills/plotting-advisor/scripts/tests/` (empty for now)

- [ ] **Step 1: Create the directory structure**

```bash
mkdir -p claude-config/skills/plotting-advisor/references
mkdir -p claude-config/skills/plotting-advisor/scripts/tests
```

Expected: directories exist (verify with `ls claude-config/skills/plotting-advisor/`).

- [ ] **Step 2: Write SKILL.md**

Create `claude-config/skills/plotting-advisor/SKILL.md` with this exact content:

````markdown
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

# Plotting Advisor

Rules engine for Python plotting. Use BEFORE writing any plotting code. Returns a structured checklist and a paste-ready YAML decision card. Library mechanics are delegated to `scientific-skills:matplotlib`, `seaborn`, or `plotly`.

## When to Use This Skill

Use this skill when:

- A user asks to plot, visualize, chart, or graph any data
- Another agent (`bioinformatician`, `calculator`, `notebook-writer`) is about to render a figure
- A figure already exists and the user asks "is this plot OK?" / "review this figure" / "make this better" → invokes the lint flow
- Intent is documented (compare, distribution, trend, composition, relation). If not, ask one question before recommending.

## When NOT to Use This Skill

Do NOT use this skill when:

- The user only needs `matplotlib`/`seaborn`/`plotly` syntax mechanics (e.g., "how do I set a log axis") → use `scientific-skills:matplotlib`
- The plot is a one-off exploratory glance during debugging (advisor overhead not justified)
- The output is a domain-specific established convention with no flexibility (e.g., a regulatory-required forest plot in clinical reporting) — note the convention and step aside
- The work is non-Python (R/ggplot, D3, etc.) — out of scope

## Workflow

Two flows, sharing the rule base in `references/`.

### Advisor flow (default)

1. **Intake** — gather from the caller:
   - Data shape (variable types, n, dimensionality, missingness)
   - Intent — one of `{compare, distribution, trend, composition, relation, ranking, geographic, network}`
   - Audience (print, slides, paper figure, dashboard)
   - Constraints (caller-specified chart type, palette, B/W output, journal style)
   - If intent is missing, ask **one** question, then proceed.

2. **Chart selection** — open `references/chart-selection.md`, walk the decision tree for (data shape × intent). Output: recommended chart + 1-2 alternatives + the rule that selected it.

3. **Palette selection** — open `references/palettes.md`:
   - Categorical ≤8 → Okabe-Ito
   - Continuous unipolar → viridis / cividis
   - Diverging around meaningful midpoint → ColorBrewer RdBu / PuOr
   - Decorative color → single neutral; reserve color for emphasis

4. **Axis & scale** (inline rules):
   - Linear axis includes zero unless truncation is justified and annotated
   - Log scale only when data spans ≥2 orders of magnitude OR the underlying process is multiplicative
   - ~5-7 major ticks per axis
   - Shared axes for small multiples
   - Date axes redundantly labeled if range crosses years

5. **Annotation & labels** (inline rules):
   - Axis labels with units in parentheses: `Time (h)`
   - Direct-label series when ≤5 (Tufte) instead of legend
   - Title states the finding; subtitle states the context
   - Annotate outliers, intervention points, baselines

6. **Accessibility floor** (inline rules, always apply):
   - Colorblind-safe palette OR redundant encoding (shape / line style / position)
   - Minimum font size: 8pt print / 14pt slide
   - Avoid red/green as the only encoding distinction
   - Test against deuteranopia simulation when categorical color count ≥3

7. **Anti-pattern vetoes** (inline, refuse outright; full list in `references/anti-patterns.md`):
   - No 3D bar/pie/surface
   - No dual y-axes
   - No rainbow palette (`jet`, `gist_rainbow`, `hsv`) on continuous data
   - No truncated baseline without annotation
   - No pie with >5 slices

8. **Output** — produce the structured checklist + decision card (formats below).

### Lint flow (when the caller hands over a figure)

1. **Input** (mutually exclusive modes):
   - `python3 scripts/style_lint.py --image figure.png`
   - `python3 scripts/style_lint.py --figure-spec figure.json` (caller produces JSON via `figure_spec.extract_spec(fig)`)
   - `python3 scripts/style_lint.py --describe "<text description of the figure>"`

2. **Inspection** — `style_lint.py` extracts properties from the input (image / spec / description) and applies the rule base.

3. **Rule check** — violations grouped by severity (critical / major / minor) with one-line fixes and reference anchors.

4. **Exit code**:
   - Default: exit 0 (advisory)
   - `--strict`: exit 2 if any critical violation present

The skill never imports matplotlib to render anything.

## Output Format

### Advisor — Checklist + Decision card

```markdown
## Plotting Advisor: [chart type recommendation, in one phrase]

### Intent & data
- Intent: compare across 4 conditions
- Data: continuous response, ~30 obs/group, balanced, no missing
- Audience: paper figure (300 dpi print)

### Recommended chart: dot plot with median bar
- Rule: Cleveland 1985 — position encodes more accurately than length/angle/area
- Alternatives considered: box plot (loses obs at n=30); violin (overstates smoothness); raincloud (complexity not justified)

### Palette: Okabe-Ito (categorical, 4 colors)
- Hex: #E69F00, #56B4E9, #009E73, #F0E442
- Why: colorblind-safe, ≤8 levels supported, perceptually balanced
- Reference: references/palettes.md#okabe-ito

### Axes & scale
- Y: linear, include zero (response is a count); 5 major ticks
- X: categorical, ordered by condition (control first)
- Units: `Response (counts/min)`

### Annotation
- Direct-label each group (n=4 ≤5 rule)
- Annotate sample size per group below x-axis
- Title states finding; subtitle states context

### Accessibility floor
- ✓ Colorblind-safe palette
- ✓ Redundant encoding (group also encoded by x-position)
- ✓ Font size 9pt (print floor 8pt)
- ✓ No red/green-only distinction

### Anti-pattern check
- ✓ No 3D, no dual axes, no rainbow, no truncated baseline, no pie

### Decision card

​```yaml
chart: dot_plot_with_median
library: seaborn
palette:
  type: categorical
  name: okabe_ito
  hex: ["#E69F00", "#56B4E9", "#009E73", "#F0E442"]
axes:
  x: {type: categorical, order: [control, t1, t2, t3]}
  y: {type: linear, include_zero: true, label: "Response (counts/min)"}
encoding:
  position: condition
  color: condition  # redundant with position — fine
annotation:
  direct_labels: true
  sample_size_per_group: true
  title: "<state finding>"
  subtitle: "<state context>"
accessibility:
  colorblind_safe: true
  min_font_pt: 9
output:
  dpi: 300
  format: pdf
delegate_to: scientific-skills:seaborn
​```
```

Every recommendation includes a rule citation (short tag like `[Cleveland 1985]` or a `references/<file>.md#anchor` pointer). Full bibliography lives in `references/SOURCES.md`.

### Lint output

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

Plus JSON to stdout for machine consumption (see `scripts/style_lint.py --help`).

## Integration with Existing Skills

- **`scientific-skills:matplotlib` / `seaborn` / `plotly`** — explicit `delegate_to` field in the decision card names the right software skill. This advisor never duplicates their syntax content.
- **`bioinformatician`** — its visualization reference can point at `references/scientific-conventions.md` (follow-up PR).
- **`notebook-writer`** — invoke this skill for every plotting cell.
- **`editor` / `consistency-auditor`** — precedent for the rules-engine-as-skill pattern.

## References

- `references/chart-selection.md` — decision tree by intent × data shape (general only)
- `references/palettes.md` — Okabe-Ito, viridis family, ColorBrewer, Tol bright
- `references/anti-patterns.md` — 12-15 anti-patterns with citations
- `references/accessibility.md` — WCAG, colorblind sim, fonts, journal overrides
- `references/interactive-adaptation.md` — how rules shift for plotly/bokeh/altair
- `references/scientific-conventions.md` — volcano, UMAP, Manhattan, forest, ROC, etc.
- `references/SOURCES.md` — full bibliography

## Scripts

- `scripts/palettes.py` — canonical hex lists, `is_colorblind_safe_categorical(hex_list)` helper
- `scripts/figure_spec.py` — `extract_spec(fig)` for matplotlib Figure → JSON; CLI mode also accepts a pickled figure
- `scripts/style_lint.py` — main CLI with three input modes; `--strict` for CI use

Run tests with:

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```
````

- [ ] **Step 3: Validate YAML frontmatter**

```bash
python3 -c "
import yaml, sys
with open('claude-config/skills/plotting-advisor/SKILL.md') as f:
    content = f.read()
assert content.startswith('---\n'), 'missing frontmatter open'
end = content.index('\n---\n', 4)
yaml.safe_load(content[4:end])
print('frontmatter OK')
"
```

Expected: prints `frontmatter OK`.

- [ ] **Step 4: Validate sync would accept the skill**

```bash
./sync-config.py push --dry-run | grep -E "plotting-advisor|ERROR" | head -20
```

Expected: shows `plotting-advisor` as a new skill to add; no `ERROR` lines.

- [ ] **Step 5: Commit**

```bash
git add claude-config/skills/plotting-advisor/SKILL.md
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): scaffold skill with SKILL.md

Establishes the plotting-advisor skill — rules engine for Python
plotting that recommends chart type, palette, axes, annotation,
and accessibility choices. Delegates library mechanics to
scientific-skills:matplotlib / seaborn / plotly.

This commit lays down only SKILL.md; reference docs and scripts
follow in subsequent commits.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: commit succeeds.

---

## Task 3: Write `references/SOURCES.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/SOURCES.md`

- [ ] **Step 1: Write the bibliography**

Create `claude-config/skills/plotting-advisor/references/SOURCES.md`:

```markdown
# Sources

Inline citations across plotting-advisor reference docs use short tags like `[Cleveland 1985]`. Full entries below.

## Core texts

- **[Tufte 1983]** Tufte, E. R. (1983). *The Visual Display of Quantitative Information*. Graphics Press.
- **[Tufte 1990]** Tufte, E. R. (1990). *Envisioning Information*. Graphics Press.
- **[Cleveland 1985]** Cleveland, W. S. (1985). *The Elements of Graphing Data*. Wadsworth.
- **[Cleveland & McGill 1984]** Cleveland, W. S. & McGill, R. (1984). Graphical perception: theory, experimentation, and application to the development of graphical methods. *Journal of the American Statistical Association* 79(387), 531-554.
- **[Wilke 2019]** Wilke, C. O. (2019). *Fundamentals of Data Visualization*. O'Reilly. <https://clauswilke.com/dataviz/>
- **[Healy 2018]** Healy, K. (2018). *Data Visualization: A Practical Introduction*. Princeton University Press. <https://socviz.co/>

## Color and accessibility

- **[Okabe & Ito 2008]** Okabe, M. & Ito, K. (2008). Color universal design — barrier-free color palette. <https://jfly.uni-koeln.de/color/>
- **[Wong 2011]** Wong, B. (2011). Color blindness. *Nature Methods* 8, 441. <https://doi.org/10.1038/nmeth.1618>
- **[Tol 2018]** Tol, P. (2018). Colour schemes. SRON Technical Note. <https://personal.sron.nl/~pault/>
- **[Borland & Taylor 2007]** Borland, D. & Taylor, R. M. (2007). Rainbow color map (still) considered harmful. *IEEE Computer Graphics and Applications* 27(2), 14-17.
- **[WCAG 2.2]** W3C Web Content Accessibility Guidelines 2.2. <https://www.w3.org/TR/WCAG22/>

## Data integrity

- **[Bergstrom & West 2020]** Bergstrom, C. T. & West, J. D. (2020). *Calling Bullshit: The Art of Skepticism in a Data-Driven World*. Random House.
- **[Weissgerber 2015]** Weissgerber, T. L., Milic, N. M., Winham, S. J. & Garovic, V. D. (2015). Beyond bar and line graphs: time for a new data presentation paradigm. *PLOS Biology* 13(4), e1002128.

## Library and methodology

- **[Hunter 2007]** Hunter, J. D. (2007). Matplotlib: a 2D graphics environment. *Computing in Science & Engineering* 9(3), 90-95.
- **[Waskom 2021]** Waskom, M. L. (2021). seaborn: statistical data visualization. *Journal of Open Source Software* 6(60), 3021.
- **[van der Walt et al. 2014]** van der Walt, S. et al. (2014). [viridis colormap]. Released with matplotlib.

## Domain-specific

- **[Ihaka 2003]** Ihaka, R. (2003). Colour for presentation graphics. *Proceedings of the 3rd International Workshop on Distributed Statistical Computing*.
- **[Robinson 1976]** Robinson, A. H. (1976). *Choropleth maps: their design and use*. (For geographic data.)
```

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/SOURCES.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add SOURCES bibliography

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Write `references/chart-selection.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/chart-selection.md`

- [ ] **Step 1: Write the decision tree**

Create the file with this content:

````markdown
# Chart Selection

Walk the tree by **intent**, then by **data shape**. Each leaf names the recommended chart + 1-2 alternatives + the rule that selected it. Domain-specific charts (volcano, UMAP, Manhattan, forest, ROC, etc.) live in `scientific-conventions.md`, not here.

## How to use

1. Identify the intent: `{compare, distribution, trend, composition, relation, ranking, geographic, network}`
2. Identify the data shape: variable types (categorical/numeric/temporal), n, dimensionality
3. Walk to the matching leaf
4. If multiple leaves match, prefer the one whose rule has the strongest perceptual grounding

## Intent: compare

### Categorical × numeric, n_groups ≤ 12, n_obs ≤ 100 → **dot plot**
- Rule: [Cleveland & McGill 1984] — position along a common scale is the most accurate visual encoding
- Alternatives: box plot (loses individual obs, hides bimodality); raincloud (more complex)
- Avoid: bar chart of means (Weissgerber 2015 "dynamite plot" — hides distribution)

### Categorical × numeric, n_obs large (>100/group) → **strip + summary** OR **box plot**
- Rule: [Wilke 2019] — at large n, individual obs overplot; show summary
- Alternatives: violin (use only if shape is the point); ridgeline if many groups

### Categorical × proportion (parts of a whole, ≤2 groups) → **grouped bar** OR **stacked bar**
- Rule: [Cleveland 1985] — bar comparison via length is reliable
- Avoid pie unless ≤5 slices AND comparison is not the point [anti-patterns.md#pie-misuse]

### Paired observations (before/after, treatment/control per subject) → **slopegraph**
- Rule: [Tufte 1990] — directly encodes the change per pair
- Alternative: paired dot plot with connecting lines

### Repeated measures over time, per subject → **spaghetti plot** (line per subject)
- Rule: [Wilke 2019] — preserves within-subject trajectory
- Add summary line (mean ± CI) over the spaghetti

## Intent: distribution

### Single distribution → **histogram with rug** OR **ECDF**
- Rule: [Wilke 2019] — ECDF preserves all data and is parameter-free
- Histogram: pick bin width via Freedman-Diaconis; show the rug for n<200

### ≤4 distributions → **overlaid ECDF** OR **small-multiple histograms**
- Rule: [Tufte 1983] — small multiples preserve comparison
- Avoid overlaid histograms (occlusion)

### >4 distributions → **ridgeline** OR **violin**
- Rule: [Wilke 2019] — ridgeline scales to many groups; violin if shape is the point
- Order ridgelines by a meaningful variable, not alphabet [anti-patterns.md#alphabet-order]

### Bivariate distribution → **hex bin** (large n) OR **scatter with marginal histograms/ECDFs**
- Rule: [Cleveland 1985] — at large n, scatter overplots; hex bin preserves density

## Intent: trend (time/ordered axis)

### Single series → **line plot**
- Rule: [Tufte 1983] — line directly encodes change
- Aspect ratio: bank to 45° for slope reading [Cleveland 1985]
- Always include zero on the y-axis unless truncation is justified and annotated

### Multiple series, ≤5 → **line plot with direct labels**
- Rule: [Tufte 1983] — direct labels beat legend at ≤5 series

### Many series (>5) → **small multiples** (one panel per series)
- Rule: [Tufte 1983] — small multiples > one busy chart
- Shared axes across panels

## Intent: composition

### Single time point → **stacked bar** (proportions) OR **treemap** (hierarchy)
- Rule: [Cleveland 1985] — stacked bar reads composition via length

### Over time → **stacked area** OR **stream graph** (only when ordering carries meaning)
- Rule: [Healy 2018] — stacked area conveys changing composition; avoid when individual series matters more

### Avoid pies → see anti-patterns.md

## Intent: relation

### Two continuous variables → **scatter plot**
- Rule: [Cleveland 1985] — scatter is the canonical bivariate display
- Add: marginal distributions; LOESS or linear fit only if the model is being claimed
- Transparency or jitter when overplotting

### Two continuous, large n → **hex bin** OR **2D density**
- Rule: [Wilke 2019] — preserves density at scales where scatter fails

### Two categorical → **mosaic plot** OR **heatmap of counts/rates**
- Rule: [Healy 2018] — mosaic shows proportional area; heatmap with order

## Intent: ranking

### Ordered categorical with one value per category → **ordered dot plot** OR **ordered bar**
- Rule: [Cleveland 1985] — order by value, not alphabet; dot plot if value range is wide

## Intent: geographic

### Rate or density across regions → **choropleth**
- Rule: [Robinson 1976] — choropleth for rates only; raw counts mislead via region size
- Use a sequential or diverging palette per `palettes.md`

### Counts per location → **proportional symbol** (not choropleth)
- Rule: [Robinson 1976] — symbol area encodes count without region-size confound

## Intent: network

### Sparse, small (n_nodes < 50) → **force-directed graph**
- Rule: [Healy 2018]

### Dense or large → **adjacency matrix** with row/column ordering by clustering
- Rule: [Healy 2018] — matrix scales where node-link does not
````

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/chart-selection.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add chart-selection decision tree

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Write `references/palettes.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/palettes.md`

- [ ] **Step 1: Write the palette catalog**

````markdown
# Palettes

Three palette families. Use the right one for the data type. Hex codes are canonical; `scripts/palettes.py` exports the same constants for the lint script.

## Categorical (discrete levels)

### Okabe-Ito (default categorical, ≤8 levels) {#okabe-ito}

| Index | Hex | Name |
|---|---|---|
| 0 | #000000 | black |
| 1 | #E69F00 | orange |
| 2 | #56B4E9 | sky blue |
| 3 | #009E73 | bluish green |
| 4 | #F0E442 | yellow |
| 5 | #0072B2 | blue |
| 6 | #D55E00 | vermillion |
| 7 | #CC79A7 | reddish purple |

- **Source**: [Okabe & Ito 2008]; popularized by [Wong 2011] *Nature Methods*
- **Why**: deuteranopia / protanopia / tritanopia safe; designed for color universal design
- **Use when**: ≤8 categorical levels; print or screen
- **Avoid when**: levels >8 — collapse categories first; do not extend the palette

### Tol bright (alt categorical, ≤7 levels) {#tol-bright}

| Index | Hex |
|---|---|
| 0 | #4477AA |
| 1 | #EE6677 |
| 2 | #228833 |
| 3 | #CCBB44 |
| 4 | #66CCEE |
| 5 | #AA3377 |
| 6 | #BBBBBB |

- **Source**: [Tol 2018]
- **Why**: print-friendly, qualitatively distinct
- **Use when**: Okabe-Ito's saturation clashes with the medium; ≤7 levels

## Sequential (continuous unipolar)

### viridis family {#viridis}

- **Members**: `viridis`, `cividis`, `plasma`, `magma`, `inferno`
- **Source**: matplotlib (van der Walt et al. 2014)
- **Why**: perceptually uniform; colorblind-safe; monotonic luminance
- **Use when**: continuous unipolar data (intensity, magnitude, density)
- **Prefer cividis** when accessibility to tritanopia matters
- **Avoid for**: diverging data with meaningful midpoint

## Diverging (continuous with meaningful midpoint)

### ColorBrewer RdBu / PuOr / BrBG {#colorbrewer-diverging}

- **Source**: ColorBrewer (Brewer et al.); included in matplotlib
- **Use when**: data has a meaningful midpoint (z-score, log2 fold-change, residual)
- **Pick RdBu** for general use; **PuOr** when red/green sensitivity is a concern
- **Avoid for**: unipolar data — implies false symmetry around zero

## Banned for continuous (anti-pattern)

`jet`, `gist_rainbow`, `hsv`, `rainbow`, `gist_ncar` — non-monotonic luminance creates false visual structure [Borland & Taylor 2007]. See `anti-patterns.md#rainbow-on-continuous`.

## Decorative color

When color encodes **nothing** (e.g., a single trace, a single distribution), use a **single neutral** (medium grey, `#444444`, or matplotlib default blue `#1f77b4`). Reserve color for encoding meaningful variables and for emphasis.

## Quick decision

| Data | Levels / range | Use |
|---|---|---|
| Categorical | ≤8 | Okabe-Ito |
| Categorical | ≤7, print | Tol bright |
| Categorical | >8 | Collapse first; if not possible, ridgeline by group + neutral fills |
| Continuous unipolar | any | viridis (or cividis) |
| Continuous diverging | any | RdBu (or PuOr) |
| No encoding (single series) | n/a | Single neutral |
````

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/palettes.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add palette catalog

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Write `references/anti-patterns.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/anti-patterns.md`

- [ ] **Step 1: Write 14 anti-pattern entries**

````markdown
# Anti-Patterns

Plotting practices to refuse, with the rule that condemns them and the fix.

## Critical (visual lies; refuse outright)

### 3D bar / pie / surface {#3d-charts}
- **Practice**: Add a Z extrusion or perspective to a 2D chart
- **Why wrong**: Occlusion + perspective distortion make values unreadable; angle/area judgments are the worst perceptual encodings [Cleveland & McGill 1984]
- **Fix**: Use the 2D version. Always.

### Dual y-axes {#dual-axes}
- **Practice**: Two series on one chart, each with its own y-axis scale
- **Why wrong**: Implies a relationship by visual alignment that is purely a scaling choice; either axis can be tweaked to invent correlation [Bergstrom & West 2020]
- **Fix**: Two adjacent panels with shared x-axis. Or normalize both series to a common unit.

### Rainbow palette on continuous data {#rainbow-on-continuous}
- **Practice**: Using `jet`, `gist_rainbow`, `hsv`, `rainbow`, or `gist_ncar` to encode a continuous variable
- **Why wrong**: Non-monotonic luminance creates false visual boundaries; not colorblind-safe [Borland & Taylor 2007]
- **Fix**: viridis family (sequential) or RdBu (diverging). See `palettes.md`.

### Truncated baseline without annotation {#truncated-baseline}
- **Practice**: Bar chart or area chart starting at y > 0 with no annotation
- **Why wrong**: Visually amplifies differences. A 2% change reads as 200% [Bergstrom & West 2020]
- **Fix**: Extend axis to zero, OR explicitly annotate the truncation with a broken-axis indicator. For line charts (where baseline conveys less), truncation is acceptable if labeled.

### Pie chart with >5 slices, or any pie when comparison matters {#pie-misuse}
- **Practice**: Pie chart with many slices, or any pie where the goal is comparison
- **Why wrong**: Angle judgments are unreliable [Cleveland & McGill 1984]
- **Fix**: Ordered bar / dot plot

## Major (misleading or hard to read)

### Bar chart of means with error bars hiding distribution ("dynamite plot") {#dynamite-plot}
- **Practice**: Show mean ± SE/SD as a bar with a whisker
- **Why wrong**: Hides distribution shape; small n produces visually identical plots regardless of underlying data; obscures bimodality and outliers [Weissgerber 2015]
- **Fix**: Dot plot, box plot, or violin showing individual observations

### Aspect ratio that flattens slope {#aspect-ratio}
- **Practice**: Line plot stretched horizontally so slopes look flat (or vertically so they look extreme)
- **Why wrong**: Slope perception depends on banking near 45° [Cleveland 1985]
- **Fix**: Bank to 45° — set figure dimensions so the median slope ≈ 1 in the rendered axis

### Legend when direct labels fit {#legend-vs-direct-labels}
- **Practice**: Legend for ≤5 series
- **Why wrong**: Forces eye to ping-pong between legend and data [Tufte 1983]
- **Fix**: Direct-label each series at its terminal point

### Gridlines competing with data {#gridline-noise}
- **Practice**: Dense gridlines, dark color, full opacity
- **Why wrong**: Reduces data-ink ratio; pulls attention from data [Tufte 1983]
- **Fix**: Major gridlines only, light grey (`#dddddd`), behind data

### Alphabet ordering when value-order conveys meaning {#alphabet-order}
- **Practice**: Categorical axis sorted A→Z when the categories have a natural rank
- **Why wrong**: Imposes an irrelevant order; reader cannot read magnitude order from position
- **Fix**: Sort by value (ascending or descending) unless time / experimental order is more meaningful

### Count when rate is the honest comparison {#count-vs-rate}
- **Practice**: Choropleth or bar chart of counts across populations of different sizes
- **Why wrong**: Confounds population size with the variable of interest [Robinson 1976]
- **Fix**: Plot rates (per capita) and note the population denominator

### Overplotting without transparency or jitter {#overplotting}
- **Practice**: Scatter or strip plot with many overlapping points at full opacity
- **Why wrong**: Density is invisible; reader sees only the outermost points
- **Fix**: Alpha 0.3-0.6, or jitter, or hex bin / 2D density at very large n

## Minor (style and convention)

### Heatmap with uninformative row/column order {#unsorted-heatmap}
- **Practice**: Heatmap with rows/columns in input order
- **Why wrong**: Patterns invisible; reader cannot see clusters
- **Fix**: Cluster (hierarchical) or sort by a meaningful aggregate (row sum, marginal value)

### Chartjunk (decorative gradients, drop shadows, icons, 3D bevels) {#chartjunk}
- **Practice**: Visual ornament with no data role
- **Why wrong**: Reduces data-ink ratio [Tufte 1983]
- **Fix**: Strip everything that isn't encoding data

### Stacked bar for time series of many categories {#stacked-bar-time}
- **Practice**: Stacked bar with many slices, repeated across time points
- **Why wrong**: Only the bottom and top categories are readable; middle categories float
- **Fix**: Small multiples per category, or stream graph if the goal is total composition

## See also

- `palettes.md` for palette-specific guidance
- `accessibility.md` for accessibility-specific anti-patterns (low contrast, color-only encoding)
````

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/anti-patterns.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add anti-patterns catalog

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Write `references/accessibility.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/accessibility.md`

- [ ] **Step 1: Write the accessibility deep-dive**

````markdown
# Accessibility

The accessibility floor in `SKILL.md` lists the 5 minimums applied to every plot. This document expands them with specifics and tools.

## Colorblind safety

### What to check

- ~8% of men of European descent and ~0.5% of women have some form of color vision deficiency (CVD)
- The three common forms: deuteranopia, protanopia, tritanopia
- Red/green-only distinctions fail for the majority of CVD types

### How to check

1. **Use a known-safe palette** by default — Okabe-Ito (categorical), viridis / cividis (sequential), RdBu / PuOr (diverging). See `palettes.md`.
2. **Simulate** if a palette is custom:
   - [Coblis](https://www.color-blindness.com/coblis-color-blindness-simulator/) — web tool
   - [daltonlens-python](https://pypi.org/project/daltonlens/) — programmatic
   - `colorspacious` Python library — `colorspacious.cspace_convert(rgb, "sRGB1", colorspacious.CIECAM02UCS)` plus CVD transform

### Redundant encoding

When color encodes a meaningful variable, **add a second channel**:

- Categorical: color + shape (marker style)
- Categorical (lines): color + line style (solid/dashed/dotted)
- Categorical (bars): color + texture or labels
- Continuous (heatmap): color + contour lines

A B/W photocopy test: does the figure still convey the encoded variable when the color is removed?

## Contrast (WCAG)

[WCAG 2.2] minimums for graphical objects:

- **Text** (axis labels, titles): contrast ratio ≥ 4.5:1 against background
- **Non-text graphical objects** (data marks, gridlines vs. background): ≥ 3:1
- **Large text** (≥18pt or ≥14pt bold): ≥ 3:1

Quick check: `colour.contrast` (npm) or any web checker; default matplotlib palette on white background passes.

## Font sizing

| Medium | Minimum body | Minimum minor labels |
|---|---|---|
| Print figure | 8pt | 6pt |
| Paper figure (≥1-column) | 10pt | 8pt |
| Slide | 18pt | 14pt |
| Poster | 24pt | 18pt |

Matplotlib default is 10pt — fine for paper figures, too small for slides. Always set `rcParams["font.size"]` explicitly when the medium isn't paper.

## Journal-specific overrides

| Journal | Font family | Body size | Notes |
|---|---|---|---|
| Nature | Sans-serif (Arial / Helvetica) | 5-7pt minor, 7pt body | Strict size limits, B/W must work |
| Cell | Sans-serif | Similar to Nature | Color OK |
| eLife | Free-form | Typical 8-10pt | More flexible |
| PLOS | Sans-serif preferred | 8pt body | Color OK |
| PNAS | Free-form | 6-8pt body | Color OK |

When the caller names a journal, check its instructions for authors; treat the table above as a starting point.

## Screen readers and figure captions

- Figures embedded in documents need **alt text** describing the data and the finding (not just "Figure 1")
- Caption should state the conclusion the figure supports, in a sentence
- Tables of underlying data should be available for screen-reader users

## B/W fallback test

Before finalizing a color figure for print:

1. Print or convert to greyscale
2. Verify all encoded distinctions are still visible
3. If not, add a second channel (shape, texture, position) per "redundant encoding" above
````

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/accessibility.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add accessibility deep-dive

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Write `references/interactive-adaptation.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/interactive-adaptation.md`

- [ ] **Step 1: Write the interactive section**

````markdown
# Interactive Adaptation

The rules in this skill assume **static** output (matplotlib / seaborn / printed figure). When the caller is using `plotly`, `bokeh`, or `altair`, some rules shift. The accessibility floor never shifts.

## Rules that relax

### Label density

- **Static**: Direct-label series when ≤5; otherwise legend; annotate outliers
- **Interactive**: Hover tooltips replace some on-canvas labels — fewer labels needed
- **Why**: The cost of a tooltip is one mouse hover, not one cognitive saccade
- **Still required**: Title, axis labels, units. Tooltips do not substitute for these.

### Color budget

- **Static**: ≤8 categorical levels (Okabe-Ito limit)
- **Interactive**: Can exceed 8 if hover disambiguates AND there is a search/filter
- **Why**: Reader can interrogate the figure rather than relying on color alone
- **Still required**: Colorblind-safe palette; redundant encoding for the primary distinction (legend / hover identifies precisely)

### Annotation density

- **Static**: Annotate outliers, intervention points, baselines on the canvas
- **Interactive**: Hover surfaces the annotation; on-canvas only for the most important
- **Why**: Static figures must self-explain; interactive figures support exploration

## Rules that shift but do not relax

### "Always include zero on a linear axis"

- **Static**: Include zero, or annotate truncation
- **Interactive (zoom enabled)**: Zoom invalidates the rule — user can zoom to any view
- **What to do**: Annotate the **default view** if it does not include zero. Provide a "reset zoom" affordance.

### Aspect ratio for slope reading

- **Static**: Bank to 45°
- **Interactive**: User can resize; banking is fleeting
- **What to do**: Choose a sensible default that banks the median slope to 45° in the default size

## Rules that do not shift

The **accessibility floor** applies to interactive figures without exception:

- Colorblind-safe palette
- Redundant encoding for meaningful color
- Font size for default rendering
- No red/green-only distinction

## Animation

If the figure animates:

- Respect `prefers-reduced-motion` (CSS media query / browser setting). Disable or simplify animation when set.
- Animation must convey information (e.g., showing change over time) — never decorative
- Always provide a non-animated fallback (static export, frame slider)

## Static export fallback

Every interactive figure should have a sensible static PNG / PDF export for:

- Print
- PDF documentation
- Email / Slack thumbnails

The static export should apply all static-figure rules. If the interactive figure cannot be expressed as a single static frame (e.g., a multi-panel exploratory view), provide a small set of static panels covering the main views.

## When to use this skill vs. just rendering interactive

Even when the output is interactive, use this skill to choose:
- Chart type (apply the chart-selection rules)
- Palette
- Encoding strategy
- Default view (what the reader sees before interacting)

Then delegate library mechanics to `scientific-skills:plotly` (or bokeh/altair equivalent).
````

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/interactive-adaptation.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add interactive-adaptation notes

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Write `references/scientific-conventions.md`

**Files:**
- Create: `claude-config/skills/plotting-advisor/references/scientific-conventions.md`

- [ ] **Step 1: Write the scientific conventions**

````markdown
# Scientific Conventions

Domain-specific chart conventions. Each entry: the convention, the canonical layout, common failure modes. Use this when the caller's intent involves one of these established forms.

## Volcano plot (differential expression)

- **Layout**: x = log2 fold-change; y = -log10(adjusted p-value); points colored / highlighted by significance threshold
- **Palette**: Grey for non-significant; Okabe-Ito orange/blue for up/down regulated
- **Annotations**: Top-N gene names; horizontal line at significance threshold (-log10 p = 1.3 for p=0.05, or per FDR); vertical lines at fold-change threshold
- **Common failures**: No threshold lines (reader cannot judge significance); too many gene labels (text occlusion); rainbow palette
- **Citation**: Convention from microarray era, now standard in RNA-seq

## Dimensionality reduction (UMAP / t-SNE / PCA)

- **Layout**: 2D scatter; points colored by a categorical variable (cluster, condition, cell type)
- **Palette**: Okabe-Ito for ≤8 clusters; for more, use distinct hues with redundant labels
- **Annotations**: Cluster labels placed at cluster centroid; axes labeled (e.g., `UMAP 1`, `UMAP 2`)
- **Common failures**:
  - Reading distances as meaningful (t-SNE/UMAP distances are not preserved; PCA distances are)
  - Rainbow palette
  - Missing color legend
  - Tiny points at high n — use alpha + small marker size
- **Rule**: Always state which method (UMAP vs t-SNE vs PCA) and parameters (perplexity, n_neighbors, min_dist) in the caption

## Manhattan plot (GWAS)

- **Layout**: x = chromosomal position (chromosomes concatenated); y = -log10(p-value); alternating colors by chromosome
- **Palette**: Two alternating muted colors (e.g., `#2c7fb8` and `#7fcdbb`) — purely visual separation, not encoding
- **Annotations**: Genome-wide significance line (-log10 p ≈ 7.3 for p = 5×10⁻⁸); suggestive line; top SNP rsIDs
- **Common failures**: Color encoding chromosome (it's just a separator, not a variable); missing significance line; y-axis truncation

## Forest plot (meta-analysis / clinical)

- **Layout**: One row per study/group; x = effect size with confidence interval; vertical line at null effect (1.0 for OR/RR, 0 for differences)
- **Annotations**: Study names; sample size; effect size + CI as text; summary diamond at bottom (if meta-analysis)
- **Common failures**: Asymmetric CIs not visually distinguished from symmetric; weighting not encoded (use point size); missing summary

## Heatmap with hierarchical clustering (omics)

- **Layout**: Cells colored by value; rows and columns reordered by hierarchical clustering; dendrograms shown
- **Palette**: Diverging (RdBu) if values are z-scored or log-fold; sequential (viridis) if values are unipolar (e.g., raw expression after log)
- **Annotations**: Row/column dendrograms; color bars for sample metadata; legend with scale
- **Common failures**:
  - Unsorted heatmap (see `anti-patterns.md#unsorted-heatmap`)
  - Wrong palette for data type (sequential on z-scores; diverging on unipolar)
  - Missing scale legend
  - Color saturation at 99th percentile not stated (outliers wash out the rest)
- **Rule**: When using diverging palette, center it on zero (`vmin=-x, vmax=+x, center=0`)

## Kaplan-Meier survival curve

- **Layout**: Stepped line per group; x = time; y = survival probability (0-1); censored events marked with ticks
- **Annotations**: Risk table below x-axis (numbers at risk per group at each time point); log-rank p-value; group labels at line ends (direct labels)
- **Common failures**: No risk table; no censored ticks; missing CI bands; log-rank p without specifying which groups compared

## ROC curve

- **Layout**: x = false positive rate; y = true positive rate; diagonal reference line; one curve per classifier
- **Annotations**: AUC value per curve in legend; chance line; ≤5 classifiers → direct labels
- **Common failures**: Missing chance diagonal; AUC reported without confidence interval; comparing classifiers with different test sets

## Bland-Altman plot (agreement between methods)

- **Layout**: x = mean of two measurements; y = difference between them; horizontal lines at mean difference and ±1.96 SD (limits of agreement)
- **Annotations**: Limits-of-agreement lines and their values; bias line; CI on bias
- **Common failures**: Showing correlation plot instead (correlation tests association, not agreement); missing limits of agreement

## Raincloud plot

- **Layout**: Half-violin (distribution) + dot strip (individual obs) + box/median bar (summary), per group
- **Palette**: Okabe-Ito for groups
- **Use when**: Comparing distributions where shape, individual obs, AND summary all matter
- **Common failures**: Overcrowded — use for ≤4 groups; otherwise simpler form

## Notes

When the caller's request matches one of these conventions:
1. Recognize it explicitly in the checklist ("Recommended chart: volcano plot")
2. Apply the convention's standard layout
3. Apply the general rules (palette, axes, annotation, accessibility) on top
4. Cite this file in the checklist

If the convention conflicts with a general rule (rare), the convention wins, but the conflict should be noted in the checklist.
````

- [ ] **Step 2: Commit**

```bash
git add claude-config/skills/plotting-advisor/references/scientific-conventions.md
git commit -m "$(cat <<'EOF'
docs(plotting-advisor): add scientific conventions catalog

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Write `scripts/palettes.py` with tests

**Files:**
- Create: `claude-config/skills/plotting-advisor/scripts/palettes.py`
- Create: `claude-config/skills/plotting-advisor/scripts/tests/test_palettes.py`

- [ ] **Step 1: Write the failing test**

Create `claude-config/skills/plotting-advisor/scripts/tests/test_palettes.py`:

```python
"""Tests for palettes.py — canonical hex constants and colorblind-safety helper."""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import palettes  # noqa: E402


class TestPaletteConstants(unittest.TestCase):
    def test_okabe_ito_has_8_colors(self):
        self.assertEqual(len(palettes.OKABE_ITO), 8)

    def test_okabe_ito_all_hex(self):
        for c in palettes.OKABE_ITO:
            self.assertRegex(c, r"^#[0-9A-Fa-f]{6}$")

    def test_tol_bright_has_7_colors(self):
        self.assertEqual(len(palettes.TOL_BRIGHT), 7)

    def test_rainbow_names_includes_jet(self):
        self.assertIn("jet", palettes.RAINBOW_NAMES)
        self.assertIn("gist_rainbow", palettes.RAINBOW_NAMES)
        self.assertIn("hsv", palettes.RAINBOW_NAMES)

    def test_viridis_names_includes_viridis_and_cividis(self):
        self.assertIn("viridis", palettes.VIRIDIS_NAMES)
        self.assertIn("cividis", palettes.VIRIDIS_NAMES)


class TestColorblindSafeHelper(unittest.TestCase):
    def test_okabe_ito_is_safe(self):
        self.assertTrue(palettes.is_colorblind_safe_categorical(palettes.OKABE_ITO))

    def test_red_green_pair_is_not_safe(self):
        # Pure red + pure green is the classic red/green confusion case
        self.assertFalse(
            palettes.is_colorblind_safe_categorical(["#FF0000", "#00FF00"])
        )

    def test_single_color_is_trivially_safe(self):
        self.assertTrue(palettes.is_colorblind_safe_categorical(["#444444"]))

    def test_empty_palette_is_trivially_safe(self):
        self.assertTrue(palettes.is_colorblind_safe_categorical([]))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 -m unittest claude-config.skills.plotting-advisor.scripts.tests.test_palettes -v
```

Expected: `ModuleNotFoundError: No module named 'palettes'` (palettes.py doesn't exist yet).

- [ ] **Step 3: Write `palettes.py`**

Create `claude-config/skills/plotting-advisor/scripts/palettes.py`:

```python
"""Canonical color palettes for plotting-advisor.

Hex constants mirror what is documented in references/palettes.md.
The colorblind-safety helper uses a whitelist of known-safe palettes
plus a coarse RGB-distance fallback. For rigorous CIELAB ΔE under
CVD simulation, callers should use the `colorspacious` library
directly — this module deliberately keeps zero non-stdlib dependencies.
"""
from __future__ import annotations

# Okabe-Ito (Okabe & Ito 2008; Wong 2011 Nature Methods)
OKABE_ITO = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]

# Tol bright (Tol 2018)
TOL_BRIGHT = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
]

# Perceptually uniform sequential maps (matplotlib built-ins)
VIRIDIS_NAMES = frozenset({"viridis", "cividis", "plasma", "magma", "inferno"})

# ColorBrewer diverging (commonly used names available in matplotlib)
COLORBREWER_DIVERGING = frozenset(
    {"RdBu", "PuOr", "BrBG", "RdYlBu", "RdYlGn", "Spectral"}
)

# Banned for continuous data (non-monotonic luminance — Borland & Taylor 2007)
RAINBOW_NAMES = frozenset(
    {"jet", "gist_rainbow", "hsv", "rainbow", "gist_ncar", "nipy_spectral"}
)

# Palettes the lint script trusts without re-checking
_KNOWN_SAFE_CATEGORICAL = (
    tuple(OKABE_ITO),
    tuple(TOL_BRIGHT),
)


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    if len(h) != 6:
        raise ValueError(f"not a 6-digit hex color: {h!r}")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _rgb_distance(a: tuple[int, int, int], b: tuple[int, int, int]) -> float:
    """Euclidean RGB distance. Coarse but cheap; used only as a fallback."""
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


def is_colorblind_safe_categorical(
    hex_list: list[str],
    *,
    rgb_distance_tolerance: float = 80.0,
) -> bool:
    """Return True if a categorical palette is colorblind-safe.

    Strategy (in order):
    1. Trivially safe if ≤1 color.
    2. Whitelist match against known-safe palettes (Okabe-Ito, Tol bright).
    3. Fallback: every pair must have RGB distance ≥ tolerance.
       This is a coarse approximation — it catches obviously bad pairs
       (pure red vs pure green at full saturation = distance ≈ 360,
       BUT under deuteranopia they appear similar). The tolerance is
       deliberately permissive of low-saturation pairs; the lint
       script's primary path is the whitelist, with this fallback
       only flagging clear failures.

    For rigorous checks, callers should use `colorspacious` with a
    CVD transform; this helper is for the zero-dependency advisory
    case.
    """
    if len(hex_list) <= 1:
        return True

    normalized = tuple(c.upper() for c in hex_list)
    for safe in _KNOWN_SAFE_CATEGORICAL:
        safe_upper = tuple(c.upper() for c in safe)
        # Subset check — any subset of a known-safe palette is safe
        if all(c in safe_upper for c in normalized):
            return True

    # Fallback: pairwise distance + red/green saturation check
    rgbs = [_hex_to_rgb(c) for c in hex_list]
    for i in range(len(rgbs)):
        for j in range(i + 1, len(rgbs)):
            ri, gi, bi = rgbs[i]
            rj, gj, bj = rgbs[j]
            # Flag explicit red/green-only confusion: one color
            # dominated by R channel, the other by G channel, both
            # at high saturation.
            high_sat_red_i = ri > 200 and gi < 80 and bi < 80
            high_sat_green_j = gj > 200 and rj < 80 and bj < 80
            high_sat_red_j = rj > 200 and gj < 80 and bj < 80
            high_sat_green_i = gi > 200 and ri < 80 and bi < 80
            if (high_sat_red_i and high_sat_green_j) or (
                high_sat_red_j and high_sat_green_i
            ):
                return False
            if _rgb_distance(rgbs[i], rgbs[j]) < rgb_distance_tolerance:
                return False
    return True
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: all tests pass (8 tests, 0 failures).

- [ ] **Step 5: Commit**

```bash
git add claude-config/skills/plotting-advisor/scripts/palettes.py \
        claude-config/skills/plotting-advisor/scripts/tests/test_palettes.py
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): add palettes.py with colorblind-safety helper

Canonical hex constants (Okabe-Ito, Tol bright, viridis names,
ColorBrewer diverging, rainbow banned-list) and a zero-dependency
is_colorblind_safe_categorical() helper using whitelist + RGB
distance fallback. unittest coverage included.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Write `scripts/figure_spec.py` with tests

**Files:**
- Create: `claude-config/skills/plotting-advisor/scripts/figure_spec.py`
- Create: `claude-config/skills/plotting-advisor/scripts/tests/test_figure_spec.py`

This script does **not** require matplotlib at import time — it only uses matplotlib in `extract_spec(fig)` if the caller passes a Figure. Tests use a hand-built mock object (duck-typed) so they run without matplotlib.

- [ ] **Step 1: Write the failing test**

Create `claude-config/skills/plotting-advisor/scripts/tests/test_figure_spec.py`:

```python
"""Tests for figure_spec.py — matplotlib Figure → JSON spec.

Tests use a minimal duck-typed mock so they run without matplotlib.
"""
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import figure_spec  # noqa: E402


def make_mock_axis(
    *,
    xlabel="x",
    ylabel="y",
    title="",
    xscale="linear",
    yscale="linear",
    xlim=(0.0, 10.0),
    ylim=(0.0, 100.0),
    line_colors=None,
    has_legend=False,
):
    ax = mock.MagicMock()
    ax.get_xlabel.return_value = xlabel
    ax.get_ylabel.return_value = ylabel
    ax.get_title.return_value = title
    ax.get_xscale.return_value = xscale
    ax.get_yscale.return_value = yscale
    ax.get_xlim.return_value = xlim
    ax.get_ylim.return_value = ylim
    ax.name = "rectilinear"  # not 3D
    # Lines
    lines = []
    for c in (line_colors or []):
        line = mock.MagicMock()
        line.get_color.return_value = c
        lines.append(line)
    ax.get_lines.return_value = lines
    # Legend
    ax.get_legend.return_value = mock.MagicMock() if has_legend else None
    # Patches (bars / pie wedges) — empty by default
    ax.patches = []
    # Images (heatmap content) — empty by default
    ax.images = []
    return ax


class TestExtractSpec(unittest.TestCase):
    def test_extracts_axis_labels(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [
            make_mock_axis(xlabel="Time (h)", ylabel="Response (counts/min)")
        ]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["xlabel"], "Time (h)")
        self.assertEqual(spec["axes"][0]["ylabel"], "Response (counts/min)")

    def test_extracts_scale(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(yscale="log")]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["yscale"], "log")

    def test_extracts_ylim_zero_inclusion(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(ylim=(42.0, 51.0))]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["ylim"], [42.0, 51.0])
        self.assertFalse(spec["axes"][0]["ylim_includes_zero"])

    def test_extracts_line_colors(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [
            make_mock_axis(line_colors=["#E69F00", "#56B4E9"])
        ]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["line_colors"], ["#E69F00", "#56B4E9"])

    def test_detects_legend(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(has_legend=True)]
        spec = figure_spec.extract_spec(fig)
        self.assertTrue(spec["axes"][0]["has_legend"])

    def test_n_axes(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(), make_mock_axis()]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["n_axes"], 2)

    def test_detects_3d_projection(self):
        fig = mock.MagicMock()
        ax = make_mock_axis()
        ax.name = "3d"
        fig.get_axes.return_value = [ax]
        spec = figure_spec.extract_spec(fig)
        self.assertTrue(spec["axes"][0]["is_3d"])


class TestCliMode(unittest.TestCase):
    def test_cli_with_pickle(self):
        # Build a mock Figure, pickle it, and verify the CLI reads it back.
        import pickle
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(xlabel="t", ylabel="y")]
        with tempfile.NamedTemporaryFile(
            suffix=".pkl", delete=False
        ) as pkl:
            pickle.dump(fig, pkl)
            pkl_path = pkl.name
        with tempfile.NamedTemporaryFile(
            suffix=".json", delete=False
        ) as out:
            out_path = out.name
        try:
            rc = figure_spec.cli_main(
                ["--pickle", pkl_path, "--out", out_path]
            )
            self.assertEqual(rc, 0)
            with open(out_path) as f:
                data = json.load(f)
            self.assertEqual(data["axes"][0]["xlabel"], "t")
        finally:
            os.unlink(pkl_path)
            os.unlink(out_path)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: `ModuleNotFoundError: No module named 'figure_spec'`.

- [ ] **Step 3: Write `figure_spec.py`**

Create `claude-config/skills/plotting-advisor/scripts/figure_spec.py`:

```python
"""Extract a JSON-serializable spec from a matplotlib Figure.

Usable two ways:
  - Import: `from figure_spec import extract_spec; extract_spec(fig)`
  - CLI:    `python3 figure_spec.py --pickle fig.pkl --out fig.json`

The lint script (style_lint.py) consumes the JSON spec via its
--figure-spec mode. This module is designed to be matplotlib-aware
but matplotlib-optional: extract_spec uses duck-typed Figure API
(.get_axes(), ax.get_xlabel(), etc.), so it works on any object
that quacks like a Figure.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys


def extract_spec(fig) -> dict:
    """Walk fig.axes and return a JSON-serializable dict.

    Captures: per-axis labels, scales, limits, line colors, legend
    presence, 3D projection flag, bar/patch count, image presence.
    Does NOT capture: tick labels (variable), font sizes (need
    rcParams), or any image pixels.
    """
    axes_spec = []
    for ax in fig.get_axes():
        xlim = tuple(ax.get_xlim())
        ylim = tuple(ax.get_ylim())
        line_colors = [_normalize_color(line.get_color()) for line in ax.get_lines()]
        axes_spec.append(
            {
                "xlabel": ax.get_xlabel() or "",
                "ylabel": ax.get_ylabel() or "",
                "title": ax.get_title() or "",
                "xscale": ax.get_xscale(),
                "yscale": ax.get_yscale(),
                "xlim": [float(xlim[0]), float(xlim[1])],
                "ylim": [float(ylim[0]), float(ylim[1])],
                "xlim_includes_zero": xlim[0] <= 0.0 <= xlim[1],
                "ylim_includes_zero": ylim[0] <= 0.0 <= ylim[1],
                "is_3d": getattr(ax, "name", "rectilinear") == "3d",
                "line_colors": line_colors,
                "has_legend": ax.get_legend() is not None,
                "n_patches": len(getattr(ax, "patches", []) or []),
                "n_images": len(getattr(ax, "images", []) or []),
            }
        )
    return {
        "n_axes": len(axes_spec),
        "axes": axes_spec,
    }


def _normalize_color(c) -> str:
    """Best-effort color → hex string. Returns '?' if unrecognized."""
    if isinstance(c, str):
        if c.startswith("#"):
            return c.upper()
        # Common named colors that matplotlib accepts
        named = {
            "k": "#000000", "black": "#000000",
            "w": "#FFFFFF", "white": "#FFFFFF",
            "r": "#FF0000", "red": "#FF0000",
            "g": "#008000", "green": "#008000",
            "b": "#0000FF", "blue": "#0000FF",
            "c": "#00FFFF", "cyan": "#00FFFF",
            "m": "#FF00FF", "magenta": "#FF00FF",
            "y": "#FFFF00", "yellow": "#FFFF00",
        }
        return named.get(c.lower(), "?")
    if isinstance(c, (tuple, list)) and len(c) >= 3:
        r, g, b = c[:3]
        # matplotlib uses 0-1 floats; convert to 0-255 hex
        if all(isinstance(v, float) and 0.0 <= v <= 1.0 for v in (r, g, b)):
            return "#{:02X}{:02X}{:02X}".format(
                int(round(r * 255)), int(round(g * 255)), int(round(b * 255))
            )
    return "?"


def cli_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract JSON spec from a pickled matplotlib Figure."
    )
    parser.add_argument(
        "--pickle", required=True, help="Path to a pickled matplotlib Figure"
    )
    parser.add_argument(
        "--out", required=True, help="Path to write the JSON spec"
    )
    args = parser.parse_args(argv)
    try:
        with open(args.pickle, "rb") as f:
            fig = pickle.load(f)
    except Exception as e:
        print(f"error: could not load pickle: {e}", file=sys.stderr)
        return 1
    spec = extract_spec(fig)
    with open(args.out, "w") as f:
        json.dump(spec, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(cli_main())
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: 8 (palettes) + 8 (figure_spec) = 16 tests, all pass.

- [ ] **Step 5: Commit**

```bash
git add claude-config/skills/plotting-advisor/scripts/figure_spec.py \
        claude-config/skills/plotting-advisor/scripts/tests/test_figure_spec.py
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): add figure_spec.py extractor

Matplotlib Figure → JSON spec. Usable as import or CLI. Captures
per-axis labels, scales, limits, line colors, legend presence,
3D flag, and bar/image counts — the inputs the lint script needs.
unittest coverage with duck-typed mocks (no matplotlib required
at test time).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Write `scripts/style_lint.py` — describe mode + critical rules

**Files:**
- Create: `claude-config/skills/plotting-advisor/scripts/style_lint.py`
- Create: `claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py`

The lint script is the largest piece. We build it incrementally across Tasks 12-15:

- Task 12 (this task): describe-mode infrastructure + critical rules (5 checks)
- Task 13: major rules (figure-spec mode)
- Task 14: minor rules + image mode skeleton
- Task 15: image mode body + `--strict` flag + integration

- [ ] **Step 1: Write the failing test for describe mode + critical rules**

Create `claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py`:

```python
"""Tests for style_lint.py — describe mode + critical-severity rules.

Build this file out incrementally across tasks 12-15.
"""
import json
import os
import subprocess
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import style_lint  # noqa: E402


def _run_describe(text: str) -> list[dict]:
    """Run the lint engine on a describe string, return parsed violations."""
    return style_lint.lint_describe(text)


def _severities(violations: list[dict]) -> list[str]:
    return [v["severity"] for v in violations]


def _rules(violations: list[dict]) -> list[str]:
    return [v["rule"] for v in violations]


class TestDescribeModeCritical(unittest.TestCase):
    def test_rainbow_on_continuous_flagged(self):
        v = _run_describe(
            "scatter plot of two continuous variables colored by a third "
            "continuous using jet colormap"
        )
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#rainbow-on-continuous", _rules(v))

    def test_3d_chart_flagged(self):
        v = _run_describe("3D bar chart of 4 categories with z axis as count")
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#3d-charts", _rules(v))

    def test_dual_axis_flagged(self):
        v = _run_describe(
            "line chart with dual y-axes, one for revenue and one for clicks"
        )
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#dual-axes", _rules(v))

    def test_truncated_baseline_flagged(self):
        v = _run_describe(
            "bar chart of company revenue, y-axis starts at 42 to highlight differences"
        )
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#truncated-baseline", _rules(v))

    def test_missing_axis_labels_flagged(self):
        v = _run_describe("scatter plot with no axis labels")
        self.assertIn("critical", _severities(v))
        self.assertIn("axes#labels-required", _rules(v))

    def test_clean_description_no_critical(self):
        v = _run_describe(
            "dot plot of 4 conditions using Okabe-Ito palette, "
            "y-axis labeled 'Response (counts/min)' from 0 to 100, "
            "x-axis labeled 'Condition'"
        )
        self.assertNotIn("critical", _severities(v))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: `ModuleNotFoundError: No module named 'style_lint'`.

- [ ] **Step 3: Write the describe-mode infrastructure + 5 critical rules**

Create `claude-config/skills/plotting-advisor/scripts/style_lint.py`:

```python
"""Plotting Advisor lint — passive checker for figures and descriptions.

Three input modes (mutually exclusive):
  --describe TEXT       : textual description of a figure
  --figure-spec PATH    : JSON spec produced by figure_spec.extract_spec(fig)
  --image PATH          : a saved PNG/PDF (image-mode lint, optional Pillow)

Default: exit 0 (advisory). --strict exits 2 on any critical violation.

This module is built incrementally:
  Task 12: describe mode + critical rules (this task)
  Task 13: major rules
  Task 14: minor rules + figure-spec mode
  Task 15: image mode + --strict flag wiring
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field, asdict

# Lazy imports: palettes is sibling
THIS_DIR = __file__.rsplit("/", 1)[0] if "/" in __file__ else "."
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)
import palettes  # noqa: E402


SEVERITIES = ("critical", "major", "minor")


@dataclass
class Violation:
    severity: str
    rule: str
    where: str = ""
    observed: object = None
    fix: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        # Drop None observed for clean JSON
        if d["observed"] is None:
            d.pop("observed")
        return d


@dataclass
class LintReport:
    violations: list[Violation] = field(default_factory=list)

    def add(self, v: Violation) -> None:
        self.violations.append(v)

    def critical(self) -> list[Violation]:
        return [v for v in self.violations if v.severity == "critical"]

    def major(self) -> list[Violation]:
        return [v for v in self.violations if v.severity == "major"]

    def minor(self) -> list[Violation]:
        return [v for v in self.violations if v.severity == "minor"]


# ----- Describe mode -----

_RAINBOW_PATTERN = re.compile(
    r"\b(jet|gist[\s_-]?rainbow|hsv|rainbow|gist[\s_-]?ncar|nipy[\s_-]?spectral)\b",
    re.IGNORECASE,
)
_3D_PATTERN = re.compile(r"\b3[\s-]?d\b|three[\s-]?dimensional", re.IGNORECASE)
_DUAL_AXIS_PATTERN = re.compile(
    r"dual[\s-]?(y[\s-]?)?ax(e|i)s|two[\s-]?y[\s-]?ax(e|i)s|secondary[\s-]?y",
    re.IGNORECASE,
)
_TRUNCATED_PATTERN = re.compile(
    r"(y[\s-]?axis\s+starts?\s+at\s+(?!0)\d|truncat\w+\s+(y[\s-]?)?ax|baseline\s+at\s+(?!0)\d)",
    re.IGNORECASE,
)
_NO_LABELS_PATTERN = re.compile(
    r"no\s+(axis|axes)\s+labels?|axes?\s+(are\s+)?unlabel\w+|missing\s+(axis|axes)\s+labels?",
    re.IGNORECASE,
)
_CONTINUOUS_HINT_PATTERN = re.compile(
    r"continuous|gradient|colormap|color\s*scale|heatmap|density|magnitude",
    re.IGNORECASE,
)


def lint_describe(text: str) -> list[dict]:
    """Run all describe-mode checks; return list of violation dicts."""
    report = LintReport()

    # --- Critical rules ---

    # Rainbow on continuous
    if _RAINBOW_PATTERN.search(text) and _CONTINUOUS_HINT_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#rainbow-on-continuous",
                where="palette",
                fix=(
                    "Use viridis / cividis (sequential) or RdBu / PuOr "
                    "(diverging). See references/palettes.md."
                ),
            )
        )

    # 3D charts
    if _3D_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#3d-charts",
                where="chart_type",
                fix="Use the 2D version. Always.",
            )
        )

    # Dual y-axes
    if _DUAL_AXIS_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#dual-axes",
                where="axes",
                fix=(
                    "Use two adjacent panels with shared x-axis, OR "
                    "normalize both series to a common unit."
                ),
            )
        )

    # Truncated baseline (look for explicit non-zero start without annotation)
    if _TRUNCATED_PATTERN.search(text) and not re.search(
        r"annotated|broken[\s-]?axis|labeled", text, re.IGNORECASE
    ):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#truncated-baseline",
                where="y_axis",
                fix=(
                    "Extend axis to zero, OR annotate the truncation with "
                    "a broken-axis indicator."
                ),
            )
        )

    # Missing axis labels
    if _NO_LABELS_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="axes#labels-required",
                where="axes",
                fix="Add labels with units in parentheses, e.g. 'Time (h)'.",
            )
        )

    # Major and minor rules added in later tasks
    return [v.to_dict() for v in report.violations]


# ----- Render output (used by CLI) -----

def render_markdown(violations: list[dict]) -> str:
    """Render violations as a markdown report (matches SKILL.md spec)."""
    by_sev = {"critical": [], "major": [], "minor": []}
    for v in violations:
        by_sev.setdefault(v["severity"], []).append(v)
    total = len(violations)
    lines = [f"## Plotting Advisor: figure check — {total} issue{'s' if total != 1 else ''} found", ""]
    for sev in SEVERITIES:
        items = by_sev.get(sev, [])
        if not items:
            continue
        lines.append(f"### {sev.capitalize()} ({len(items)})")
        for v in items:
            head = v.get("rule", "<unknown>")
            where = v.get("where", "")
            fix = v.get("fix", "")
            line1 = f"- **{head}**"
            if where:
                line1 += f" — at `{where}`"
            lines.append(line1)
            if fix:
                lines.append(f"  - Fix: {fix}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


# ----- CLI -----

def cli_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plotting Advisor lint — passive figure checker."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--describe", help="Textual description of the figure")
    group.add_argument(
        "--figure-spec", help="Path to JSON spec from figure_spec.extract_spec()"
    )
    group.add_argument("--image", help="Path to a saved PNG/PDF")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit 2 on any critical violation (default: always exit 0)",
    )
    parser.add_argument(
        "--json", action="store_true", help="Output JSON to stdout instead of markdown"
    )
    args = parser.parse_args(argv)

    if args.describe:
        violations = lint_describe(args.describe)
    elif args.figure_spec:
        # Task 13/14 will implement this.
        print("error: --figure-spec mode not yet implemented", file=sys.stderr)
        return 1
    elif args.image:
        # Task 15 will implement this.
        print("error: --image mode not yet implemented", file=sys.stderr)
        return 1
    else:
        parser.print_help()
        return 1

    if args.json:
        print(json.dumps(violations, indent=2))
    else:
        print(render_markdown(violations), end="")

    if args.strict and any(v["severity"] == "critical" for v in violations):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(cli_main())
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: 16 (prev) + 6 (describe critical) = 22 tests, all pass.

- [ ] **Step 5: Smoke-test the CLI**

```bash
python3 claude-config/skills/plotting-advisor/scripts/style_lint.py \
  --describe "3D bar chart of revenue using jet colormap, no axis labels"
```

Expected: markdown output with 3 critical violations listed (3d-charts, rainbow-on-continuous, labels-required).

- [ ] **Step 6: Commit**

```bash
git add claude-config/skills/plotting-advisor/scripts/style_lint.py \
        claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): style_lint describe-mode + critical rules

Five critical-severity checks (rainbow on continuous, 3D charts,
dual y-axes, truncated baseline, missing axis labels) plus the
CLI scaffolding, Violation/LintReport dataclasses, and markdown/
JSON output rendering. Major and minor rules + figure-spec/image
modes follow in the next tasks.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: Add major rules + figure-spec mode

**Files:**
- Modify: `claude-config/skills/plotting-advisor/scripts/style_lint.py`
- Modify: `claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py`

- [ ] **Step 1: Add failing tests for major rules and figure-spec mode**

Append the following classes to `test_style_lint.py` (before the `if __name__ == "__main__":` line):

```python
class TestDescribeModeMajor(unittest.TestCase):
    def test_pie_too_many_slices(self):
        v = _run_describe("pie chart with 8 slices showing market share")
        self.assertIn("anti-patterns#pie-misuse", _rules(v))
        # severity major
        for vv in v:
            if vv["rule"] == "anti-patterns#pie-misuse":
                self.assertEqual(vv["severity"], "major")

    def test_dynamite_plot_flagged(self):
        v = _run_describe(
            "bar chart of means with error bars showing standard error "
            "across 5 conditions"
        )
        self.assertIn("anti-patterns#dynamite-plot", _rules(v))

    def test_red_green_only_palette_flagged(self):
        v = _run_describe(
            "scatter plot of two groups colored red and green"
        )
        self.assertIn("accessibility#red-green-only", _rules(v))

    def test_alphabet_order_flagged(self):
        v = _run_describe(
            "bar chart of revenue per state sorted alphabetically"
        )
        self.assertIn("anti-patterns#alphabet-order", _rules(v))


class TestFigureSpecMode(unittest.TestCase):
    def _lint_spec(self, spec: dict) -> list[dict]:
        return style_lint.lint_figure_spec(spec)

    def test_missing_axis_labels_critical(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "",
                    "ylabel": "",
                    "title": "Sales",
                    "xscale": "linear",
                    "yscale": "linear",
                    "xlim": [0, 10],
                    "ylim": [0, 100],
                    "xlim_includes_zero": True,
                    "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": ["#1f77b4"],
                    "has_legend": False,
                    "n_patches": 0,
                    "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("axes#labels-required", _rules(v))

    def test_3d_critical(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x",
                    "ylabel": "y",
                    "title": "",
                    "xscale": "linear",
                    "yscale": "linear",
                    "xlim": [0, 10],
                    "ylim": [0, 100],
                    "xlim_includes_zero": True,
                    "ylim_includes_zero": True,
                    "is_3d": True,
                    "line_colors": [],
                    "has_legend": False,
                    "n_patches": 0,
                    "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("anti-patterns#3d-charts", _rules(v))

    def test_truncated_baseline_critical(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x",
                    "ylabel": "y",
                    "title": "",
                    "xscale": "linear",
                    "yscale": "linear",
                    "xlim": [0, 10],
                    "ylim": [42, 51],
                    "xlim_includes_zero": True,
                    "ylim_includes_zero": False,
                    "is_3d": False,
                    "line_colors": [],
                    "has_legend": False,
                    "n_patches": 4,  # bars present
                    "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("anti-patterns#truncated-baseline", _rules(v))

    def test_dual_axes_critical(self):
        spec = {"n_axes": 2, "axes": [
            {"xlabel": "t", "ylabel": "rev", "title": "", "xscale": "linear",
             "yscale": "linear", "xlim": [0, 10], "ylim": [0, 100],
             "xlim_includes_zero": True, "ylim_includes_zero": True,
             "is_3d": False, "line_colors": ["#1f77b4"], "has_legend": False,
             "n_patches": 0, "n_images": 0},
            # Twin axis: matplotlib's ax.twinx() puts another axis with the
            # same xlim but different ylim and a different ylabel
            {"xlabel": "t", "ylabel": "clicks", "title": "", "xscale": "linear",
             "yscale": "linear", "xlim": [0, 10], "ylim": [0, 5000],
             "xlim_includes_zero": True, "ylim_includes_zero": True,
             "is_3d": False, "line_colors": ["#ff7f0e"], "has_legend": False,
             "n_patches": 0, "n_images": 0},
        ]}
        v = self._lint_spec(spec)
        self.assertIn("anti-patterns#dual-axes", _rules(v))

    def test_too_many_categorical_colors_major(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "linear",
                    "xlim": [0, 10], "ylim": [0, 100],
                    "xlim_includes_zero": True, "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": [
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                        "#bcbd22", "#17becf",  # 10 colors > 8 limit
                    ],
                    "has_legend": True,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("color#too-many-categorical", _rules(v))


class TestRedGreenPaletteIsCaught(unittest.TestCase):
    def test_red_green_pair_flagged(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "linear",
                    "xlim": [0, 10], "ylim": [0, 100],
                    "xlim_includes_zero": True, "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": ["#FF0000", "#00FF00"],
                    "has_legend": True,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = style_lint.lint_figure_spec(spec)
        self.assertIn("accessibility#colorblind-unsafe", _rules(v))
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v 2>&1 | tail -30
```

Expected: 10 failures / errors (new tests reference `lint_figure_spec` which doesn't exist; major describe rules not yet implemented).

- [ ] **Step 3: Add major describe-mode rules**

In `style_lint.py`, before the `# ----- Render output` section, add these regexes and append to `lint_describe`:

```python
_PIE_SLICES_PATTERN = re.compile(
    r"pie\s+chart\s+with\s+(\d+)\s+slices?|pie\s+chart.*?(\d+)\s+slices?",
    re.IGNORECASE,
)
_DYNAMITE_PATTERN = re.compile(
    r"bar\s+chart\s+of\s+means?\s+with\s+error\s+bars?",
    re.IGNORECASE,
)
_RED_GREEN_PALETTE_PATTERN = re.compile(
    r"colored?\s+red\s+and\s+green|red[\s/]green\s+(palette|colors?|encoding)",
    re.IGNORECASE,
)
_ALPHABET_ORDER_PATTERN = re.compile(
    r"sorted\s+alphabetic\w*|alphabetic\w*\s+order(ed)?",
    re.IGNORECASE,
)
```

And inside `lint_describe`, **before `return [v.to_dict() ...]`**, add:

```python
    # --- Major rules ---

    m = _PIE_SLICES_PATTERN.search(text)
    if m:
        n_str = next((g for g in m.groups() if g), None)
        if n_str and int(n_str) > 5:
            report.add(
                Violation(
                    severity="major",
                    rule="anti-patterns#pie-misuse",
                    where="chart_type",
                    observed=int(n_str),
                    fix="Use an ordered bar or dot plot. Pie charts require ≤5 slices.",
                )
            )

    if _DYNAMITE_PATTERN.search(text):
        report.add(
            Violation(
                severity="major",
                rule="anti-patterns#dynamite-plot",
                where="chart_type",
                fix=(
                    "Show individual observations — dot plot, box plot, or violin. "
                    "See Weissgerber 2015."
                ),
            )
        )

    if _RED_GREEN_PALETTE_PATTERN.search(text):
        report.add(
            Violation(
                severity="major",
                rule="accessibility#red-green-only",
                where="palette",
                fix=(
                    "Use Okabe-Ito (e.g., #E69F00 orange and #56B4E9 sky blue). "
                    "Or add a shape/line encoding so color isn't load-bearing."
                ),
            )
        )

    if _ALPHABET_ORDER_PATTERN.search(text):
        report.add(
            Violation(
                severity="major",
                rule="anti-patterns#alphabet-order",
                where="x_axis_order",
                fix=(
                    "Sort by value (ascending or descending), unless time or "
                    "experimental order is more meaningful."
                ),
            )
        )
```

- [ ] **Step 4: Add `lint_figure_spec` function**

In `style_lint.py`, add this function after `lint_describe`:

```python
def lint_figure_spec(spec: dict) -> list[dict]:
    """Run all rules against a JSON spec produced by figure_spec.extract_spec."""
    report = LintReport()
    axes = spec.get("axes", [])
    n_axes = spec.get("n_axes", len(axes))

    # --- Critical rules (per axis) ---
    for i, ax in enumerate(axes):
        where = f"axes[{i}]"

        # Labels required
        if not ax.get("xlabel") or not ax.get("ylabel"):
            report.add(
                Violation(
                    severity="critical",
                    rule="axes#labels-required",
                    where=where,
                    observed={
                        "xlabel": ax.get("xlabel", ""),
                        "ylabel": ax.get("ylabel", ""),
                    },
                    fix="Add labels with units in parentheses, e.g. 'Time (h)'.",
                )
            )

        # 3D
        if ax.get("is_3d"):
            report.add(
                Violation(
                    severity="critical",
                    rule="anti-patterns#3d-charts",
                    where=where,
                    fix="Use the 2D version. Always.",
                )
            )

        # Truncated baseline (only meaningful for bar/area charts — proxy: any patches)
        if (
            ax.get("yscale") == "linear"
            and not ax.get("ylim_includes_zero")
            and ax.get("n_patches", 0) > 0
        ):
            report.add(
                Violation(
                    severity="critical",
                    rule="anti-patterns#truncated-baseline",
                    where=f"{where}.ylim",
                    observed=ax.get("ylim"),
                    fix=(
                        "Extend axis to zero, OR annotate the truncation "
                        "with a broken-axis indicator."
                    ),
                )
            )

    # Dual y-axes: two axes that share xlim but differ on ylabel/ylim
    if n_axes == 2 and len(axes) == 2:
        a, b = axes[0], axes[1]
        if (
            a.get("xlim") == b.get("xlim")
            and a.get("ylabel") != b.get("ylabel")
            and (a.get("ylabel") and b.get("ylabel"))
            and a.get("ylim") != b.get("ylim")
        ):
            report.add(
                Violation(
                    severity="critical",
                    rule="anti-patterns#dual-axes",
                    where="axes[0,1]",
                    fix=(
                        "Use two adjacent panels with shared x-axis, OR "
                        "normalize both series to a common unit."
                    ),
                )
            )

    # --- Major rules (per axis) ---
    for i, ax in enumerate(axes):
        where = f"axes[{i}]"
        colors = [c for c in ax.get("line_colors", []) if c and c != "?"]

        # Too many categorical colors
        unique_colors = list(dict.fromkeys(colors))  # preserve order, dedupe
        if len(unique_colors) > 8:
            report.add(
                Violation(
                    severity="major",
                    rule="color#too-many-categorical",
                    where=f"{where}.line_colors",
                    observed=len(unique_colors),
                    fix=(
                        "≤8 categorical levels. Collapse categories or use "
                        "small multiples."
                    ),
                )
            )

        # Palette colorblind safety
        if len(unique_colors) >= 2:
            if not palettes.is_colorblind_safe_categorical(unique_colors):
                report.add(
                    Violation(
                        severity="major",
                        rule="accessibility#colorblind-unsafe",
                        where=f"{where}.line_colors",
                        observed=unique_colors,
                        fix=(
                            "Use Okabe-Ito categorical palette, or add a "
                            "redundant encoding (shape/line style)."
                        ),
                    )
                )

    return [v.to_dict() for v in report.violations]
```

- [ ] **Step 5: Wire `--figure-spec` mode in `cli_main`**

In `style_lint.py`, replace the `elif args.figure_spec:` block in `cli_main` with:

```python
    elif args.figure_spec:
        try:
            with open(args.figure_spec) as f:
                spec = json.load(f)
        except Exception as e:
            print(f"error: could not read figure-spec JSON: {e}", file=sys.stderr)
            return 1
        violations = lint_figure_spec(spec)
```

- [ ] **Step 6: Run tests to verify pass**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: 22 (prev) + 10 (this task) = 32 tests, all pass.

- [ ] **Step 7: Smoke-test figure-spec mode**

```bash
cat > /tmp/spec.json <<'EOF'
{
  "n_axes": 1,
  "axes": [
    {"xlabel": "", "ylabel": "", "title": "Revenue", "xscale": "linear",
     "yscale": "linear", "xlim": [0, 10], "ylim": [42, 51],
     "xlim_includes_zero": true, "ylim_includes_zero": false,
     "is_3d": false, "line_colors": [], "has_legend": false,
     "n_patches": 4, "n_images": 0}
  ]
}
EOF
python3 claude-config/skills/plotting-advisor/scripts/style_lint.py \
  --figure-spec /tmp/spec.json
```

Expected: markdown report with two critical violations (labels-required, truncated-baseline).

- [ ] **Step 8: Commit**

```bash
git add claude-config/skills/plotting-advisor/scripts/style_lint.py \
        claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): style_lint major rules + figure-spec mode

Major describe-mode rules (pie>5, dynamite plot, red/green-only,
alphabet order). Full figure-spec mode with rules for missing
labels, 3D, truncated baseline, dual axes, too-many-colors, and
colorblind-unsafe palette (delegated to palettes.is_colorblind
_safe_categorical).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 14: Add minor rules + remaining major rules

**Files:**
- Modify: `claude-config/skills/plotting-advisor/scripts/style_lint.py`
- Modify: `claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py`

- [ ] **Step 1: Add failing tests for minor and remaining rules**

Append to `test_style_lint.py`:

```python
class TestMinorRules(unittest.TestCase):
    def test_legend_with_few_series_flagged_describe(self):
        v = _run_describe(
            "line plot with 3 series using a legend; "
            "x-axis 'time (s)', y-axis 'voltage (V)'"
        )
        # Minor: legend when ≤5 series — suggest direct labels
        self.assertIn("annotation#legend-when-direct-labels-fit", _rules(v))
        for vv in v:
            if vv["rule"] == "annotation#legend-when-direct-labels-fit":
                self.assertEqual(vv["severity"], "minor")

    def test_legend_with_3_series_flagged_spec(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "linear",
                    "xlim": [0, 10], "ylim": [0, 100],
                    "xlim_includes_zero": True, "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": ["#E69F00", "#56B4E9", "#009E73"],
                    "has_legend": True,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = style_lint.lint_figure_spec(spec)
        self.assertIn("annotation#legend-when-direct-labels-fit", _rules(v))

    def test_log_scale_with_narrow_range_flagged(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "log",
                    "xlim": [0, 10], "ylim": [10, 30],  # less than 1 decade
                    "xlim_includes_zero": True, "ylim_includes_zero": False,
                    "is_3d": False,
                    "line_colors": ["#1f77b4"],
                    "has_legend": False,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = style_lint.lint_figure_spec(spec)
        self.assertIn("axes#log-narrow-range", _rules(v))


class TestStrictFlag(unittest.TestCase):
    def test_strict_exits_2_on_critical(self):
        # Use a separate process via cli_main
        rc = style_lint.cli_main([
            "--describe",
            "3D bar chart of revenue",
            "--strict",
        ])
        self.assertEqual(rc, 2)

    def test_no_strict_exits_0_on_critical(self):
        rc = style_lint.cli_main([
            "--describe",
            "3D bar chart of revenue",
        ])
        self.assertEqual(rc, 0)

    def test_strict_exits_0_on_no_critical(self):
        rc = style_lint.cli_main([
            "--describe",
            "dot plot of 4 conditions with Okabe-Ito palette, "
            "y axis 'Response (counts/min)' from 0 to 100, "
            "x axis 'Condition'",
            "--strict",
        ])
        self.assertEqual(rc, 0)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v 2>&1 | tail -20
```

Expected: 5 failures / errors.

- [ ] **Step 3: Add minor rules to `lint_describe` and `lint_figure_spec`**

In `style_lint.py`, add this regex with the other describe patterns:

```python
_LEGEND_FEW_SERIES_PATTERN = re.compile(
    r"(\d+)\s+series.*?legend|legend.*?(\d+)\s+series",
    re.IGNORECASE,
)
```

And append to the end of `lint_describe` (before the return):

```python
    # --- Minor rules ---

    m = _LEGEND_FEW_SERIES_PATTERN.search(text)
    if m:
        n_str = next((g for g in m.groups() if g), None)
        if n_str and 2 <= int(n_str) <= 5:
            report.add(
                Violation(
                    severity="minor",
                    rule="annotation#legend-when-direct-labels-fit",
                    where="annotation",
                    observed=int(n_str),
                    fix=(
                        "Direct-label each series at its terminal point "
                        "(Tufte). Legend is unnecessary at this count."
                    ),
                )
            )
```

In `lint_figure_spec`, add **after** the major-rules per-axis loop:

```python
    # --- Minor rules (per axis) ---
    for i, ax in enumerate(axes):
        where = f"axes[{i}]"
        colors = [c for c in ax.get("line_colors", []) if c and c != "?"]
        n_series = len(colors)

        # Legend when ≤5 series
        if ax.get("has_legend") and 2 <= n_series <= 5:
            report.add(
                Violation(
                    severity="minor",
                    rule="annotation#legend-when-direct-labels-fit",
                    where=f"{where}.legend",
                    observed=n_series,
                    fix="Direct-label each series at its terminal point (Tufte).",
                )
            )

        # Log scale with narrow range (<2 decades)
        if ax.get("yscale") == "log":
            ylo, yhi = ax.get("ylim", [None, None])
            if ylo and yhi and ylo > 0 and yhi > 0:
                import math
                decades = math.log10(yhi) - math.log10(ylo)
                if decades < 2.0:
                    report.add(
                        Violation(
                            severity="major",
                            rule="axes#log-narrow-range",
                            where=f"{where}.yaxis",
                            observed=round(decades, 2),
                            fix=(
                                "Log scale is appropriate only when range "
                                "spans ≥2 decades, or for multiplicative "
                                "processes. Use linear if not."
                            ),
                        )
                    )
        if ax.get("xscale") == "log":
            xlo, xhi = ax.get("xlim", [None, None])
            if xlo and xhi and xlo > 0 and xhi > 0:
                import math
                decades = math.log10(xhi) - math.log10(xlo)
                if decades < 2.0:
                    report.add(
                        Violation(
                            severity="major",
                            rule="axes#log-narrow-range",
                            where=f"{where}.xaxis",
                            observed=round(decades, 2),
                            fix=(
                                "Log scale is appropriate only when range "
                                "spans ≥2 decades, or for multiplicative "
                                "processes. Use linear if not."
                            ),
                        )
                    )
```

- [ ] **Step 4: Run tests to verify pass**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: 32 (prev) + 5 (this task) = 37 tests, all pass.

- [ ] **Step 5: Commit**

```bash
git add claude-config/skills/plotting-advisor/scripts/style_lint.py \
        claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): style_lint minor rules + log-range + --strict tests

Minor rules: legend-when-direct-labels-fit (both modes). Bumps
log-scale-with-narrow-range to major (axes#log-narrow-range).
Tests for --strict exit-code behavior (exits 2 on critical when
--strict; 0 otherwise; 0 when no critical).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 15: Add image mode (Pillow optional, fail-graceful)

**Files:**
- Modify: `claude-config/skills/plotting-advisor/scripts/style_lint.py`
- Modify: `claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py`

The image mode is the most heuristic — it cannot read axis labels reliably without OCR. It checks what it CAN see from raw pixels: dominant colors, color count, image dimensions. It fails gracefully when Pillow is unavailable.

- [ ] **Step 1: Add failing tests for image mode**

Append to `test_style_lint.py`:

```python
class TestImageMode(unittest.TestCase):
    def setUp(self):
        try:
            from PIL import Image  # noqa: F401
            self.have_pillow = True
        except ImportError:
            self.have_pillow = False

    def test_missing_image_returns_error_violation(self):
        violations = style_lint.lint_image("/nonexistent/path.png")
        self.assertTrue(any(v["severity"] == "error" for v in violations))

    def test_image_mode_skips_gracefully_without_pillow(self):
        if self.have_pillow:
            self.skipTest("Pillow is installed; cannot test no-Pillow path here")
        violations = style_lint.lint_image("/any/path.png")
        # Should return a single 'error' violation explaining Pillow is missing
        self.assertTrue(any(
            v["severity"] == "error" and "pillow" in v.get("fix", "").lower()
            for v in violations
        ))
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v 2>&1 | tail -10
```

Expected: 2 failures (lint_image doesn't exist).

- [ ] **Step 3: Add `lint_image` function**

In `style_lint.py`, after `lint_figure_spec`:

```python
def lint_image(path: str) -> list[dict]:
    """Inspect a saved PNG/PDF/JPG for rule violations.

    This is the coarsest mode. Without OCR, we cannot reliably read
    axis labels; we focus on what raw pixels reveal: dominant colors,
    rainbow-palette signature (continuous hue sweep), 3D perspective
    is not detectable from a raster, so we only catch palette issues.

    Returns 'error' severity violations on failure (file missing,
    Pillow unavailable, corrupt image) — never raises.
    """
    report = LintReport()

    try:
        from PIL import Image  # type: ignore
    except ImportError:
        report.add(
            Violation(
                severity="error",
                rule="lint#dependency-missing",
                where="image-mode",
                fix=(
                    "Install Pillow for image-mode lint: pip install Pillow. "
                    "Or use --figure-spec mode (no dependencies)."
                ),
            )
        )
        return [v.to_dict() for v in report.violations]

    import os
    if not os.path.exists(path):
        report.add(
            Violation(
                severity="error",
                rule="lint#file-not-found",
                where=path,
                fix="Check the path. Lint is advisory and never blocks.",
            )
        )
        return [v.to_dict() for v in report.violations]

    try:
        img = Image.open(path).convert("RGB")
    except Exception as e:
        report.add(
            Violation(
                severity="error",
                rule="lint#image-unreadable",
                where=path,
                observed=str(e),
                fix="Pillow could not open the file. Verify it's a valid PNG/JPG/etc.",
            )
        )
        return [v.to_dict() for v in report.violations]

    # Quantize to a small palette to catch dominant colors despite anti-aliasing.
    quantized = img.quantize(colors=16, method=Image.Quantize.FASTOCTREE)
    pal = quantized.getpalette() or []
    # Get color counts
    counts = quantized.getcolors(maxcolors=16) or []
    counts.sort(reverse=True)

    # Convert top non-background colors to hex
    # Heuristic: background is the most-frequent color; skip it.
    dominant_hexes: list[str] = []
    for i, (cnt, idx) in enumerate(counts):
        r, g, b = pal[idx * 3 : idx * 3 + 3]
        hex_code = "#{:02X}{:02X}{:02X}".format(r, g, b)
        # Skip near-white and near-black backgrounds for palette judgment
        is_near_white = r > 240 and g > 240 and b > 240
        is_near_black = r < 20 and g < 20 and b < 20
        if i == 0 and (is_near_white or is_near_black):
            continue
        dominant_hexes.append(hex_code)

    # Rule: rainbow signature — too many hues spread across the spectrum
    # is a soft indicator of a rainbow palette on continuous data.
    # We only flag when there are many distinct hues AND no clear primary.
    if len(dominant_hexes) >= 8:
        # Quick hue spread check
        hues = []
        for hex_code in dominant_hexes[:10]:
            r = int(hex_code[1:3], 16)
            g = int(hex_code[3:5], 16)
            b = int(hex_code[5:7], 16)
            mx = max(r, g, b)
            mn = min(r, g, b)
            if mx == mn:
                continue
            if mx == r:
                h = 60 * ((g - b) / (mx - mn)) % 360
            elif mx == g:
                h = 60 * ((b - r) / (mx - mn)) + 120
            else:
                h = 60 * ((r - g) / (mx - mn)) + 240
            hues.append(h % 360)
        if hues:
            hue_range = max(hues) - min(hues)
            if hue_range > 270:
                report.add(
                    Violation(
                        severity="critical",
                        rule="anti-patterns#rainbow-on-continuous",
                        where="image-palette",
                        observed={
                            "n_dominant_colors": len(dominant_hexes),
                            "hue_range_deg": round(hue_range, 1),
                        },
                        fix=(
                            "Dominant hues span >270° — consistent with a "
                            "rainbow palette. Use viridis / cividis (sequential) "
                            "or RdBu (diverging)."
                        ),
                    )
                )

    # Rule: colorblind-safety check on the top few dominant non-background colors
    if 2 <= len(dominant_hexes) <= 8:
        if not palettes.is_colorblind_safe_categorical(dominant_hexes[:8]):
            report.add(
                Violation(
                    severity="major",
                    rule="accessibility#colorblind-unsafe",
                    where="image-palette",
                    observed=dominant_hexes[:8],
                    fix=(
                        "Dominant palette is not colorblind-safe. Switch to "
                        "Okabe-Ito or Tol bright."
                    ),
                )
            )

    return [v.to_dict() for v in report.violations]
```

- [ ] **Step 4: Wire `--image` mode in `cli_main`**

In `style_lint.py`, replace the `elif args.image:` block in `cli_main` with:

```python
    elif args.image:
        violations = lint_image(args.image)
```

- [ ] **Step 5: Run tests to verify pass**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```

Expected: 37 (prev) + 2 (this task) = 39 tests, all pass.

- [ ] **Step 6: Smoke-test image mode with a missing file**

```bash
python3 claude-config/skills/plotting-advisor/scripts/style_lint.py \
  --image /tmp/does-not-exist.png
```

Expected: markdown report listing one `error`-severity violation (`lint#file-not-found` if Pillow installed, or `lint#dependency-missing` if not). Exit code 0.

- [ ] **Step 7: Commit**

```bash
git add claude-config/skills/plotting-advisor/scripts/style_lint.py \
        claude-config/skills/plotting-advisor/scripts/tests/test_style_lint.py
git commit -m "$(cat <<'EOF'
feat(plotting-advisor): style_lint image mode

Image mode uses Pillow (optional) to extract dominant colors after
16-color quantization, flags rainbow-palette signatures (hue span
>270°) and colorblind-unsafe palettes. Fails gracefully when
Pillow is unavailable or the file is missing/corrupt — emits
'error' severity violations rather than raising.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 16: Integration smoke + sync push + verify

**Files:**
- All previously created files; no new files

- [ ] **Step 1: Run the full test suite one more time**

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v 2>&1 | tail -10
```

Expected: 39 tests, all pass.

- [ ] **Step 2: Validate sync dry-run**

```bash
./sync-config.py push --dry-run 2>&1 | grep -E "plotting-advisor|ERROR|WARN" | head -30
```

Expected: shows `plotting-advisor` files being added to `~/.claude/skills/`; no ERROR or unexpected WARN lines.

- [ ] **Step 3: Run sync push**

```bash
./sync-config.py push
```

Expected: confirms push, no errors.

- [ ] **Step 4: Verify the skill is live at `~/.claude/`**

```bash
test -f ~/.claude/skills/plotting-advisor/SKILL.md && echo "SKILL.md present"
ls ~/.claude/skills/plotting-advisor/references/ | wc -l
ls ~/.claude/skills/plotting-advisor/scripts/ | wc -l
python3 -c "
import yaml
with open('$HOME/.claude/skills/plotting-advisor/SKILL.md') as f:
    content = f.read()
end = content.index('\n---\n', 4)
meta = yaml.safe_load(content[4:end])
assert meta['name'] == 'plotting-advisor', meta
print('YAML OK:', meta['name'], 'v' + str(meta.get('version', '?')))
"
```

Expected:
- `SKILL.md present`
- `7` (reference files count)
- `4` (scripts: palettes.py + figure_spec.py + style_lint.py + tests dir)
- `YAML OK: plotting-advisor v1.0`

- [ ] **Step 5: Integration smoke — describe mode**

```bash
python3 ~/.claude/skills/plotting-advisor/scripts/style_lint.py \
  --describe "3D bar chart with dual y-axes using jet colormap, no axis labels" \
  --json | python3 -m json.tool
```

Expected: JSON array containing critical violations: `anti-patterns#3d-charts`, `anti-patterns#dual-axes`, `anti-patterns#rainbow-on-continuous`, `axes#labels-required`.

- [ ] **Step 6: Integration smoke — clean description**

```bash
python3 ~/.claude/skills/plotting-advisor/scripts/style_lint.py \
  --describe "dot plot of 4 conditions using Okabe-Ito palette; x-axis 'Condition' (categorical); y-axis 'Response (counts/min)' from 0 to 100; direct-labeled groups"
```

Expected: `## Plotting Advisor: figure check — 0 issues found`.

- [ ] **Step 7: Verify `./sync-config.py status` is clean**

```bash
./sync-config.py status
```

Expected: no divergence between `claude-config/` and `~/.claude/`.

- [ ] **Step 8: Commit any incidental changes (typically nothing — sync is read-only on the repo)**

```bash
git status
```

If clean: proceed. If files changed (unlikely): inspect and commit separately.

---

## Task 17: Update planning entry with SHAs

**Files:**
- Modify: `planning/mac/2026-05-24-plotting-advisor-skill.md`

- [ ] **Step 1: Gather commit SHAs**

```bash
git log --oneline --since="2026-05-23" | grep -i "plotting-advisor\|plan plotting" | head -20
```

Expected: a list of ~15 commits (one per task, plus the original planning-entry commit).

- [ ] **Step 2: Update the planning entry**

Edit `planning/mac/2026-05-24-plotting-advisor-skill.md`:

- Set `**Status**: Success` (or `Partial` / `Failed` per actual outcome)
- Mark all `Changes Planned` checkboxes as done
- Fill in `Actual Outcome` with one paragraph summarizing what shipped
- Fill in `Assessment` (Result, Improvements, Issues, Lessons Learned)
- Replace the `[pending]` line under `Related Commits` with the actual SHAs from Step 1, one per line, with the commit message excerpt:

```markdown
## Related Commits

- abc1234: chore(planning): plan plotting-advisor skill
- def5678: feat(plotting-advisor): scaffold skill with SKILL.md
- ...
```

- [ ] **Step 3: Commit the planning update**

```bash
git add planning/mac/2026-05-24-plotting-advisor-skill.md
git commit -m "$(cat <<'EOF'
chore(planning): record SHAs on plotting-advisor entry

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 4: Final verification**

```bash
git log --oneline -20 | head -20
./sync-config.py status
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v 2>&1 | tail -5
```

Expected: all green. Skill is live, planning entry is complete, history is clean.

---

## Self-Review

### 1. Spec coverage

| Spec section | Task(s) covering it |
|---|---|
| Skill identity (frontmatter) | Task 2 |
| When to use / when not to use | Task 2 (in SKILL.md body) |
| Workflow — advisor flow | Task 2 |
| Workflow — lint flow | Task 2 + Tasks 12-15 |
| Inline rules (axis, annotation, accessibility floor, anti-pattern vetoes) | Task 2 |
| Output format — checklist + decision card | Task 2 (in SKILL.md) |
| Output format — lint | Tasks 12-15 (render_markdown in style_lint.py) |
| chart-selection.md | Task 4 |
| palettes.md | Task 5 |
| anti-patterns.md | Task 6 |
| accessibility.md | Task 7 |
| interactive-adaptation.md | Task 8 |
| scientific-conventions.md | Task 9 |
| SOURCES.md | Task 3 |
| palettes.py | Task 10 |
| figure_spec.py (import + CLI) | Task 11 |
| style_lint.py — describe mode | Task 12 |
| style_lint.py — figure-spec mode | Task 13 |
| style_lint.py — image mode | Task 15 |
| --strict flag | Task 14 (tests) + Task 12 (impl in cli_main) |
| Integration with software skills (delegate_to in card) | Task 2 (in SKILL.md output example) |
| Planning entry creation | Task 1 |
| Planning entry SHA update | Task 17 |
| sync-config.py push --dry-run + push | Task 16 |
| Smoke-test invocation | Task 16 |
| pytest smoke tests for each input mode + each severity-bearing rule | Tasks 10-15 (unittest, not pytest — stdlib choice) |

**Gaps**: None. The spec's "Build sequence" 8 steps map to Tasks 1, 2, 3-9, 10-11, 12-15, 16, 16, 17 respectively.

**Deviation from spec**: Spec said "smoke tests" generically. Implementation uses `unittest` (stdlib) rather than pytest to avoid adding a non-stdlib test dependency. Tests cover every input mode and every severity-bearing rule, per spec requirement.

### 2. Placeholder scan

- No `TBD`, `TODO`, `implement later`, or `fill in details`
- No "Add appropriate error handling" — every step shows the actual code
- No "Write tests for the above" without test code — every test block is concrete
- No "Similar to Task N" — each task is self-contained
- Reference docs use full content rather than templates with placeholders

### 3. Type / API consistency

- `Violation` dataclass: defined Task 12, used identically in Tasks 13-15
- `lint_describe(text) -> list[dict]`: Task 12, called consistently in tests
- `lint_figure_spec(spec) -> list[dict]`: Task 13, called consistently
- `lint_image(path) -> list[dict]`: Task 15, called consistently
- `palettes.is_colorblind_safe_categorical(hex_list)`: Task 10, called from Tasks 13 + 15 with identical signature
- `figure_spec.extract_spec(fig) -> dict`: Task 11, JSON shape consumed by Task 13 with matching field names (`xlabel`, `ylabel`, `xscale`, `yscale`, `xlim`, `ylim`, `xlim_includes_zero`, `ylim_includes_zero`, `is_3d`, `line_colors`, `has_legend`, `n_patches`, `n_images`, `title`)
- Rule identifiers (`anti-patterns#3d-charts`, etc.) match between SKILL.md anchors, reference doc anchors, and lint output

No drift between tasks.

---

## Plan complete and saved to `docs/superpowers/plans/2026-05-24-plotting-advisor.md`. Two execution options:

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints.

**Which approach?**
