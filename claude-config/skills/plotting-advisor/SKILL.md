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

```yaml
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
```
```

Every recommendation includes a rule citation (short tag like `[Cleveland 1985]` or a `references/<file>.md#anchor` pointer). Full bibliography lives in `references/sources.md`.

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
- `references/sources.md` — full bibliography

## Scripts

- `scripts/palettes.py` — canonical hex lists, `is_colorblind_safe_categorical(hex_list)` helper
- `scripts/figure_spec.py` — `extract_spec(fig)` for matplotlib Figure → JSON; CLI mode also accepts a pickled figure
- `scripts/style_lint.py` — main CLI with three input modes; `--strict` for CI use

Run tests with:

```bash
python3 -m unittest discover -s claude-config/skills/plotting-advisor/scripts/tests -v
```
