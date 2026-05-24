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
