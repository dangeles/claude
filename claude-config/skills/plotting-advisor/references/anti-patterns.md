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
