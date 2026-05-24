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
