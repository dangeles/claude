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
