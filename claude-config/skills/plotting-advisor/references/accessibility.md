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
