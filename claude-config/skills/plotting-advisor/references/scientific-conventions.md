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
