# Figure Legend Examples

Publication-quality figure legends for common bioinformatics visualizations.

## General Structure

```
Figure [N]. [One-sentence description of main finding].
([A-Z]) [Panel-specific description]. [Methods summary]. [Statistical test]. [Sample size]. [Key observations].
```

## Example 1: Heatmap (Gene Expression)

**Figure 1. Gene expression varies across cell types.**
(A) Heatmap showing scaled expression (z-score) of top 50 variable genes across 27 cell types from single-cell RNA-seq (n = 40,000 cells). Rows: Genes (hierarchical clustering, Euclidean distance). Columns: Cell types (grouped by tissue). Color scale: blue (low expression) to red (high expression). Sensory receptor genes (blue sidebar) cluster with specialized cell types. (B) Number of genes detected per cell type (≥ 10 counts). Box plots show median (center line), quartiles (box), and 1.5×IQR (whiskers). Specialized cells express significantly more receptors than non-specialized cells (median 127 vs. 68, Wilcoxon p = 3.2×10⁻⁸).

## Example 2: Volcano Plot (Differential Expression)

**Figure 2. Receptor genes are enriched in specialized tissues.**
Volcano plot showing differential expression of 1,300+ genes between specialized cells (n = 22,000 cells) and non-specialized cells (n = 19,000 cells). X-axis: log2(fold-change). Y-axis: -log10(FDR-adjusted p-value). Horizontal dashed line: FDR = 0.05 threshold. Vertical dashed lines: |log2FC| = 1 threshold. Red points: specialized-enriched genes (n = 198). Blue points: non-specialized-enriched genes (n = 129). Gray points: not significant. Statistical test: Wilcoxon rank-sum test with Benjamini-Hochberg correction. Labeled genes: top 10 by effect size. Notable specialized-enriched genes include GENE-1 (log2FC = 5.2), GENE-2 (log2FC = 4.8), and GENE-3 (log2FC = 4.3).

## Example 3: PCA Plot (Dimensionality Reduction)

**Figure 3. Cell types cluster by tissue of origin.**
Principal component analysis (PCA) of 42,035 cells based on expression of 1,341 genes. Each point represents a cell, colored by annotated cell type. PC1 (23.4% variance) separates neurons (left) from non-neurons (right). PC2 (14.7% variance) separates muscle/intestine (bottom) from hypodermis/glia (top). Ellipses: 95% confidence intervals per cell type. PCA performed on log-normalized counts. This separation indicates gene expression profiles are cell-type-specific and reflect developmental lineage.

## Example 4: Bar Plot with Error Bars

**Figure 4. Receptor gene families show tissue-specific expression.**
Number of genes per family significantly enriched in specialized tissue (FDR < 0.05, log2FC > 1). Bars: mean ± SEM across 5 specialized cell types. Family A (n = 23 genes) shows highest enrichment, followed by Family B (n = 18) and Family C (n = 12). Family D (n = 6) shows modest enrichment. Statistical test: One-way ANOVA with Tukey post-hoc (F(4,60) = 18.3, p = 1.2×10⁻⁹). ***: p < 0.001, **: p < 0.01, *: p < 0.05 vs. Family D.

## Example 5: Scatter Plot with Regression

**Figure 5. Gene count correlates with cell type diversity.**
Scatter plot showing correlation between number of genes expressed per cell type (x-axis) and number of unique cell states within that type (y-axis). Each point: one of 27 cell types. Line: linear regression (R² = 0.64, p = 1.3×10⁻⁶). Shaded region: 95% confidence interval. Specialized cells (red points) express more genes (median 127) than non-specialized cells (blue points, median 68) and show greater within-type heterogeneity. This correlation suggests these genes contribute to cell state specification. Outliers (labeled): certain cells express few genes (52) despite high heterogeneity (8 states).

## Example 6: Box Plot (Distribution Comparison)

**Figure 6. Olfactory genes show bimodal expression in sensory neurons.**
Box plots showing expression distribution (log2(CPM+1)) of three gene classes across specialized sensory neurons (cell types A-E; n = 1,247 cells). Sensory receptor genes (red, n = 89) show bimodal distribution with distinct ON/OFF states (Hartigan's dip test, p = 0.002). Metabolic genes (blue, n = 52) show unimodal low expression (p = 0.34). Signaling receptors (green, n = 71) show uniform low expression. Box plots: center line (median), box (IQR), whiskers (1.5×IQR), points (outliers). Statistical comparison: Kruskal-Wallis test (H(2) = 127.3, p = 1.1×10⁻²⁸) with Dunn's post-hoc correction.

## Example 7: Venn Diagram (Set Overlap)

**Figure 7. Cell-type-specific genes show minimal overlap.**
Venn diagram showing overlap of genes enriched in three major tissues: neurons (red, n = 198 genes), muscle (blue, n = 43), and intestine (green, n = 67). Enrichment criteria: log2FC > 1, FDR < 0.05 vs. all other cell types. Neurons and intestine share 12 genes (5.2%), likely reflecting neuroendocrine signaling. Only 3 genes are enriched in all three tissues (gene-47, gene-118, gene-203), suggesting broad signaling roles. 183 genes (92.4%) are neuron-specific, indicating high tissue specialization.

## Example 8: Network Diagram (Gene Relationships)

**Figure 8. Co-expressed genes form functional modules.**
Network showing co-expression relationships among top 100 variable genes. Nodes: genes (size = mean expression, color = dominant cell type). Edges: Pearson correlation > 0.7 (p < 0.01). Node clustering: Louvain algorithm identified 5 modules (outlined). Module 1 (red): sensory receptor genes (n = 27), expressed in specialized neurons. Module 2 (blue): signaling receptors (n = 18), broadly neuronal. Module 3 (green): metabolic receptors (n = 12), epithelial tissues. Modules 4-5: small clusters (n = 8, 5). Hub genes (degree > 10): GENE-1, GENE-2, GENE-3. This modular structure suggests coordinated regulation of functionally related genes.

## Example 9: Time Course (Line Plot)

**Figure 9. gene expression dynamics during development.**
Line plots showing mean expression (log2(CPM+1)) of four gene families across five developmental stages (embryo, L1, L2, L3, adult). Lines: mean ± SEM (n = 3 biological replicates per stage). Chemosensory genes (red) increase from L1 to adult (linear model slope = 0.34, p = 0.003), coinciding with sensory neuron maturation. Neuropeptide receptors (blue) remain stable (slope = 0.02, p = 0.72). Metabolic genes (green) peak at L3 (one-way ANOVA, F(4,10) = 12.4, p = 0.001), aligning with rapid growth phase. Germline genes (purple) are adult-specific (undetectable before L4). Shaded regions: 95% confidence intervals.

## Example 10: Dot Plot (Expression + Prevalence)

**Figure 10. Cell-type marker genes show restricted expression.**
Dot plot showing expression of 15 candidate marker genes across 27 cell types. Dot size: percentage of cells expressing gene (detection threshold: > 0 counts). Dot color: mean expression in positive cells (log2(CPM+1), blue = low, red = high). Each row: one gene. Each column: one cell type. Candidate markers selected by: (1) high cell type specificity (entropy < 1.5), (2) high expression in target cell type (log2FC > 2), (3) prevalence in target type (> 50% cells). GENE-1 is cell-type-A-specific (98% type A cells, <2% other cells). GENE-2 is broadly neuronal. GENE-3 is epithelial-specific. These markers enable cell type identification in uncharacterized samples.

## Key Elements of Good Figure Legends

### 1. Start with Main Finding
> ❌ "Heatmap of gene expression"
> ✅ "gene expression varies across cell types"

### 2. Define All Visual Elements
- Axes (what, units, scale)
- Colors (what they represent, scale bar)
- Symbols (shapes, sizes, meanings)
- Lines (what they connect, statistical fits)

### 3. Report Sample Sizes
> "n = 42,035 cells"
> "3 biological replicates per condition"
> "27 cell types"

### 4. State Statistical Tests
> "Wilcoxon rank-sum test with Benjamini-Hochberg correction"
> "Linear regression (R² = 0.64, p = 1.3×10⁻⁶)"

### 5. Describe Key Methods
> "PCA performed on log-normalized counts"
> "Hierarchical clustering, Euclidean distance"

### 6. Highlight Key Observations
> "Neurons express significantly more genes than non-neuronal cells (median 127 vs. 68)"

### 7. Define Thresholds
> "Detection threshold: > 10 counts"
> "FDR < 0.05, |log2FC| > 1"

### 8. Explain Abbreviations
> First use: "G protein-coupled receptors (genes)"
> Subsequent uses: "genes"

### 9. Note Data Transformations
> "Scaled expression (z-score)"
> "log2(counts per million + 1)"

### 10. Cite Panel Labels
> "(A) Heatmap... (B) Box plot..."

## Common Mistakes to Avoid

### ❌ Too Vague
> "Figure shows gene expression is different"

### ✅ Specific
> "198 genes are neuron-enriched (log2FC > 1, FDR < 0.05)"

### ❌ Missing Methods
> "PCA of cells"

### ✅ Complete Methods
> "PCA of 42,035 cells based on expression of 1,341 genes. PCA performed on log-normalized counts."

### ❌ No Statistics
> "Expression was higher in neurons"

### ✅ Quantified
> "Expression was higher in neurons (median 127 vs. 68 genes, Wilcoxon p = 3.2×10⁻⁸)"

### ❌ Undefined Visual Elements
> "Red points show significant genes"

### ✅ Defined
> "Red points: neuron-enriched genes (n = 198, log2FC > 1, FDR < 0.05)"
