# Analysis Plan Template

## Background
[2-3 sentences providing context for the analysis]

Example:
> Cell surface receptors mediate cellular responses to diverse stimuli including neurotransmitters, hormones, and environmental signals. In model organisms, these receptors play critical roles in sensory perception, development, and homeostasis. Understanding receptor expression patterns across tissues provides insights into their physiological functions.

## Research Question
[One clear, specific question]

Example:
> How do receptor expression patterns differ between neuronal and non-neuronal tissues?

## Hypothesis
[Testable prediction based on biological knowledge]

Example:
> Sensory receptors will be enriched in neuronal tissues, while receptors mediating systemic processes (e.g., metabolism, reproduction) will be enriched in non-neuronal tissues.

## Methods

### Data
- [Describe dataset, source, sample size]
- [Specify organism, tissue types, experimental conditions]

Example:
> Single-cell RNA-seq data from published study: 40,000+ cells from developmental stage, annotated into 25+ cell types including neurons, muscle, epithelial, and connective tissues.

### Analysis Approach
[Bullet points describing computational steps]

Example:
- Filter genes of interest from expression matrix (n genes based on annotation)
- Calculate mean expression per cell type
- Perform differential expression analysis: tissue A vs. tissue B (Wilcoxon rank-sum test)
- Correct for multiple testing (Benjamini-Hochberg FDR < 0.05)
- Cluster genes by expression pattern (hierarchical clustering)
- Identify tissue-specific genes (fold-change > 2, adjusted p < 0.05)

### Statistical Tests
[Specify tests and significance thresholds]

Example:
- Differential expression: Wilcoxon rank-sum test (non-parametric, appropriate for count data)
- Multiple testing correction: Benjamini-Hochberg FDR < 0.05
- Effect size threshold: log2(fold-change) > 1

### Visualizations
[List planned figures]

Example:
1. Heatmap: Gene expression across cell types (top 50 variable genes)
2. Volcano plot: Differential expression (tissue A vs. tissue B)
3. Bar plot: Number of tissue-specific genes per tissue
4. Dot plot: Expression of candidate marker genes across cell types

## Expected Outcomes

### Predicted Results
[What you expect to find based on biological knowledge]

Example:
> Sensory receptor genes will show tissue-specific expression, particularly in specialized sensory cell types. Signaling receptors will be broadly expressed across cell types. Non-specialized tissues will express genes involved in systemic functions (e.g., metabolism, development).

### Alternative Scenarios
[What other results might mean]

Example:
> If many sensory genes are expressed in non-specialized tissues, this could indicate non-canonical roles or promiscuous expression patterns. Broad gene expression across cell types would suggest functional redundancy or pleiotropic roles.

## Quality Control Checks
[How to validate analysis quality]

Example:
- Positive controls: Known tissue-specific genes should show expected patterns
- Negative controls: Housekeeping genes should show uniform expression
- Biological validation: Results should match known expression patterns from literature
- Technical validation: Check for batch effects, outlier cells, doublets

## Success Criteria
[How to know if analysis answered the question]

Example:
- Clear separation of tissue-specific gene expression patterns
- Identification of tissue-specific marker genes
- Recapitulation of known biology (positive control validation)
- Novel insights into understudied gene families
