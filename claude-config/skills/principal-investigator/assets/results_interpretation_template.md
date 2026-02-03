# Results Interpretation Template

Use this structure to interpret analysis results in notebook markdown cells or manuscript results sections.

## Pattern: Observation → Evidence → Meaning → Caveats

### 1. State the Observation
[What did you find? Objective, factual description]

Example:
> We identified 327 genes differentially expressed between tissue A and tissue B (FDR < 0.05, |log2FC| > 1).

### 2. Provide Statistical Evidence
[Quantify the finding with statistics]

Example:
> Of these, 198 genes (60.6%) showed enrichment in tissue A (median log2FC = 2.3), while 129 genes (39.4%) showed enrichment in tissue B (median log2FC = -1.8). The largest effect sizes were observed for sensory receptor genes: GENE-1 (log2FC = 5.2, FDR = 1.3×10⁻⁴⁵) and GENE-2 (log2FC = 4.8, FDR = 2.1×10⁻³⁸).

### 3. Explain Biological Meaning
[Why does this matter? What does it tell us biologically?]

Example:
> This enrichment of sensory receptors in tissue A aligns with their known roles in detecting environmental cues. The high expression of GENE-1 and GENE-2 in specialized cell types suggests functional specialization. Conversely, tissue B enrichment of metabolic genes reflects their roles in systemic processes.

### 4. Acknowledge Caveats
[What are the limitations? What could explain alternative interpretations?]

Example:
> These findings are based on L2-stage larvae and may not reflect expression patterns in other developmental stages. Additionally, low-abundance genes may be underestimated due to sequencing depth limitations. Cell type annotations are computational predictions and require experimental validation.

## Examples for Common Result Types

### Differential Expression Results

**Template**:
> [N] genes were differentially expressed between [condition A] and [condition B] ([test], [threshold]). [X%] were upregulated in [condition A] (median [effect size]), while [Y%] were upregulated in [condition B] (median [effect size]). Notable genes include [specific examples with statistics]. These results suggest [biological interpretation], consistent with [known biology or previous studies]. However, [caveats about technical/biological limitations].

### Clustering Results

**Template**:
> Unsupervised clustering identified [N] distinct clusters ([method], [parameters]). Cluster [X] contained [N genes/cells] characterized by [defining features]. Gene ontology enrichment analysis revealed [pathway/process] (FDR = [value]). This pattern suggests [biological interpretation]. The separation of [cluster A] from [cluster B] indicates [meaning], though [caveat about interpretation].

### Correlation Results

**Template**:
> [Gene A] and [Gene B] showed [strength] correlation (r = [value], p = [p-value]). This association is consistent with [biological relationship], as [explanation]. However, correlation does not imply causation, and [alternative explanations].

### Positive Control Validation

**Template**:
> As a positive control, we examined [known gene/pattern]. As expected, [gene X] showed [expected pattern] ([statistics]), confirming the validity of our analysis. This recapitulation of known biology increases confidence in novel findings.

### Unexpected Findings

**Template**:
> Unexpectedly, we observed [surprising result] ([statistics]). This finding contrasts with [previous expectation], potentially because [hypothetical explanation]. Alternative explanations include [technical artifact possibility] or [novel biological mechanism]. Further experiments are needed to [validation approach].

## Figure Interpretation Integration

When referring to figures:

> **Figure 1A shows [what the figure displays]**. [Describe key patterns]. [Interpret biological meaning]. See figure legend for complete methods and statistics.

Example:
> **Figure 2C shows gene expression across 27 cell types**. Sensory receptor genes (blue cluster) are highly enriched in specialized sensory neurons, while signaling receptors (green cluster) are broadly expressed across neuronal subtypes. This tissue-specific expression pattern reflects functional specialization of gene families. See Figure 2 legend for clustering methods and statistical tests.

## Writing for Different Audiences

### For Lab Notebooks
- Emphasize observations and raw findings
- Include negative results and troubleshooting notes
- Document decisions and reasoning

### For Manuscripts
- Lead with biological significance
- Emphasize novelty and impact
- Minimize technical details (move to methods)
- Connect to broader context

### For Grant Proposals
- Highlight preliminary success
- Emphasize feasibility
- Connect to proposed aims
- Show expertise

## Common Mistakes to Avoid

### ❌ Over-interpretation
> "gene-123 **causes** increased lifespan"
(correlation ≠ causation)

### ✅ Appropriate Interpretation
> "gene-123 expression **correlates with** increased lifespan, suggesting a potential role in longevity regulation"

### ❌ Hiding Negative Results
> "We identified 327 differentially expressed genes"
(omitting that 1,014 showed no difference)

### ✅ Complete Reporting
> "Of 1,341 genes examined, 327 (24.4%) were differentially expressed, while 1,014 (75.6%) showed similar expression across tissues"

### ❌ Unsupported Claims
> "This proves genes regulate aging"
(single correlation doesn't prove mechanism)

### ✅ Evidence-Based Claims
> "This association suggests genes may contribute to aging regulation, consistent with previous genetic studies (Smith et al., 2020)"

### ❌ Ignoring Effect Size
> "gene-5 was significantly different (p < 0.05)"
(p-value without magnitude)

### ✅ Reporting Effect and Significance
> "gene-5 showed modest enrichment in neurons (log2FC = 0.8, FDR = 0.03), suggesting weak tissue preference"

## Integration with Analysis Plan

**Always refer back to your analysis plan**:
- Did results match expectations?
- Were hypotheses supported or rejected?
- What new questions emerged?

Example:
> Our hypothesis that chemosensory genes would be neuron-enriched was strongly supported (198/327 neuron-enriched genes, 60.6%). However, we unexpectedly found 42 genes with mixed expression across neuronal and non-neuronal tissues, suggesting more complex regulation than anticipated. This motivates future investigation of cell-type-specific alternative splicing.
