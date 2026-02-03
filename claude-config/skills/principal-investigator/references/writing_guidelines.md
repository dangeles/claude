# Scientific Writing Guidelines

Comprehensive guide for writing publication-quality scientific text.

## Journal-Specific Style Guides

### Nature Family (Nature, Nature Methods, Nature Communications)
- **Abstract**: 150-200 words, structured (Background, Methods, Results, Conclusions)
- **Introduction**: 3-4 paragraphs, emphasize novelty and significance
- **Results**: Present in logical order, not chronological
- **Methods**: Brief in main text, detailed in supplementary
- **Tense**: Past for specific findings, present for established knowledge
- **Citations**: Author (year) format, e.g., "Smith et al. (2020) showed..."
- **Figures**: Max 6-8 main figures, unlimited supplementary

### Cell Family (Cell, Cell Reports, iScience)
- **Abstract**: 150 words, unstructured single paragraph
- **Highlights**: 3-4 bullet points (85 characters each)
- **eTOC**: One sentence summary (40-50 words)
- **Results**: Combined Results and Discussion
- **Methods**: At end or supplementary
- **Tense**: Present tense for figures ("Figure 1A shows...")
- **Citations**: (Author, year) format, e.g., "(Smith et al., 2020)"

### Science Family (Science, Science Advances)
- **Abstract**: 125 words, structured (one-sentence summaries per section)
- **Main text**: ~3000 words total (very concise)
- **Subheadings**: Required in Results
- **Methods**: Supplementary materials only
- **Tense**: Present for data ("Data show...")
- **Citations**: Numbered (1, 2, 3)

### PLOS Family (PLOS ONE, PLOS Biology, PLOS Genetics)
- **Abstract**: 300 words, unstructured
- **Introduction**: No formal Abstract section, background integrated
- **Results**: Separate from Discussion
- **Methods**: Can be at end of main text
- **Figures**: Unlimited
- **Open Access**: CC-BY license, author retains copyright

### eLife
- **Abstract**: 150 words, single paragraph
- **eLife Digest**: 200-300 word lay summary (written by staff)
- **Structure**: Flexible, author decides organization
- **Methods**: Integrated in Results or separate
- **Transparency**: Data availability required

## Tense Usage

### Past Tense
Use for specific findings from your study:
> "We identified 327 differentially expressed genes"
> "Expression analysis revealed tissue-specific patterns"
> "RNA-seq data were collected from developmental stage X"

### Present Tense
Use for established facts and general truths:
> "Receptor proteins are essential for cell signaling"
> "The genome contains thousands of annotated genes"
> "Figure 2A shows gene expression patterns"

### Present Perfect
Use to connect past research to current state of knowledge:
> "Previous studies have shown that sensory receptors..."
> "Gene expression profiling has revealed..."

### Future Tense
Use sparingly, mainly in Discussion:
> "Future experiments will determine..."
> "This approach will enable..."

## Active vs. Passive Voice

### Modern Preference: Active Voice
Clearer, more direct, easier to read:
> ✅ "We analyzed 1,341 genes"
> ❌ "1,341 genes were analyzed"

> ✅ "We performed differential expression analysis"
> ❌ "Differential expression analysis was performed"

### When Passive is Appropriate
Focus on action/object rather than actor:
> ✅ "Samples were collected at L2 stage" (who collected is irrelevant)
> ✅ "Data were normalized using DESeq2" (method matters, not who)

## Common Phrases to Avoid

### Avoid Hedge Words (Unless Necessary)
| ❌ Weak | ✅ Strong |
|---------|-----------|
| "It appears that genes may be involved" | "Genes regulate signaling" |
| "We believe that our results suggest" | "Our results demonstrate" |
| "It is possible that GENE-X might" | "GENE-X regulates..." |

**When to hedge**: Speculative interpretations, unexpected findings
> ✅ "These results suggest GENE-X may regulate longevity" (hypothesis, not proven)

### Avoid Redundancy
| ❌ Redundant | ✅ Concise |
|--------------|------------|
| "We first began by analyzing" | "We analyzed" |
| "The results obtained showed" | "Results showed" |
| "A total of 327 genes" | "327 genes" |
| "Green in color" | "Green" |
| "Completely eliminate" | "Eliminate" |

### Avoid Vague Language
| ❌ Vague | ✅ Specific |
|----------|-------------|
| "Many genes" | "198 genes (60.6%)" |
| "Significantly different" | "Higher in tissue A (log2FC = 2.3, p < 0.001)" |
| "Recent studies" | "Smith et al. (2023)" |
| "Fairly consistent" | "Correlated (r = 0.87)" |

### Avoid Informal Language
| ❌ Informal | ✅ Formal |
|-------------|-----------|
| "Lots of genes" | "Numerous genes" or "198 genes" |
| "Turns out" | "Analysis revealed" |
| "Get" | "Obtain", "identify", "observe" |
| "Pretty significant" | "Significant (p < 0.001)" |

### Avoid Anthropomorphism
| ❌ Anthropomorphism | ✅ Appropriate |
|---------------------|-----------------|
| "Cells want to express genes" | "Cells express genes" |
| "Cells decide to differentiate" | "Cells differentiate in response to signals" |
| "Receptors try to bind ligands" | "Receptors bind ligands" |

## Paragraph Structure

### Standard Pattern
1. **Topic sentence**: Main point of paragraph
2. **Evidence**: Data, citations, examples
3. **Analysis**: Interpretation of evidence
4. **Transition**: Link to next paragraph

### Example
> Chemosensory genes are enriched in neuronal tissues (topic). Of 1,341 genes examined, 198 (60.6%) showed neuron-specific expression (log2FC > 1, FDR < 0.05) (evidence). This enrichment reflects the specialized role of neurons in detecting environmental cues through gene-mediated signaling (analysis). In contrast, non-neuronal tissues express genes involved in systemic functions (transition).

## Section-Specific Guidelines

### Abstract
**Structure (if journal requires)**:
1. Background (1-2 sentences): Context and knowledge gap
2. Methods (1-2 sentences): Approach used
3. Results (2-3 sentences): Key findings with numbers
4. Conclusions (1 sentence): Significance and impact

**Example**:
> Cell surface receptors mediate cellular responses to diverse stimuli, but tissue-specific expression patterns remain poorly characterized in many organisms. We analyzed single-cell RNA-seq data from 40,000+ cells to profile receptor genes across 27 cell types. We identified 327 genes with tissue-specific expression (FDR < 0.05), including 198 specialized-cell-enriched sensory receptors and 129 genes mediating systemic functions. These findings reveal receptor expression is highly cell-type-specific and provide a resource for functional studies.

### Introduction
**Funnel Structure**: Broad → Narrow → Study
1. **Paragraph 1**: General importance of topic
2. **Paragraph 2**: What is known, narrowing to specific area
3. **Paragraph 3**: Knowledge gap or unsolved problem
4. **Paragraph 4**: Your approach and objectives

**Avoid**:
- Overly broad openings ("Since the beginning of life...")
- Extensive methodology (save for Methods)
- Detailed results (save for Results)

### Results
**Organization**:
- Logical flow, not chronological order
- One main finding per section
- Subheadings help readability
- Refer to figures explicitly

**Example Structure**:
```
Gene expression is cell-type-specific
    ↓
Sensory genes are specialized-cell-enriched
    ↓
Gene families show distinct tissue patterns
    ↓
Co-expressed genes form functional modules
```

**Writing Pattern**:
1. State finding
2. Reference figure
3. Provide statistics
4. Interpret briefly

> We identified 327 genes with tissue-specific expression (Figure 2A). Of these, 198 (60.6%) were enriched in specialized cells (median log2FC = 2.3, FDR < 10⁻⁵), while 129 (39.4%) were enriched in non-specialized tissues. This enrichment suggests functional specialization of gene families across tissues.

### Discussion
**Structure**:
1. **Paragraph 1**: Restate main findings (no new data)
2. **Paragraphs 2-4**: Interpret findings in context of literature
3. **Paragraph 5**: Address limitations and caveats
4. **Paragraph 6**: Future directions
5. **Paragraph 7**: Concluding statement

**Key Principles**:
- Start with your findings, not others'
- Compare to published work (agreement/disagreement)
- Propose mechanisms
- Acknowledge limitations honestly
- End with impact/significance

### Methods
**Organization**:
- Chronological or by technique
- Subheadings for each method
- Enough detail for reproduction
- Cite established protocols

**Example Subheadings**:
- Data Acquisition
- Quality Control
- Differential Expression Analysis
- Clustering and Visualization
- Statistical Analysis

**Key Information to Include**:
- Software versions (DESeq2 v1.34.0)
- Statistical tests
- Significance thresholds
- Parameters used
- Data sources (accession numbers)

## Numbers and Statistics

### When to Use Numerals vs. Words

**Use numerals**:
- All measurements: "3 replicates", "5 mL", "10 µM"
- Statistics: "p = 0.003", "n = 42,035"
- Counts ≥10: "327 genes", "15 cell types"

**Use words**:
- Numbers <10 in narrative: "three major cell types"
- Start of sentence: "Twenty-seven cell types were analyzed" (or rephrase: "We analyzed 27 cell types")

### Reporting Statistics

**Always include**:
1. Test used
2. Test statistic
3. P-value
4. Sample size
5. Effect size

**Example**:
> Neurons expressed more genes than non-neurons (median 127 vs. 68, Wilcoxon W = 2.1×10⁸, p = 3.2×10⁻⁸, n = 22,418 and 19,617 cells, Cohen's d = 1.4).

**P-value Formatting**:
- p < 0.001 (not p = 0.000)
- p = 0.03 (not p < 0.05)
- p = 1.3×10⁻⁴⁵ (for very small values)

**Never**:
- "p = NS" (report actual value)
- "highly significant" (report p-value)
- p-value alone without effect size

## Abbreviations

### First Use
Define on first use in Abstract, main text, and each figure legend:
> "G protein-coupled receptors (genes)"

### Standard Abbreviations (No Definition Needed)
- DNA, RNA, ATP, GTP
- ANOVA, PCA, SEM, SD
- Fig., vs., e.g., i.e.

### Avoid Over-Abbreviating
If term appears <5 times, don't abbreviate:
> ❌ "single-cell RNA sequencing (scRNA-seq)" then used twice
> ✅ Just write "single-cell RNA sequencing" each time

## Inclusive Language

### Gender
- Avoid: "mailman", "chairman", "manpower"
- Use: "mail carrier", "chair", "workforce"

### Person-First Language
- Avoid: "diabetic patients"
- Use: "patients with diabetes"

### Age
- Avoid: "elderly subjects"
- Use: "older adults" or specific age range "adults aged 65-80"

## Citations

### How Many to Use
- **Introduction**: Cite key findings and reviews
- **Discussion**: Cite papers you compare to
- **Methods**: Cite original method descriptions

### What to Cite
- ✅ Published findings
- ✅ Software/tools
- ✅ Databases
- ❌ General knowledge ("DNA is double-stranded")
- ❌ Your own unpublished data (describe in methods)

### Citation Style Examples

**Nature (Author-year)**:
> Smith et al.¹ showed that genes regulate aging.

**Cell (Parenthetical)**:
> genes regulate aging in model organisms (Smith et al., 2020).

**Science (Numbered)**:
> genes regulate aging (1, 2).

## Revision Checklist

### Content
- [ ] Does abstract accurately reflect paper?
- [ ] Is introduction focused and concise?
- [ ] Are results presented logically?
- [ ] Does discussion interpret findings?
- [ ] Are conclusions supported by data?

### Clarity
- [ ] Is every sentence necessary?
- [ ] Are methods reproducible?
- [ ] Are figures referenced in text?
- [ ] Are all abbreviations defined?
- [ ] Are statistics complete?

### Style
- [ ] Consistent tense?
- [ ] Active voice preferred?
- [ ] No redundant phrases?
- [ ] Specific, not vague?
- [ ] Appropriate formality?

### Formatting
- [ ] Follows journal guidelines?
- [ ] Figures high quality?
- [ ] References formatted correctly?
- [ ] Supplementary materials complete?

## Common Errors and Fixes

| Error | Fix |
|-------|-----|
| "Data is" | "Data are" (plural) |
| "Between each sample" | "Among samples" (>2 items) |
| "Significant decrease" | "Significant decrease (p = 0.03)" |
| "As shown in Figure 1A" | "Figure 1A shows" or "(Figure 1A)" |
| "Higher as compared to controls" | "Higher than controls" |
| "Fold increase of 2 times" | "2-fold increase" or "increase of 2-fold" |
| "Correlated to" | "Correlated with" |
| "Comprised of" | "Composed of" or "comprises" |
