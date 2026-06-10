# Output Templates: Biological Context & Expert Handoff

Verbatim markdown templates the bioinformatician emits into notebooks. Copy and fill in the bracketed placeholders. See `SKILL.md` for when to use each (sanity-check framework, handoff workflow).

## Template: Differential Expression Analysis

```markdown
## Biological Context
Comparing [condition A] vs [condition B] to identify genes involved in [biological process]. Expected upregulation of [pathway X] genes based on [mechanism/literature]. Positive controls: [gene1, gene2]. Expected log2FC range: [X-Y] based on [citation].

## Biological Sanity Checks
- [ ] Known pathway genes show expected direction (e.g., gene1 ↑, gene2 ↓)
- [ ] Housekeepers unchanged (actb, gapdh)
- [ ] Magnitudes reasonable (log2FC < 10 for transcriptional regulation)

## Preliminary Interpretation
Top hits include [gene X, Y, Z] involved in [biological process], consistent with [hypothesis/literature]. [Gene W] unexpected - requires expert validation.

**Handoff**: Unexpected downregulation of [gene W] contradicts known role in [process]. Biologist-commentator needed for mechanism assessment.
```

## Template: Single-Cell Clustering

```markdown
## Biological Context
Clustering [tissue] cells to identify cell types. Expected populations: [celltype1 (markers: a,b,c), celltype2 (markers: d,e,f)]. Reference atlas: [citation if available].

## Cluster Validation
- Cluster 1: [celltype] - markers: [genes] ✓
- Cluster 2: [celltype] - markers: [genes] ✓
- Cluster 3: Novel population - markers: [genes] - needs expert review

**Handoff**: Cluster 3 shows unexpected marker combination [X+Y+Z-]. Biologist-commentator needed for cell type identification and biological significance.
```

## Template: Expert Handoff Format

Use this concise format when escalating to biologist-commentator:

```markdown
## Expert Interpretation Needed

**Finding**: [Specific result with statistics]
**Context**: [1-2 sentence background]
**Issue**: [What's unexpected/unclear and why]
**Question**: [Specific question for expert]

**Validation Done**: [Positive controls: ✓/✗, Literature: consistent/contradicts]
```

**Example**:
```markdown
## Expert Interpretation Needed

**Finding**: Gene X shows 8-fold upregulation (padj<0.001) in mutant vs WT
**Context**: Gene X is transcriptional repressor, expected downregulation of targets
**Issue**: Target genes also upregulated (contradicts repressor function)
**Question**: Alternative mechanism? Post-transcriptional regulation? Data artifact?

**Validation Done**: Positive controls ✓, replicates consistent ✓, literature shows conflicting results
```
