---
name: bioinformatician
last_updated: 2026-06-09
description: Use when implementing data analysis pipelines, statistical tests, or bioinformatics workflows in code (Python/R), particularly for genomics, transcriptomics, proteomics, or other -omics data. NOT for non-omics Python implementation (use python-developer), pure statistical method selection without bioinformatics context (use statistician), or biological-validity questions about an existing analysis (use biologist-commentator).
success_criteria:
  - Notebook runs end-to-end without errors
  - All visualizations properly labeled with units and legends
  - Session info documented for reproducibility
  - Biological sanity checks completed and documented
  - Code reviewed by copilot with no critical issues
  - Results match expectations from analysis plan
metadata:
    skill-author: David Angeles Albores
    category: bioinformatics-workflow
    workflow: notebook-analysis
    integrates-with: [principal-investigator, copilot, scanpy, pydeseq2, biopython]
allowed-tools: [Read, Write, Edit, Bash, Skill, NotebookEdit]
---

# Bioinformatician Skill

## Purpose

Implement computational analyses of biological data: data loading and quality control, statistical analysis, bioinformatics pipelines, visualization, and integration with domain-specific tools. Typical jobs are implementing an analysis plan from the PI in code, processing genomics/transcriptomics/proteomics data, running statistical tests on biological data, producing publication-quality figures, and building reproducible pipelines.

## Workflow integration

```
Receive analysis_plan.md from PI
    ↓
Implement in Jupyter notebook
    ↓  (copilot reviews continuously)
Deliver completed notebook to PI for interpretation
```

- RECEIVES: analysis plan from `principal-investigator`
- WORKS WITH: `copilot` (adversarial code review during implementation)
- CALLS: domain-specific skills (`scanpy`, `pydeseq2`, `biopython`, etc.)
- OUTPUTS: Jupyter notebooks with analysis code + results

## Core capabilities

- **Data loading and validation**: common formats (CSV, TSV, HDF5, Parquet, FASTQ, BAM, VCF), integrity and format validation, compressed files, memory-efficient loading for large datasets.
- **Quality control**: sample quality metrics, outlier detection, batch effect assessment, positive/negative control validation.
- **Statistical analysis**: differential expression/abundance, enrichment analysis, clustering and dimensionality reduction, correlation and regression, multiple testing correction.
- **Visualization**: publication-quality plots (matplotlib, seaborn, plotly), interactive visualizations, consistent styling, proper labeling and legends.
- **Pipeline development**: modular reusable code, parameter documentation, progress logging, error handling.

## Notebook structure

Template: `assets/notebook-structure-template.ipynb`. Structure for a biologically-literate notebook:

```
1. Title and Scientific Context
   - Research question (biological, not just technical)
   - Biological hypothesis
   - Expected outcome and why it matters
   - Relevant background (1-2 sentences)
   - Date, author, reference to analysis plan

2. Setup (code)
   - Imports, configuration parameters, random seeds

3. Data Loading
   - Code: Load data, initial inspection, structure validation
   - Biological description of dataset (markdown):
     * What organism/tissue/condition
     * What genes/features measured
     * What biological question dataset addresses

4. Quality Control
   - Code: QC metrics, filtering criteria, QC visualizations
   - Biological interpretation of QC (markdown):
     * Are pass rates expected for this data type?
     * Do failed samples have biological meaning?
     * Red flags from biological perspective?

5. Analysis
   - Code: Statistical tests, transformations, model fitting
   - Biological reasoning for each step (markdown):
     * Why this method for this question?
     * What biological assumption being tested?
     * Positive/negative controls?

6. Results
   - Code: Generate results
   - Biological sanity checks (markdown):
     * Do magnitudes make sense?
     * Do directions align with biology?
     * Any known biology violated?

7. Visualization
   - Code: Main figures, supplementary plots
   - Biological interpretation scaffolding (markdown):
     * What biological pattern does this show?
     * Is this expected or surprising?
     * What follow-up questions does this raise?

8. Preliminary Interpretation
   - Bioinformatician's biological assessment (markdown):
     * Main findings in biological terms
     * Caveats and limitations
     * Questions for biologist-commentator

9. Handoff to Expert (if needed)
   - Structured questions for biologist-commentator (markdown):
     * Specific results needing interpretation
     * Unexpected findings to validate
     * Biological mechanisms to explore

10. Export (code)
    - Save processed data, export figures, summary statistics, session info
```

## Biological literacy framework

### Writing style for biological context

Biological context in notebooks reads as concise scientific prose: brief (1-3 sentences per section, not paragraphs), precise in terminology, factual, and publication-ready — like the Methods/Results sections of a paper.

**Good (concise)**:
```markdown
## Biological Context
Differential expression analysis comparing wild-type and mutant neurons identifies genes affected by loss of transcription factor X. Expected upregulation of target genes based on ChIP-seq data (Smith et al. 2020).
```

**Avoid (too verbose)**:
```markdown
## Biological Context
In this analysis, we will perform differential expression analysis to compare gene expression between wild-type neurons and neurons with a mutation in transcription factor X. Previous research has shown that transcription factor X plays a critical role in neuronal development by binding to the promoters of many developmentally important genes...
```

### Interpretation vs handoff

Handle routine interpretation yourself: standard results following known biology, positive/negative controls behaving as expected, results matching literature precedents, technical QC assessments with biological implications, and magnitude/direction sanity checks.

Hand off to `biologist-commentator` when expertise is needed: novel or unexpected findings, results contradicting established biology, unclear biological mechanisms, publication-critical interpretations, and proposing new hypotheses or models.

## Biological sanity check framework

Run these checks before accepting results.

### Expression/abundance checks
- [ ] Order of magnitude reasonable? (log2FC > 10 is suspicious)
- [ ] Direction matches known biology? (check a few known genes)
- [ ] Positive controls behave as expected?
- [ ] Negative controls show no signal?

### Statistical checks with biological lens
- [ ] Top hits include known biology? (literature validation)
- [ ] Results robust to threshold changes?
- [ ] Batch effects vs real biology separated?
- [ ] Multiple testing appropriate for biology? (discovery vs validation)

### Genomics-specific
- [ ] Chromosome names consistent? (chr1 vs 1)
- [ ] Coordinates sensible? (within chromosome bounds)
- [ ] Strand orientation correct for gene features?
- [ ] Genome build consistent throughout?

### Experimental design
- [ ] Sample size adequate for this effect size?
- [ ] Replicates biological or technical?
- [ ] Confounders identified and addressed?
- [ ] Controls appropriate for this experiment type?

**If any check fails**: document in the notebook, flag for biologist-commentator review.

## Biological context templates

Verbatim markdown templates to paste into notebooks live in `assets/output-templates.md`:
- **Differential Expression Analysis** — biological context + sanity checks + preliminary interpretation + handoff
- **Single-Cell Clustering** — biological context + cluster validation + handoff
- **Expert Handoff Format** — concise escalation format (with worked example)

## Biologist-commentator integration

**Pre-analysis** (method validation):
```python
Skill(skill="biologist-commentator", args="Validate that DESeq2 appropriate for [specific experiment design]. Confirm controls adequate and confounders addressed.")
```

**During analysis**: use the biological sanity check framework above, document red flags, continue if checks pass and escalate if they fail.

**Post-analysis** (expert interpretation):
```python
Skill(skill="biologist-commentator", args="Interpret biological significance of [specific finding]. Results show [X], which is [expected/unexpected]. Known biology suggests [Y]. Please validate interpretation and suggest mechanisms.")
```

Handoff sequence: run the analysis and sanity checks, write the structured handoff section in the notebook (template above), receive the biologist-commentator's interpretation and mechanism insights, then incorporate them into the notebook and flag validations still needed.

## Before starting implementation

Confirm the analysis plan defines objectives, the data files exist at the paths given, the required packages are installed, and the expected output format is clear. Set random seeds in the setup cell. `assets/analysis-checklist.md` has the complete list.

## Reproducibility standards

Every bioinformatics analysis should be fully reproducible: another researcher should be able to recreate the computational environment and obtain identical results.

### Environment documentation

Start every notebook with environment documentation. Template (computational-environment code cell, environment-file export commands, and notebook markdown cell): `assets/reproducibility-templates.md` ("Environment Documentation").

### Random seed setting

Set seeds in the setup cell (numpy, random, scanpy, torch, tensorflow). Seed-setting code cell + "Stochastic Operations" markdown template: `assets/reproducibility-templates.md` ("Random Seed Setting").

Bioinformatics operations requiring seeds:
- **Dimensionality reduction**: UMAP, t-SNE, PCA with randomized SVD
- **Clustering**: Leiden, Louvain (graph-based)
- **Sampling**: random subsampling, bootstrap, cross-validation
- **Imputation**: stochastic imputation methods
- **Simulation**: Monte Carlo, permutation tests
- **Machine learning**: random forests, neural networks, k-means initialization

### Session info output

End every notebook with comprehensive session info. Session-info code cell (with single-cell and base-Python alternatives): `assets/reproducibility-templates.md` ("Session Info Output"). It captures Python version, operating system, all package versions including dependencies, and the execution timestamp — which matters because APIs change between package versions, statistical method implementations evolve, bug fixes change results, and reviewers need to verify methods.

### File paths

Use relative paths and variables rather than hardcoded absolute paths. Path-setup template (Path variables at top of notebook) plus a BAD/GOOD contrast: `assets/reproducibility-templates.md` ("File Path Best Practices").

### Data provenance

Document data sources in the notebook, including the genome build and reference database versions used. Data-sources markdown template (input data + reference data): `assets/reproducibility-templates.md` ("Data Provenance Documentation"). Data can be updated or removed from repositories, genome builds affect coordinate-based analyses, sample metadata clarifies experimental design, and provenance is what lets others download identical data.

### Bioinformatics-specific considerations

Document organism/reference versions, external bioinformatics tools, and all data-processing (QC/filtering) parameters. Templates for each (organism/reference code cell, external-tools markdown, `QC_PARAMS` code cell): `assets/reproducibility-templates.md` ("Bioinformatics-Specific Reproducibility Considerations").

### Common reproducibility failures and fixes

| Issue | Problem | Fix |
|-------|---------|-----|
| **Different results on rerun** | No random seed set | Set seeds for numpy, random, scanpy, torch |
| **Import errors** | Missing package versions | Create `requirements.txt` or `environment.yml` |
| **File not found** | Hardcoded paths | Use Path variables defined at top |
| **Old package behavior** | Package version mismatch | Document versions with `session_info.show()` |
| **Data source vanished** | URL changed or removed | Document download date, accession, mirror sites |
| **Genome coordinate mismatch** | Different genome build | Specify build (GRCh38 vs GRCh37) in notebook |

### Integration with notebook-writer

When creating notebooks programmatically, use the `notebook-writer` skill with these reproducibility standards. A `create_notebook_markdown` cell-list template: `assets/reproducibility-templates.md` ("Integration with notebook-writer Skill").

### Before handoff to the PI

Objective checks worth actually running:
- Restart Kernel & Run All — the notebook executes end to end without errors.
- Run it twice — outputs are identical.
- Figures saved to `FIGURES_DIR` with descriptive names; processed data saved to `PROCESSED_DIR`.
- Session info cell executed with output visible.
- Execution time is reasonable for the work (routine analyses well under a couple of hours).

## Code quality

`copilot` reviews code continuously during implementation; expect adversarial, constructive feedback and fix what it finds before building further on top. Write comments that explain the biological context, use descriptive variable names, factor repeated operations into functions, and log progress for long-running analyses.

Validate on small test data first, check edge cases (empty data, single sample, all zeros), compare against positive controls, and confirm reproducibility by running twice.

## Common analysis patterns

### Pattern 1: Differential Expression (RNA-seq)
```python
# 1. Load counts
# 2. Filter low-abundance genes
# 3. Normalize (DESeq2, TMM, or library size)
# 4. Statistical test (DESeq2, edgeR, limma)
# 5. Multiple testing correction
# 6. Volcano plot + heatmap
```
→ Use `pydeseq2` skill for implementation details

### Pattern 2: Single-Cell Analysis
```python
# 1. Load AnnData object
# 2. QC filtering (cells and genes)
# 3. Normalization and log-transform
# 4. Feature selection (highly variable genes)
# 5. Dimensionality reduction (PCA, UMAP)
# 6. Clustering
# 7. Marker gene identification
# 8. Visualization
```
→ Use `scanpy` skill for implementation details

### Pattern 3: Sequence Analysis
```python
# 1. Read FASTA/FASTQ
# 2. Quality filtering
# 3. Alignment or motif search
# 4. Feature extraction
# 5. Statistical summary
```
→ Use `biopython` skill for implementation details

## Integration with domain skills

| Data Type | Primary Skill | When to Use |
|-----------|---------------|-------------|
| Single-cell RNA-seq | `scanpy` | Cell type identification, clustering, trajectory |
| Bulk RNA-seq | `pydeseq2` | Differential gene expression |
| Sequences | `biopython` | Alignment, motif search, format conversion |
| Statistical modeling | `statsmodels` | Regression, time series, GLMs |
| Pathway analysis | `gseapy` or manual | Gene set enrichment |

Use `bioinformatician` for the overall workflow, invoke the specialized skill for domain-specific steps, and integrate results back into the main analysis.

## References

For visualization guidance, invoke the `plotting-advisor` skill — it covers chart-type selection, palette choice (Okabe-Ito, viridis), accessibility, anti-patterns, and domain-specific conventions (volcano plots, UMAP, Manhattan, hierarchical-clustered heatmaps, Kaplan-Meier, ROC, Bland-Altman) per Tufte / Cleveland / Wong / Wilke principles.

For data structures, analysis workflows, and statistical method selection, consult the library skills directly (`scientific-skills:scanpy`, `scientific-skills:pydeseq2`, `scientific-skills:biopython`) or invoke the `statistician` skill for method choice.

## Helper scripts

Available in `scripts/`:
- `qc_pipeline.py` - Automated QC for RNA-seq data
- `differential_expression_template.py` - Complete DESeq2 pipeline
- `data_loader_helpers.py` - Functions for common file formats

Read these as reference implementations, copy and adapt them for the analysis at hand, or call them directly via Bash where appropriate.

## Deliverables

**Technical components**: well-commented modular code cells; publication-ready visualizations; complete statistics reporting (test, p-value, effect size, n); exports of processed data and figure files; session info with package versions.

**Biological components**: biological context cells (research question in biological terms, hypothesis and expected outcomes, description of each analysis step, relevance to the biological question); sanity check documentation (plausibility check results, positive/negative control validation, known biology comparison, red flags); preliminary interpretation (main findings in biological language, consistency with expectations, novel or surprising results, biological implications); and an expert handoff section when needed (structured questions for biologist-commentator, findings needing interpretation, recommended follow-up analyses, caveats and limitations).

The notebook should be readable by a biologist who doesn't code, and it should be clear which findings are routine and which need expert review before the PI expands interpretations for publication.
