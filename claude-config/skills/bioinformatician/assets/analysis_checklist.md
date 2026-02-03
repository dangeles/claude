# Analysis Pre-Flight Checklist

Complete this checklist before starting any analysis implementation.

## 1. Requirements Review

- [ ] **Analysis plan received** from PI or clearly defined
- [ ] **Research question** understood
- [ ] **Expected outputs** defined (figures, tables, statistics)
- [ ] **Success criteria** clear

## 2. Data Validation

- [ ] **Data files exist** at specified paths
- [ ] **File formats correct** (CSV, HDF5, FASTQ, etc.)
- [ ] **File sizes reasonable** (not corrupted/truncated)
- [ ] **Data structure matches expectations** (columns, samples, genes)
- [ ] **Metadata available** (sample annotations, conditions)
- [ ] **Positive controls present** (if applicable)
- [ ] **Negative controls present** (if applicable)

## 3. Environment Setup

- [ ] **Required packages installed** (check environment.yml or requirements.txt)
- [ ] **Package versions compatible** (no known conflicts)
- [ ] **Random seeds set** for reproducibility
- [ ] **Working directory correct**
- [ ] **Output directories created**

## 4. Analysis Design

- [ ] **Statistical test appropriate** for data type (consult `references/statistical_methods.md`)
- [ ] **Sample size adequate** for test power
- [ ] **Normalization method selected** (if needed)
- [ ] **Multiple testing correction planned** (Bonferroni, BH, etc.)
- [ ] **Thresholds defined** (p-value, fold-change, FDR)
- [ ] **Edge cases considered** (empty data, missing values, zeros)

## 5. Reproducibility

- [ ] **Random seed(s) set** (numpy, random, torch if applicable)
- [ ] **Package versions documented** (will add session info at end)
- [ ] **Parameters explicit** (not hardcoded, defined at top)
- [ ] **File paths relative or configurable** (not absolute hardcoded paths)

## 6. Quality Control Strategy

- [ ] **QC metrics defined** (what to check)
- [ ] **Filtering criteria established** (thresholds for exclusion)
- [ ] **Outlier detection method** chosen
- [ ] **Batch effect assessment** planned (if multi-batch data)
- [ ] **Technical replicates** handled appropriately

## 7. Visualization Plan

- [ ] **Figure types chosen** appropriate for data
- [ ] **Color schemes accessible** (colorblind-friendly if publication)
- [ ] **Figure dimensions** suitable for publication
- [ ] **File formats** defined (PNG for web, PDF for publication)

## 8. Integration Points

- [ ] **Specialized skills identified** (scanpy, pydeseq2, biopython)
- [ ] **Helper scripts available** (in `scripts/` if needed)
- [ ] **Copilot review** expected during implementation

## 9. Before Execution

- [ ] **Test on subset first** (small sample to verify workflow)
- [ ] **Estimate compute time** (avoid surprises)
- [ ] **Memory requirements** considered (for large datasets)
- [ ] **Error handling** planned (what to do if step fails)

## 10. After Execution

- [ ] **All code executed successfully**
- [ ] **Results make biological sense** (sanity check)
- [ ] **Positive controls validated** (expected results reproduced)
- [ ] **Figures properly labeled** (axes, legends, titles)
- [ ] **Statistics complete** (test, statistic, p-value, n, effect size)
- [ ] **Data exported** (processed files saved)
- [ ] **Figures exported** (high-resolution files saved)
- [ ] **Session info added** (package versions)
- [ ] **Notebook ready for PI interpretation**

## Common Gotchas to Check

### Data Loading
- [ ] Correct delimiter (comma vs tab vs space)
- [ ] Header row present/absent
- [ ] Index column specified
- [ ] Data types correct (numeric not read as strings)
- [ ] Missing values handled (NaN, NA, blank)

### Genomics-Specific
- [ ] **Coordinate system** (0-based vs 1-based indexing)
- [ ] **Chromosome naming** (chr1 vs 1)
- [ ] **Strand orientation** (+ vs -)
- [ ] **Genome version** (hg38, mm10, WS289, etc.)

### Statistics
- [ ] **Assumptions checked** (normality, equal variance, independence)
- [ ] **One-tailed vs two-tailed** test appropriate
- [ ] **Paired vs unpaired** test correct
- [ ] **Multiple testing correction** applied
- [ ] **Effect size** reported (not just p-value)

### Reproducibility
- [ ] **Set random seeds** before any stochastic operation
  ```python
  import numpy as np
  import random
  np.random.seed(42)
  random.seed(42)
  ```
- [ ] **Sort before operations** that may have undefined order
- [ ] **Version dependencies** documented

### Performance
- [ ] **Memory efficient** for large data (chunking, lazy loading)
- [ ] **Parallelization** used appropriately (if beneficial)
- [ ] **Intermediate results saved** (for long pipelines)
- [ ] **Progress logging** for long-running steps

## Decision Trees

### Which normalization method?

**RNA-seq (bulk)**:
- DESeq2 median-of-ratios → Best for differential expression
- TMM (edgeR) → Alternative, similar performance
- TPM/FPKM → For visualization, comparison across genes

**RNA-seq (single-cell)**:
- Library size normalization → Basic
- SCTransform → Handles technical variation
- scran pooling → Accounts for composition bias

**Microarray**:
- Quantile normalization → Standard
- RMA → For Affymetrix

**Proteomics**:
- Total intensity normalization → Basic
- Median normalization → Robust to outliers
- Quantile normalization → Strong assumption

### Which test for differential expression?

**Data characteristics**:
- Count data (RNA-seq) → DESeq2, edgeR, limma-voom
- Continuous (microarray, proteomics) → limma, t-test
- Single-cell → Wilcoxon, MAST, DESeq2 (pseudobulk)

**Sample size**:
- n < 3 per group → Avoid statistics, descriptive only
- n = 3-5 → DESeq2 (shrinkage helps), Wilcoxon
- n > 5 → Any appropriate test

**Replicates**:
- Biological replicates → Standard tests
- Technical replicates → Average first, then test
- No replicates → Can't test differential expression (descriptive only)

### Which clustering method?

**Data type**:
- Gene expression → Hierarchical or k-means
- Single-cell → Louvain, Leiden (graph-based)
- Spatial → Graph-based with spatial constraint

**Number of clusters**:
- Known a priori → k-means with k specified
- Unknown → Hierarchical (dendogram), or Louvain (resolution parameter)

**Distance metric**:
- Euclidean → General purpose
- Correlation → Gene expression patterns
- Cosine → High-dimensional sparse data

## Troubleshooting

### Issue: Code runs but results look wrong

**Check**:
1. Positive controls (known biology reproduced?)
2. Scale of effects (fold-changes biologically plausible?)
3. Data loaded correctly (plot raw data)
4. Normalization applied (before and after)
5. Filtering too aggressive (lost important features?)

### Issue: Statistical test returns no significant results

**Check**:
1. Sample size adequate? (power analysis)
2. Effect size large enough? (biological vs statistical significance)
3. Multiple testing correction too stringent? (try FDR instead of Bonferroni)
4. Variance too high? (outliers, batch effects)
5. Test assumptions met? (normality, equal variance)

### Issue: Memory error with large dataset

**Solutions**:
1. Load in chunks (pandas `chunksize` parameter)
2. Use Dask for lazy evaluation
3. Downsample for exploratory analysis
4. Switch to HDF5 or Parquet format
5. Filter early (remove low-abundance features)

### Issue: Results not reproducible

**Check**:
1. Random seed set before all stochastic operations
2. Data not sorted (undefined order in hash tables, sets)
3. Package versions different
4. Parallel processing (non-deterministic order)

## Final Validation

Before delivering notebook to PI:

- [ ] **Run notebook top to bottom** (Restart kernel & run all)
- [ ] **No errors or warnings** (or explicitly handled)
- [ ] **Figures render correctly**
- [ ] **Numbers match text** (no copy-paste errors)
- [ ] **Copilot approved** (no outstanding critical issues)
- [ ] **Ready for biological interpretation** by PI
