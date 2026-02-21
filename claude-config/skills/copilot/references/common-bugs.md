# Common Bioinformatics Bugs

Catalog of bugs frequently found in bioinformatics code, organized by category.

## Genomic Coordinates

### Off-By-One Errors (0-based vs 1-based)

**Problem**: Different file formats use different indexing systems.

| Format | Indexing | Example |
|--------|----------|---------|
| BED | 0-based, half-open [start, end) | chrA:0-100 = first 100 bp |
| GFF/GTF | 1-based, closed [start, end] | chrA:1-100 = first 100 bp |
| VCF | 1-based | POS=100 is 100th base |
| SAM/BAM | 1-based | POS=100 is 100th base |
| Python slicing | 0-based, half-open | seq[0:100] = first 100 chars |

**Bug Example**:
```python
# 🔴 Reading GFF (1-based) and using as Python slice (0-based)
gff = pd.read_csv('genes.gff', sep='\t')
gene_seq = genome_seq[gff['start'][0]:gff['end'][0]]  # Wrong! Off by 1
```

**Fix**:
```python
# ✅ Convert 1-based to 0-based for Python
start_0based = gff['start'][0] - 1  # GFF starts are 1-based
end_0based = gff['end'][0]          # GFF ends are 1-based inclusive
gene_seq = genome_seq[start_0based:end_0based]
```

### Chromosome Naming Inconsistency

**Problem**: Some files use "chr" prefix, others don't.

**Bug Example**:
```python
# 🔴 Coordinate mismatch
bam_reads = reads[reads['chr'] == 'chrA']  # BAM uses 'chrA'
gff_genes = genes[genes['seqid'] == 'A']   # GFF uses 'A'
# These won't overlap!
```

**Fix**:
```python
# ✅ Standardize chromosome names
def standardize_chr(chr_name):
    """Ensure chromosome names have 'chr' prefix."""
    if not str(chr_name).startswith('chr'):
        return f'chr{chr_name}'
    return chr_name

bam_reads['chr'] = bam_reads['chr'].apply(standardize_chr)
gff_genes['seqid'] = gff_genes['seqid'].apply(standardize_chr)
```

### Strand Confusion

**Problem**: Forgetting to reverse-complement for minus strand.

**Bug Example**:
```python
# 🔴 Extract sequence without considering strand
gene_seq = genome[start:end]
# If gene is on minus strand, this is wrong!
```

**Fix**:
```python
# ✅ Handle strand orientation
from Bio.Seq import Seq
gene_seq = genome[start:end]
if strand == '-':
    gene_seq = str(Seq(gene_seq).reverse_complement())
```

## Mathematical Errors

### Division by Zero

**Problem**: Normalizing or calculating ratios without checking denominator.

**Bug Examples**:
```python
# 🔴 Normalize counts
cpm = (counts / counts.sum()) * 1e6  # Fails if all counts are 0

# 🔴 Fold-change
fold_change = treatment / control  # Fails if control == 0

# 🔴 Coefficient of variation
cv = std / mean  # Fails if mean == 0
```

**Fixes**:
```python
# ✅ Check for zero before division
total = counts.sum()
cpm = (counts / total * 1e6) if total > 0 else np.zeros_like(counts)

# ✅ Add pseudocount
fold_change = (treatment + 1) / (control + 1)

# ✅ Handle zero mean
cv = std / mean if mean > 0 else np.nan
```

### Log of Zero

**Problem**: Taking logarithm of zero or negative values.

**Bug Examples**:
```python
# 🔴 Log transform counts
log_counts = np.log(counts)  # log(0) = -inf

# 🔴 Log2 fold-change
log2fc = np.log2(treatment / control)  # Two problems: 0/0 and log(0)
```

**Fixes**:
```python
# ✅ Use log1p (log(1 + x))
log_counts = np.log1p(counts)  # log1p(0) = 0

# ✅ Add pseudocount before log
log2fc = np.log2((treatment + 1) / (control + 1))
```

### Integer Overflow

**Problem**: Genomic positions exceed int32 range.

**Bug Example**:
```python
# 🔴 Position stored as int32
positions = np.array([123456789, 234567890], dtype=np.int32)
# Human chr1 is 248,956,422 bp - exceeds int32 max (2,147,483,647)
```

**Fix**:
```python
# ✅ Use int64 for genomic positions
positions = np.array([123456789, 234567890], dtype=np.int64)
```

### Floating Point Precision

**Problem**: Comparing floats with ==.

**Bug Example**:
```python
# 🔴 Float comparison
if p_value == 0.05:  # May never be exactly 0.05
    print("Significant")
```

**Fix**:
```python
# ✅ Use threshold
if p_value < 0.05:
    print("Significant")

# ✅ Or use np.isclose() for equality
if np.isclose(p_value, 0.05):
    print("Marginal")
```

## Statistical Errors

### No Multiple Testing Correction

**Problem**: Testing thousands of genes but using p < 0.05 directly.

**Bug Example**:
```python
# 🔴 No correction for multiple tests
sig_genes = genes[genes['p_value'] < 0.05]
# With 20,000 genes, expect ~1,000 false positives!
```

**Fix**:
```python
# ✅ Apply FDR correction
from statsmodels.stats.multitest import multipletests
_, genes['p_adj'], _, _ = multipletests(genes['p_value'], method='fdr_bh')
sig_genes = genes[genes['p_adj'] < 0.05]
```

### Wrong Test for Data Type

**Problem**: Using t-test on count data or non-parametric test unnecessarily.

**Bug Examples**:
```python
# 🔴 t-test on RNA-seq counts (violates normality assumption)
from scipy.stats import ttest_ind
t_stat, p_val = ttest_ind(treatment_counts, control_counts)

# 🔴 Non-parametric test on normally distributed data (loses power)
from scipy.stats import mannwhitneyu
u_stat, p_val = mannwhitneyu(normal_data_A, normal_data_B)
```

**Fixes**:
```python
# ✅ DESeq2 for RNA-seq counts
# Use pydeseq2 skill

# ✅ t-test for continuous normally distributed data
from scipy.stats import ttest_ind
t_stat, p_val = ttest_ind(continuous_A, continuous_B)
```

### Ignoring Effect Size

**Problem**: Reporting only p-values without magnitude of difference.

**Bug Example**:
```python
# 🔴 Only reporting significance
if p_value < 0.05:
    print(f"Gene X is significant (p = {p_value})")
# Might be p = 0.049 with log2FC = 0.01 (biologically meaningless)
```

**Fix**:
```python
# ✅ Report effect size + significance
if p_adj < 0.05 and abs(log2_fc) > 1:
    print(f"Gene X: log2FC = {log2_fc:.2f}, padj = {p_adj:.2e}")
```

### Paired vs Unpaired Tests

**Problem**: Using unpaired test on paired data or vice versa.

**Bug Example**:
```python
# 🔴 Unpaired test on paired samples (before/after treatment)
ttest_ind(before, after)  # Wrong! These are paired measurements
```

**Fix**:
```python
# ✅ Paired test
from scipy.stats import ttest_rel
ttest_rel(before, after)
```

## Data Handling Errors

### Index/Column Mismatch

**Problem**: Assuming count matrix columns match metadata rows.

**Bug Example**:
```python
# 🔴 Assume order matches
counts = pd.read_csv('counts.csv', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)
# Metadata might be in different order!
```

**Fix**:
```python
# ✅ Align explicitly
common_samples = counts.columns.intersection(metadata.index)
counts = counts[common_samples]
metadata = metadata.loc[common_samples]
assert (counts.columns == metadata.index).all()
```

### Missing Value Handling

**Problem**: Not handling NaN/missing values explicitly.

**Bug Example**:
```python
# 🔴 Correlation with NaN values
corr = np.corrcoef(gene_A, gene_B)  # Silently returns NaN if any value is NaN
```

**Fix**:
```python
# ✅ Explicit handling
# Option 1: Drop missing
valid = ~(np.isnan(gene_A) | np.isnan(gene_B))
corr = np.corrcoef(gene_A[valid], gene_B[valid])

# Option 2: Use pandas (handles NaN automatically)
corr = pd.Series(gene_A).corr(pd.Series(gene_B))
```

### String vs Numeric Confusion

**Problem**: Gene IDs read as numbers lose leading zeros.

**Bug Example**:
```python
# 🔴 Gene IDs like "0001" become 1
genes = pd.read_csv('genes.csv')  # Reads "0001" as integer 1
```

**Fix**:
```python
# ✅ Force string dtype
genes = pd.read_csv('genes.csv', dtype={'gene_id': str})
```

## Normalization Errors

### Wrong Normalization Order

**Problem**: Normalizing before filtering or vice versa.

**Bug Example**:
```python
# 🔴 Normalize then filter
cpm = (counts / counts.sum()) * 1e6
filtered = cpm[cpm.mean(axis=1) > 10]  # Wrong! Should filter counts, then normalize
```

**Fix**:
```python
# ✅ Filter then normalize
filtered_counts = counts[counts.mean(axis=1) > 10]
cpm = (filtered_counts / filtered_counts.sum()) * 1e6
```

### Batch-Specific Normalization

**Problem**: Normalizing batches separately instead of together.

**Bug Example**:
```python
# 🔴 Normalize each batch separately
batch1_norm = normalize(batch1)
batch2_norm = normalize(batch2)
combined = pd.concat([batch1_norm, batch2_norm])
# Batches now on different scales!
```

**Fix**:
```python
# ✅ Normalize together
combined = pd.concat([batch1, batch2], axis=1)
combined_norm = normalize(combined)
```

### Library Size vs Total Count Normalization

**Problem**: Using total counts when library size normalization is appropriate.

**Bug Example**:
```python
# 🔴 Simple total count normalization (doesn't account for composition bias)
normalized = counts / counts.sum(axis=0)
```

**Fix**:
```python
# ✅ Use proper normalization (DESeq2, TMM)
# For RNA-seq, use DESeq2 median-of-ratios or edgeR TMM
# See pydeseq2 skill for implementation
```

## Memory and Performance Errors

### Loading Entire File into Memory

**Problem**: Reading multi-GB files without chunking.

**Bug Example**:
```python
# 🔴 Load 100GB BAM file
reads = pd.read_csv('huge_file.bam')  # Out of memory!
```

**Fix**:
```python
# ✅ Process in chunks
for chunk in pd.read_csv('huge_file.csv', chunksize=10000):
    process(chunk)
```

### Unnecessary Copies

**Problem**: Creating DataFrame copies unnecessarily.

**Bug Example**:
```python
# 🔴 Unnecessary copy
df_copy = df.copy()
df_copy['new_col'] = df_copy['old_col'] * 2
```

**Fix**:
```python
# ✅ Modify in place (if appropriate)
df['new_col'] = df['old_col'] * 2
```

### Loop Instead of Vectorization

**Problem**: Iterating over DataFrame rows instead of using vectorized operations.

**Bug Example**:
```python
# 🔴 Slow loop
for i in range(len(df)):
    df.loc[i, 'log_expr'] = np.log1p(df.loc[i, 'expression'])
```

**Fix**:
```python
# ✅ Vectorized (100-1000x faster)
df['log_expr'] = np.log1p(df['expression'])
```

## Reproducibility Errors

### Missing Random Seed

**Problem**: Stochastic operations not reproducible.

**Bug Example**:
```python
# 🔴 No random seed
from sklearn.decomposition import PCA
pca = PCA(n_components=10)
transformed = pca.fit_transform(data)
# Different result each time if algorithm is stochastic
```

**Fix**:
```python
# ✅ Set random seed
import numpy as np
import random
np.random.seed(42)
random.seed(42)

from sklearn.decomposition import PCA
pca = PCA(n_components=10, random_state=42)
transformed = pca.fit_transform(data)
```

### Undefined Ordering

**Problem**: Relying on dictionary/set order (pre-Python 3.7).

**Bug Example**:
```python
# 🔴 Undefined order
gene_set = set(['GENE1', 'GENE2', 'GENE3'])
for gene in gene_set:  # Order may vary
    process(gene)
```

**Fix**:
```python
# ✅ Explicit ordering
gene_list = sorted(['GENE1', 'GENE2', 'GENE3'])
for gene in gene_list:
    process(gene)
```

### Hardcoded Paths

**Problem**: Absolute paths that won't work on other systems.

**Bug Example**:
```python
# 🔴 Hardcoded absolute path
data = pd.read_csv('/Users/alice/data/counts.csv')
```

**Fix**:
```python
# ✅ Relative paths or Path objects
from pathlib import Path
data_dir = Path('data')
data = pd.read_csv(data_dir / 'counts.csv')
```

## Visualization Errors

### Misleading Axes

**Problem**: Y-axis doesn't start at zero for bar plots.

**Bug Example**:
```python
# 🔴 Bar plot with truncated Y-axis
plt.bar(['A', 'B'], [98, 100])
plt.ylim(95, 101)  # Makes 2% difference look huge
```

**Fix**:
```python
# ✅ Start at zero for bar plots
plt.bar(['A', 'B'], [98, 100])
plt.ylim(0, 110)  # Honest representation
```

### Overplotting

**Problem**: Thousands of points plotted as opaque dots.

**Bug Example**:
```python
# 🔴 10,000 opaque points
plt.scatter(x, y)  # Can't see density
```

**Fix**:
```python
# ✅ Use transparency or density plot
plt.scatter(x, y, alpha=0.1)
# Or hexbin
plt.hexbin(x, y, gridsize=50)
```

### Missing Labels

**Problem**: No axis labels, legend, or title.

**Bug Example**:
```python
# 🔴 Unlabeled plot
plt.scatter(x, y)
plt.show()
```

**Fix**:
```python
# ✅ Complete labeling
plt.scatter(x, y, label='Sample A')
plt.xlabel('Gene Expression (log2 CPM)')
plt.ylabel('Protein Abundance (log2)')
plt.title('Gene-Protein Correlation')
plt.legend()
plt.show()
```

## Version Compatibility Errors

### Breaking API Changes

**Problem**: Code written for old package version fails with new version.

**Bug Example**:
```python
# 🔴 Code written for pandas 0.25
df.append(new_row)  # Removed in pandas 2.0
```

**Fix**:
```python
# ✅ Check package versions and use compatible methods
import pandas as pd
print(f"pandas version: {pd.__version__}")

# For pandas >= 2.0
df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
```

### Python 2 vs 3

**Problem**: Code written for Python 2 fails in Python 3.

**Bug Examples**:
```python
# 🔴 Python 2 division
result = 5 / 2  # Python 2: 2, Python 3: 2.5

# 🔴 Python 2 print statement
print "Hello"  # SyntaxError in Python 3
```

**Fix**:
```python
# ✅ Explicit division
result = 5 // 2  # Integer division: 2 in both versions
result = 5 / 2   # Float division: 2.5 in Python 3

# ✅ print() function
print("Hello")  # Works in both versions
```

## Summary Checklist

When reviewing code, check for:
- [ ] Coordinate system clearly documented (0-based vs 1-based)
- [ ] Chromosome names consistent across files
- [ ] Strand orientation handled for sequences
- [ ] Division by zero prevented
- [ ] Log of zero/negative prevented
- [ ] Multiple testing correction applied
- [ ] Appropriate statistical test for data type
- [ ] Effect size reported with p-value
- [ ] Sample order matches between counts and metadata
- [ ] Missing values handled explicitly
- [ ] Normalization order correct (filter then normalize)
- [ ] Random seed set for reproducibility
- [ ] No hardcoded absolute paths
- [ ] Memory-efficient for large data
- [ ] Plots properly labeled
- [ ] Package versions documented
