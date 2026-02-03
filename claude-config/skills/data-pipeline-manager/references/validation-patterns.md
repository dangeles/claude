# Validation Patterns Reference

This document provides comprehensive validation patterns for common data pipeline scenarios, organized by validation type and data format.

## Table of Contents

1. [File System Validation](#file-system-validation)
2. [Format Validation](#format-validation)
3. [Schema Validation](#schema-validation)
4. [Data Quality Validation](#data-quality-validation)
5. [Bioinformatics-Specific Validation](#bioinformatics-specific-validation)
6. [Cross-Stage Validation](#cross-stage-validation)

## File System Validation

### File Existence and Accessibility

```python
import os
import stat

def validate_file_access(file_path, require_read=True, require_write=False):
    """Comprehensive file access validation"""

    # Check existence
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    # Check it's a file (not directory)
    if not os.path.isfile(file_path):
        raise ValueError(f"Path is not a file: {file_path}")

    # Check read permission
    if require_read and not os.access(file_path, os.R_OK):
        raise PermissionError(f"File not readable: {file_path}")

    # Check write permission
    if require_write and not os.access(file_path, os.W_OK):
        raise PermissionError(f"File not writable: {file_path}")

    return True
```

### Directory Structure Validation

```python
def validate_directory_structure(base_dir, required_structure):
    """
    Validate directory structure matches expected layout

    Args:
        base_dir: Base directory path
        required_structure: Dict mapping paths to expected file counts

    Example:
        required = {
            'raw_data': {'min_files': 1, 'type': 'directory'},
            'raw_data/*.fastq.gz': {'min_files': 2, 'type': 'file'},
            'metadata.csv': {'type': 'file'}
        }
    """
    import glob

    errors = []

    for path_pattern, requirements in required_structure.items():
        full_pattern = os.path.join(base_dir, path_pattern)

        if '*' in path_pattern:
            # Pattern matching
            matches = glob.glob(full_pattern)
            if len(matches) < requirements.get('min_files', 1):
                errors.append(f"Expected at least {requirements['min_files']} files "
                            f"matching {path_pattern}, found {len(matches)}")
        else:
            # Exact path
            full_path = os.path.join(base_dir, path_pattern)
            if not os.path.exists(full_path):
                errors.append(f"Missing required path: {path_pattern}")
            elif requirements['type'] == 'directory' and not os.path.isdir(full_path):
                errors.append(f"Expected directory: {path_pattern}")
            elif requirements['type'] == 'file' and not os.path.isfile(full_path):
                errors.append(f"Expected file: {path_pattern}")

    if errors:
        raise ValueError(f"Directory structure validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return True
```

### Disk Space Validation

```python
import shutil

def validate_disk_space(path, required_gb, buffer_factor=1.2):
    """
    Check sufficient disk space is available

    Args:
        path: Path to check
        required_gb: Required space in GB
        buffer_factor: Safety factor (1.2 = 20% buffer)
    """
    usage = shutil.disk_usage(path)
    available_gb = usage.free / (1024**3)
    required_with_buffer = required_gb * buffer_factor

    if available_gb < required_with_buffer:
        raise ValueError(
            f"Insufficient disk space at {path}:\n"
            f"  Available: {available_gb:.1f} GB\n"
            f"  Required: {required_with_buffer:.1f} GB (including {buffer_factor}x buffer)"
        )

    return available_gb
```

## Format Validation

### FASTQ Validation

```python
import gzip

def validate_fastq_format(file_path, check_n_records=1000):
    """
    Validate FASTQ format structure

    Checks:
    - Proper 4-line structure
    - Header lines start with @
    - Plus lines start with +
    - Quality and sequence length match
    """

    opener = gzip.open if file_path.endswith('.gz') else open
    errors = []

    with opener(file_path, 'rt') as f:
        record_num = 0

        while record_num < check_n_records:
            # Read one record (4 lines)
            try:
                header = f.readline()
                if not header:  # EOF
                    break

                seq = f.readline()
                plus = f.readline()
                qual = f.readline()

                line_num = record_num * 4 + 1

                # Check header
                if not header.startswith('@'):
                    errors.append(f"Line {line_num}: Header must start with '@', got: {header[:20]}")

                # Check plus line
                if not plus.startswith('+'):
                    errors.append(f"Line {line_num + 2}: Plus line must start with '+', got: {plus[:20]}")

                # Check sequence and quality length match
                seq_len = len(seq.strip())
                qual_len = len(qual.strip())
                if seq_len != qual_len:
                    errors.append(f"Line {line_num}: Sequence length ({seq_len}) != "
                                f"quality length ({qual_len})")

                # Check sequence contains valid bases
                valid_bases = set('ACGTNacgtn')
                seq_bases = set(seq.strip())
                invalid = seq_bases - valid_bases
                if invalid:
                    errors.append(f"Line {line_num + 1}: Invalid bases: {invalid}")

                record_num += 1

            except Exception as e:
                errors.append(f"Error reading record {record_num}: {e}")
                break

    if record_num == 0:
        raise ValueError(f"No valid FASTQ records found in {file_path}")

    if errors:
        raise ValueError(f"FASTQ validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors[:10]))  # Limit error output

    return record_num
```

### CSV Validation

```python
import pandas as pd

def validate_csv_format(file_path, expected_columns=None, required_columns=None,
                       column_types=None, max_rows_check=None):
    """
    Comprehensive CSV validation

    Args:
        file_path: Path to CSV file
        expected_columns: Expected column names (order matters)
        required_columns: Columns that must be present (order doesn't matter)
        column_types: Dict mapping column names to expected dtypes
        max_rows_check: Only check first N rows for performance
    """

    # Read CSV
    try:
        df = pd.read_csv(file_path, nrows=max_rows_check)
    except Exception as e:
        raise ValueError(f"Cannot read CSV file {file_path}: {e}")

    errors = []

    # Check if empty
    if len(df) == 0:
        raise ValueError(f"CSV file is empty: {file_path}")

    # Check expected columns (exact match, order matters)
    if expected_columns is not None:
        if list(df.columns) != expected_columns:
            errors.append(f"Column mismatch.\n"
                         f"  Expected: {expected_columns}\n"
                         f"  Found: {list(df.columns)}")

    # Check required columns (subset, order doesn't matter)
    if required_columns is not None:
        missing = set(required_columns) - set(df.columns)
        if missing:
            errors.append(f"Missing required columns: {missing}")

    # Check column types
    if column_types is not None:
        for col, expected_type in column_types.items():
            if col not in df.columns:
                continue  # Will be caught by required_columns check

            actual_type = df[col].dtype
            if not pd.api.types.is_dtype_equal(actual_type, expected_type):
                # Try to convert
                try:
                    df[col] = df[col].astype(expected_type)
                except Exception:
                    errors.append(f"Column '{col}' has type {actual_type}, "
                                f"expected {expected_type} and cannot convert")

    if errors:
        raise ValueError(f"CSV validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return df
```

### BAM/SAM Validation

```python
import subprocess

def validate_bam_file(bam_path, check_sorted=True, check_indexed=False):
    """
    Validate BAM/SAM file integrity

    Uses samtools for validation
    """

    # Check file exists and is not empty
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    file_size = os.path.getsize(bam_path)
    if file_size < 1000:
        raise ValueError(f"BAM file too small ({file_size} bytes): {bam_path}")

    # Check BAM integrity with samtools quickcheck
    result = subprocess.run(
        ['samtools', 'quickcheck', bam_path],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise ValueError(f"BAM file is corrupted: {bam_path}\n{result.stderr}")

    # Check if sorted (if required)
    if check_sorted:
        result = subprocess.run(
            ['samtools', 'view', '-H', bam_path],
            capture_output=True,
            text=True
        )

        header = result.stdout
        if 'SO:coordinate' not in header and 'SO:queryname' not in header:
            raise ValueError(f"BAM file is not sorted: {bam_path}")

    # Check if indexed (if required)
    if check_indexed:
        index_path = f"{bam_path}.bai"
        if not os.path.exists(index_path):
            raise ValueError(f"BAM index not found: {index_path}")

    # Get basic statistics
    result = subprocess.run(
        ['samtools', 'flagstat', bam_path],
        capture_output=True,
        text=True
    )

    stats = result.stdout

    # Extract read count
    total_reads = int(stats.split('\n')[0].split()[0])
    if total_reads == 0:
        raise ValueError(f"BAM file contains no reads: {bam_path}")

    return {'total_reads': total_reads, 'file_size_gb': file_size / (1024**3)}
```

## Schema Validation

### Pandas DataFrame Schema Validation

```python
def validate_dataframe_schema(df, schema):
    """
    Validate DataFrame against a schema definition

    Args:
        df: pandas DataFrame
        schema: Dict with column definitions

    Example schema:
        {
            'sample_id': {
                'dtype': 'object',
                'nullable': False,
                'unique': True,
                'pattern': r'^sample_\d+$'
            },
            'age': {
                'dtype': 'int64',
                'nullable': False,
                'min': 0,
                'max': 120
            },
            'score': {
                'dtype': 'float64',
                'nullable': True,
                'min': 0.0,
                'max': 1.0
            }
        }
    """
    import re

    errors = []

    for col, constraints in schema.items():
        # Check column exists
        if col not in df.columns:
            errors.append(f"Missing column: {col}")
            continue

        series = df[col]

        # Check dtype
        if 'dtype' in constraints:
            expected_dtype = constraints['dtype']
            if not pd.api.types.is_dtype_equal(series.dtype, expected_dtype):
                errors.append(f"Column '{col}': expected dtype {expected_dtype}, "
                            f"got {series.dtype}")

        # Check nullable
        if not constraints.get('nullable', True):
            null_count = series.isna().sum()
            if null_count > 0:
                errors.append(f"Column '{col}': {null_count} null values found "
                            f"(null not allowed)")

        # Check unique
        if constraints.get('unique', False):
            duplicates = series.duplicated().sum()
            if duplicates > 0:
                errors.append(f"Column '{col}': {duplicates} duplicate values found "
                            f"(must be unique)")

        # Check min/max for numeric columns
        if 'min' in constraints:
            valid_values = series.dropna()
            below_min = (valid_values < constraints['min']).sum()
            if below_min > 0:
                errors.append(f"Column '{col}': {below_min} values below minimum "
                            f"{constraints['min']}")

        if 'max' in constraints:
            valid_values = series.dropna()
            above_max = (valid_values > constraints['max']).sum()
            if above_max > 0:
                errors.append(f"Column '{col}': {above_max} values above maximum "
                            f"{constraints['max']}")

        # Check pattern for string columns
        if 'pattern' in constraints:
            valid_values = series.dropna()
            pattern = re.compile(constraints['pattern'])
            non_matching = ~valid_values.astype(str).str.match(pattern)
            if non_matching.any():
                count = non_matching.sum()
                errors.append(f"Column '{col}': {count} values don't match pattern "
                            f"{constraints['pattern']}")

        # Check allowed values
        if 'allowed_values' in constraints:
            valid_values = series.dropna()
            allowed = set(constraints['allowed_values'])
            invalid = set(valid_values.unique()) - allowed
            if invalid:
                errors.append(f"Column '{col}': invalid values {invalid}, "
                            f"allowed: {allowed}")

    if errors:
        raise ValueError(f"Schema validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return True
```

## Data Quality Validation

### Missing Data Detection

```python
import numpy as np

def validate_missing_data(df, max_missing_per_column=0.5, max_missing_per_row=0.8):
    """
    Validate missing data doesn't exceed thresholds

    Args:
        df: pandas DataFrame
        max_missing_per_column: Maximum fraction of missing values per column
        max_missing_per_row: Maximum fraction of missing values per row
    """

    errors = []
    warnings = []

    # Check missing per column
    missing_per_col = df.isna().sum() / len(df)
    problematic_cols = missing_per_col[missing_per_col > max_missing_per_column]

    if len(problematic_cols) > 0:
        errors.append(f"Columns with >{max_missing_per_column:.0%} missing data:")
        for col, pct in problematic_cols.items():
            errors.append(f"  - {col}: {pct:.1%}")

    # Check missing per row
    missing_per_row = df.isna().sum(axis=1) / len(df.columns)
    problematic_rows = (missing_per_row > max_missing_per_row).sum()

    if problematic_rows > 0:
        warnings.append(f"{problematic_rows} rows with >{max_missing_per_row:.0%} missing data")

    if errors:
        raise ValueError(f"Missing data validation failed:\n" +
                        "\n".join(errors))

    if warnings:
        print("⚠ Warnings:\n" + "\n".join(f"  - {w}" for w in warnings))

    return True
```

### Numeric Range Validation

```python
def validate_numeric_ranges(df, range_specs):
    """
    Validate numeric columns fall within expected ranges

    Args:
        range_specs: Dict mapping column names to (min, max) tuples

    Example:
        ranges = {
            'age': (0, 120),
            'temperature': (35.0, 42.0),
            'count': (0, None)  # None means no upper limit
        }
    """

    errors = []

    for col, (min_val, max_val) in range_specs.items():
        if col not in df.columns:
            errors.append(f"Column not found: {col}")
            continue

        series = df[col].dropna()

        if min_val is not None:
            below = (series < min_val).sum()
            if below > 0:
                actual_min = series.min()
                errors.append(f"Column '{col}': {below} values below {min_val} "
                            f"(min value: {actual_min})")

        if max_val is not None:
            above = (series > max_val).sum()
            if above > 0:
                actual_max = series.max()
                errors.append(f"Column '{col}': {above} values above {max_val} "
                            f"(max value: {actual_max})")

    if errors:
        raise ValueError(f"Range validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return True
```

### Distribution Validation

```python
def validate_distribution(data, expected_mean=None, expected_std=None,
                         tolerance=0.2, check_normality=False):
    """
    Validate data distribution meets expectations

    Args:
        data: Numeric array or Series
        expected_mean: Expected mean (None to skip)
        expected_std: Expected standard deviation (None to skip)
        tolerance: Relative tolerance (0.2 = 20%)
        check_normality: Whether to test for normality
    """
    from scipy import stats

    data = np.array(data)
    data = data[~np.isnan(data)]  # Remove NaN

    if len(data) == 0:
        raise ValueError("No valid data points")

    actual_mean = np.mean(data)
    actual_std = np.std(data)

    errors = []

    # Check mean
    if expected_mean is not None:
        if abs(actual_mean - expected_mean) / expected_mean > tolerance:
            errors.append(f"Mean {actual_mean:.2f} differs from expected "
                         f"{expected_mean:.2f} by >{tolerance:.0%}")

    # Check standard deviation
    if expected_std is not None:
        if abs(actual_std - expected_std) / expected_std > tolerance:
            errors.append(f"Std {actual_std:.2f} differs from expected "
                         f"{expected_std:.2f} by >{tolerance:.0%}")

    # Check for normality
    if check_normality and len(data) >= 20:
        statistic, p_value = stats.shapiro(data[:5000])  # Limit for performance
        if p_value < 0.05:
            errors.append(f"Data not normally distributed (Shapiro-Wilk p={p_value:.4f})")

    if errors:
        raise ValueError(f"Distribution validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return {'mean': actual_mean, 'std': actual_std}
```

## Bioinformatics-Specific Validation

### Gene Count Matrix Validation

```python
def validate_count_matrix(counts_df, min_total_counts=1e6, min_detected_genes=10000,
                         max_zero_genes_frac=0.5):
    """
    Validate gene count matrix quality

    Args:
        counts_df: DataFrame with genes as rows, samples as columns
        min_total_counts: Minimum total counts per sample
        min_detected_genes: Minimum detected genes per sample
        max_zero_genes_frac: Maximum fraction of all-zero genes
    """

    errors = []
    warnings = []

    n_genes, n_samples = counts_df.shape

    # Check counts are non-negative
    if (counts_df < 0).any().any():
        errors.append("Count matrix contains negative values")

    # Check for non-numeric values
    if not np.issubdtype(counts_df.values.dtype, np.number):
        errors.append(f"Count matrix contains non-numeric values: {counts_df.values.dtype}")

    # Check for NaN or Inf
    if counts_df.isna().any().any():
        nan_counts = counts_df.isna().sum().sum()
        errors.append(f"Count matrix contains {nan_counts} NaN values")

    if np.isinf(counts_df.values).any():
        inf_counts = np.isinf(counts_df.values).sum()
        errors.append(f"Count matrix contains {inf_counts} Inf values")

    # Per-sample validation
    for sample in counts_df.columns:
        sample_counts = counts_df[sample]

        # Total counts
        total = sample_counts.sum()
        if total < min_total_counts:
            errors.append(f"Sample '{sample}': low total counts {total:,.0f} "
                         f"(expected >{min_total_counts:,.0f})")

        # Detected genes
        detected = (sample_counts > 0).sum()
        if detected < min_detected_genes:
            warnings.append(f"Sample '{sample}': low detected genes {detected:,} "
                          f"(expected >{min_detected_genes:,})")

    # All-zero genes
    all_zero = (counts_df.sum(axis=1) == 0).sum()
    zero_frac = all_zero / n_genes
    if zero_frac > max_zero_genes_frac:
        warnings.append(f"{all_zero:,} genes ({zero_frac:.1%}) have zero counts "
                       f"across all samples")

    if errors:
        raise ValueError(f"Count matrix validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    if warnings:
        print("⚠ Warnings:\n" + "\n".join(f"  - {w}" for w in warnings))

    return {
        'n_genes': n_genes,
        'n_samples': n_samples,
        'total_counts_per_sample': counts_df.sum(axis=0).to_dict(),
        'detected_genes_per_sample': (counts_df > 0).sum(axis=0).to_dict()
    }
```

### Genome Coordinate Validation

```python
def validate_genomic_coordinates(df, chromosome_col='chr', start_col='start',
                                end_col='end', reference_chroms=None):
    """
    Validate genomic coordinates are well-formed

    Args:
        df: DataFrame with genomic coordinates
        chromosome_col: Column name for chromosome
        start_col: Column name for start position
        end_col: Column name for end position
        reference_chroms: Set of valid chromosome names (None to skip check)
    """

    errors = []

    # Check required columns exist
    for col in [chromosome_col, start_col, end_col]:
        if col not in df.columns:
            errors.append(f"Missing required column: {col}")

    if errors:
        raise ValueError(f"Coordinate validation failed:\n" + "\n".join(errors))

    # Check start < end
    invalid_ranges = df[start_col] >= df[end_col]
    if invalid_ranges.any():
        count = invalid_ranges.sum()
        errors.append(f"{count} records with start >= end")

    # Check coordinates are non-negative
    if (df[start_col] < 0).any():
        errors.append(f"Negative start coordinates found")
    if (df[end_col] < 0).any():
        errors.append(f"Negative end coordinates found")

    # Check chromosomes are valid
    if reference_chroms is not None:
        invalid_chroms = set(df[chromosome_col].unique()) - set(reference_chroms)
        if invalid_chroms:
            errors.append(f"Invalid chromosome names: {invalid_chroms}")

    # Check for suspiciously large ranges
    ranges = df[end_col] - df[start_col]
    if (ranges > 10_000_000).any():  # >10 Mb
        count = (ranges > 10_000_000).sum()
        errors.append(f"{count} records with suspiciously large ranges (>10 Mb)")

    if errors:
        raise ValueError(f"Coordinate validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return True
```

## Cross-Stage Validation

### Sample Consistency Validation

```python
def validate_sample_consistency(*dataframes, sample_id_col='sample_id'):
    """
    Validate sample IDs are consistent across multiple datasets

    Args:
        *dataframes: Multiple DataFrames to check
        sample_id_col: Column name containing sample IDs
    """

    if len(dataframes) < 2:
        raise ValueError("Need at least 2 dataframes to check consistency")

    # Extract sample sets
    sample_sets = []
    for i, df in enumerate(dataframes):
        if sample_id_col not in df.columns:
            raise ValueError(f"DataFrame {i} missing column: {sample_id_col}")
        sample_sets.append(set(df[sample_id_col].unique()))

    # Find intersection (common samples)
    common_samples = sample_sets[0]
    for s in sample_sets[1:]:
        common_samples = common_samples.intersection(s)

    # Report differences
    errors = []
    for i, sample_set in enumerate(sample_sets):
        missing = sample_set - common_samples
        if missing:
            errors.append(f"Dataset {i} has {len(missing)} samples not in all datasets: "
                         f"{list(missing)[:5]}")

    if errors:
        raise ValueError(f"Sample consistency validation failed:\n" +
                        "\n".join(f"  - {e}" for e in errors))

    return common_samples
```

### Dimension Consistency Validation

```python
def validate_dimension_consistency(data_dict):
    """
    Validate dimensions are consistent across related datasets

    Args:
        data_dict: Dict mapping names to DataFrames with expected relationships

    Example:
        data = {
            'counts': counts_df,  # genes x samples
            'metadata': metadata_df,  # samples x features
            'gene_info': gene_info_df  # genes x info
        }

        relationships = {
            ('counts', 'metadata'): ('columns', 'index'),  # counts cols = metadata rows
            ('counts', 'gene_info'): ('index', 'index')    # counts rows = gene_info rows
        }
    """

    def validate_relationship(data_dict, name1, dim1, name2, dim2):
        """Validate one dimension relationship"""

        df1 = data_dict[name1]
        df2 = data_dict[name2]

        # Get dimensions
        if dim1 == 'index':
            ids1 = set(df1.index)
        elif dim1 == 'columns':
            ids1 = set(df1.columns)
        else:
            raise ValueError(f"Unknown dimension: {dim1}")

        if dim2 == 'index':
            ids2 = set(df2.index)
        elif dim2 == 'columns':
            ids2 = set(df2.columns)
        else:
            raise ValueError(f"Unknown dimension: {dim2}")

        # Check consistency
        missing_from_2 = ids1 - ids2
        missing_from_1 = ids2 - ids1

        errors = []
        if missing_from_2:
            errors.append(f"{name1}.{dim1} has {len(missing_from_2)} IDs not in {name2}.{dim2}")
        if missing_from_1:
            errors.append(f"{name2}.{dim2} has {len(missing_from_1)} IDs not in {name1}.{dim1}")

        return errors

    # Example usage would define relationships and validate them
    # This is a framework - specific relationships depend on use case

    return True
```

## Summary

This reference provides reusable validation patterns for common pipeline scenarios. Key principles:

1. **Validate early** - Check inputs before expensive processing
2. **Be specific** - Provide actionable error messages
3. **Be comprehensive** - Check multiple aspects (existence, format, content, quality)
4. **Use appropriate thresholds** - Based on domain knowledge and historical data
5. **Distinguish errors from warnings** - Not all issues are fatal
6. **Log validation results** - Track metrics over time

These patterns can be adapted and combined for specific pipeline requirements.
