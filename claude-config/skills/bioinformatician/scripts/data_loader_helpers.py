#!/usr/bin/env python3
"""
Data Loading Helper Functions

Provides standardized functions for loading common bioinformatics file formats.
These are reference implementations that Claude can:
1. Read to understand proper data loading patterns
2. Copy/adapt for specific analyses
3. Call directly if appropriate

Author: David Angeles Albores
Date: 2026-01-28
"""

import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import logging
from typing import Union, Optional, Dict, List

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def load_count_matrix(
    filepath: Union[str, Path],
    sep: str = '\t',
    index_col: int = 0,
    compression: Optional[str] = 'infer'
) -> pd.DataFrame:
    """
    Load gene expression count matrix (rows = genes, columns = samples).

    Parameters
    ----------
    filepath : str or Path
        Path to count matrix file (CSV, TSV, or compressed)
    sep : str, default '\t'
        Delimiter (tab, comma, space)
    index_col : int, default 0
        Column to use as row index (gene names)
    compression : str or None, default 'infer'
        Compression type ('gzip', 'bz2', 'zip', 'xz', or 'infer')

    Returns
    -------
    pd.DataFrame
        Count matrix with genes as rows, samples as columns

    Examples
    --------
    >>> counts = load_count_matrix('counts.tsv')
    >>> counts = load_count_matrix('counts.csv.gz', sep=',')
    """
    logger.info(f"Loading count matrix from {filepath}")

    # Validate file exists
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    # Load data
    try:
        df = pd.read_csv(
            filepath,
            sep=sep,
            index_col=index_col,
            compression=compression
        )

        # Validate numeric data
        if not df.apply(lambda col: pd.api.types.is_numeric_dtype(col)).all():
            logger.warning("Non-numeric columns detected. Check data format.")

        logger.info(f"Loaded {df.shape[0]} genes × {df.shape[1]} samples")
        return df

    except Exception as e:
        logger.error(f"Failed to load count matrix: {e}")
        raise


def load_metadata(
    filepath: Union[str, Path],
    index_col: int = 0,
    required_columns: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Load sample metadata (experimental design).

    Parameters
    ----------
    filepath : str or Path
        Path to metadata file (usually CSV or TSV)
    index_col : int, default 0
        Column containing sample IDs (must match count matrix columns)
    required_columns : list of str, optional
        Columns that must be present

    Returns
    -------
    pd.DataFrame
        Metadata with samples as rows

    Examples
    --------
    >>> metadata = load_metadata('samples.csv', required_columns=['condition', 'batch'])
    """
    logger.info(f"Loading metadata from {filepath}")

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    # Infer separator
    sep = ',' if filepath.suffix == '.csv' else '\t'

    metadata = pd.read_csv(filepath, sep=sep, index_col=index_col)

    # Validate required columns
    if required_columns:
        missing = set(required_columns) - set(metadata.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

    logger.info(f"Loaded metadata for {len(metadata)} samples")
    logger.info(f"Columns: {list(metadata.columns)}")

    return metadata


def validate_count_metadata_match(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    strict: bool = True
) -> Dict[str, any]:
    """
    Validate that count matrix and metadata have matching samples.

    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (samples as columns)
    metadata : pd.DataFrame
        Metadata (samples as rows)
    strict : bool, default True
        If True, require exact match. If False, subset to common samples.

    Returns
    -------
    dict
        Validation results and alignment info

    Raises
    ------
    ValueError
        If samples don't match and strict=True

    Examples
    --------
    >>> validation = validate_count_metadata_match(counts, metadata)
    >>> if not validation['exact_match']:
    ...     print(validation['missing_in_counts'])
    """
    count_samples = set(counts.columns)
    metadata_samples = set(metadata.index)

    common = count_samples & metadata_samples
    only_counts = count_samples - metadata_samples
    only_metadata = metadata_samples - count_samples

    exact_match = (count_samples == metadata_samples)

    logger.info(f"Samples in counts: {len(count_samples)}")
    logger.info(f"Samples in metadata: {len(metadata_samples)}")
    logger.info(f"Common samples: {len(common)}")

    if not exact_match:
        logger.warning(f"Samples only in counts: {len(only_counts)}")
        logger.warning(f"Samples only in metadata: {len(only_metadata)}")

        if strict:
            raise ValueError(
                f"Sample mismatch! "
                f"Only in counts: {only_counts}, "
                f"Only in metadata: {only_metadata}"
            )

    return {
        'exact_match': exact_match,
        'common_samples': sorted(common),
        'n_common': len(common),
        'missing_in_counts': sorted(only_metadata),
        'missing_in_metadata': sorted(only_counts),
    }


def load_gene_annotations(
    filepath: Union[str, Path],
    gene_id_col: str = 'gene_id',
    gene_name_col: Optional[str] = 'gene_name'
) -> pd.DataFrame:
    """
    Load gene annotations (ID, name, biotype, chromosome, etc.).

    Parameters
    ----------
    filepath : str or Path
        Path to annotation file
    gene_id_col : str, default 'gene_id'
        Column containing gene IDs (used as index)
    gene_name_col : str or None, default 'gene_name'
        Column containing gene symbols

    Returns
    -------
    pd.DataFrame
        Gene annotations with gene IDs as index

    Examples
    --------
    >>> annotations = load_gene_annotations('genes.tsv')
    >>> gpcr_genes = annotations[annotations['gene_name'].str.contains('GPR')]
    """
    logger.info(f"Loading gene annotations from {filepath}")

    filepath = Path(filepath)
    sep = ',' if filepath.suffix == '.csv' else '\t'

    annotations = pd.read_csv(filepath, sep=sep)

    # Validate required columns
    if gene_id_col not in annotations.columns:
        raise ValueError(f"Column '{gene_id_col}' not found in annotations")

    # Set index
    annotations = annotations.set_index(gene_id_col)

    logger.info(f"Loaded annotations for {len(annotations)} genes")
    logger.info(f"Columns: {list(annotations.columns)}")

    return annotations


def load_fasta(filepath: Union[str, Path]) -> Dict[str, str]:
    """
    Load sequences from FASTA file (plain or gzipped).

    Parameters
    ----------
    filepath : str or Path
        Path to FASTA file

    Returns
    -------
    dict
        Dictionary mapping sequence IDs to sequences

    Examples
    --------
    >>> sequences = load_fasta('genome.fa.gz')
    >>> chrI_seq = sequences['chrI']

    Notes
    -----
    For large genomes, consider using Biopython's SeqIO.index() instead
    for memory-efficient random access.
    """
    logger.info(f"Loading FASTA from {filepath}")

    filepath = Path(filepath)

    # Determine if gzipped
    open_func = gzip.open if filepath.suffix == '.gz' else open
    mode = 'rt' if filepath.suffix == '.gz' else 'r'

    sequences = {}
    current_id = None
    current_seq = []

    with open_func(filepath, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)

                # Start new sequence
                current_id = line[1:].split()[0]  # First word after '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    logger.info(f"Loaded {len(sequences)} sequences")

    return sequences


def load_gff(
    filepath: Union[str, Path],
    feature_type: Optional[str] = None
) -> pd.DataFrame:
    """
    Load genomic features from GFF/GTF file.

    Parameters
    ----------
    filepath : str or Path
        Path to GFF/GTF file
    feature_type : str, optional
        Filter for specific feature type (e.g., 'gene', 'exon', 'CDS')

    Returns
    -------
    pd.DataFrame
        Genomic features with standard GFF columns

    Examples
    --------
    >>> genes = load_gff('annotations.gff', feature_type='gene')
    >>> genes_on_chrX = genes[genes['seqid'] == 'chrX']
    """
    logger.info(f"Loading GFF from {filepath}")

    filepath = Path(filepath)
    open_func = gzip.open if filepath.suffix == '.gz' else open
    mode = 'rt' if filepath.suffix == '.gz' else 'r'

    # GFF column names
    columns = [
        'seqid', 'source', 'type', 'start', 'end',
        'score', 'strand', 'phase', 'attributes'
    ]

    rows = []
    with open_func(filepath, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip comments

            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue  # Skip malformed lines

            # Filter by feature type if specified
            if feature_type and fields[2] != feature_type:
                continue

            rows.append(fields)

    gff = pd.DataFrame(rows, columns=columns)

    # Convert numeric columns
    gff['start'] = pd.to_numeric(gff['start'])
    gff['end'] = pd.to_numeric(gff['end'])

    logger.info(f"Loaded {len(gff)} features")
    if feature_type:
        logger.info(f"Filtered for feature type: {feature_type}")

    return gff


def quick_data_summary(df: pd.DataFrame) -> None:
    """
    Print quick summary of DataFrame for initial inspection.

    Parameters
    ----------
    df : pd.DataFrame
        Data to summarize

    Examples
    --------
    >>> counts = load_count_matrix('counts.tsv')
    >>> quick_data_summary(counts)
    """
    print("=" * 60)
    print("DATA SUMMARY")
    print("=" * 60)
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"\nFirst 5 row names:\n{df.index[:5].tolist()}")
    print(f"\nFirst 5 column names:\n{df.columns[:5].tolist()}")
    print(f"\nData types:\n{df.dtypes.value_counts()}")
    print(f"\nMissing values:\n{df.isnull().sum().sum()} total ({df.isnull().sum().sum() / df.size * 100:.2f}%)")
    print(f"\nNumeric summary:\n{df.describe()}")
    print("=" * 60)


# Example usage (for Claude to reference)
if __name__ == '__main__':
    """
    Example workflow for loading and validating data.
    Claude can adapt this pattern for specific analyses.
    """

    # Example file paths (user would provide actual paths)
    count_file = 'data/counts.tsv'
    metadata_file = 'data/samples.csv'

    # Load data
    counts = load_count_matrix(count_file)
    metadata = load_metadata(metadata_file, required_columns=['condition', 'replicate'])

    # Validate match
    validation = validate_count_metadata_match(counts, metadata, strict=False)

    if not validation['exact_match']:
        # Subset to common samples
        common = validation['common_samples']
        counts = counts[common]
        metadata = metadata.loc[common]
        logger.info(f"Subsetted to {len(common)} common samples")

    # Quick inspection
    quick_data_summary(counts)
    quick_data_summary(metadata)

    logger.info("Data loading complete!")
