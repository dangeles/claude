# RNA-seq Pipeline Design Example

## Overview

This example demonstrates a complete RNA-seq analysis pipeline design following Data Pipeline Manager best practices. The pipeline processes raw FASTQ files through quality control, alignment, quantification, and produces a count matrix ready for differential expression analysis.

## Pipeline Architecture

### Pipeline Pattern
Hybrid: Sequential stages with parallel sample processing

### Input Specifications

**Required Files:**
- Raw FASTQ files (paired-end): `{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz`
- Sample metadata: `samples.csv`
- Reference genome: `genome.fa`
- Gene annotation: `genes.gtf`

**Sample Metadata Format:**
```csv
sample_id,condition,batch,replicate,sequencing_run
sample_001,control,batch1,rep1,run_2026_01_15
sample_002,control,batch1,rep2,run_2026_01_15
sample_003,treatment,batch1,rep1,run_2026_01_15
sample_004,treatment,batch1,rep2,run_2026_01_15
```

### Output Specifications

**Final Outputs:**
- Count matrix: `counts_matrix.csv` (genes × samples)
- QC report: `multiqc_report.html`
- Quality metrics: `qc_metrics.json`
- Processing log: `pipeline.log`
- Sample statistics: `sample_stats.csv`

### Resource Requirements

**Per Sample:**
- CPU: 8 cores
- Memory: 32 GB
- Disk: 50 GB temporary space
- Time: ~2 hours per sample

**Total for 20 samples:**
- Processing time: 4-5 hours (parallel processing)
- Storage: 1 TB temporary + 200 GB permanent

## Detailed Pipeline Stages

### Stage 0: Pre-Flight Validation

**Purpose:** Validate all inputs before starting expensive processing

**Validation Checks:**

```python
import os
import pandas as pd
import gzip

def validate_inputs(config):
    """Comprehensive input validation"""
    errors = []

    # 1. Check metadata file exists and is valid
    if not os.path.exists(config['metadata_file']):
        errors.append(f"Metadata file not found: {config['metadata_file']}")
    else:
        try:
            metadata = pd.read_csv(config['metadata_file'])
            required_cols = ['sample_id', 'condition', 'batch']
            missing = set(required_cols) - set(metadata.columns)
            if missing:
                errors.append(f"Missing metadata columns: {missing}")

            # Check for duplicate sample IDs
            if metadata['sample_id'].duplicated().any():
                dups = metadata[metadata['sample_id'].duplicated()]['sample_id'].tolist()
                errors.append(f"Duplicate sample IDs: {dups}")
        except Exception as e:
            errors.append(f"Cannot read metadata: {e}")

    # 2. Check reference files exist
    if not os.path.exists(config['genome_fa']):
        errors.append(f"Reference genome not found: {config['genome_fa']}")
    if not os.path.exists(config['annotation_gtf']):
        errors.append(f"Annotation file not found: {config['annotation_gtf']}")

    # 3. Check FASTQ files for each sample
    for sample_id in metadata['sample_id']:
        r1_file = f"{config['fastq_dir']}/{sample_id}_R1.fastq.gz"
        r2_file = f"{config['fastq_dir']}/{sample_id}_R2.fastq.gz"

        for fastq_file in [r1_file, r2_file]:
            if not os.path.exists(fastq_file):
                errors.append(f"FASTQ file not found: {fastq_file}")
            else:
                # Check file is not empty
                file_size = os.path.getsize(fastq_file)
                if file_size < 1000:  # Very small, likely empty
                    errors.append(f"FASTQ file too small: {fastq_file} ({file_size} bytes)")

                # Check file is valid gzip
                try:
                    with gzip.open(fastq_file, 'rt') as f:
                        first_line = f.readline()
                        if not first_line.startswith('@'):
                            errors.append(f"Invalid FASTQ format: {fastq_file}")
                except Exception as e:
                    errors.append(f"Cannot read FASTQ file {fastq_file}: {e}")

    # 4. Check output directory is writable
    if not os.access(config['output_dir'], os.W_OK):
        errors.append(f"Output directory not writable: {config['output_dir']}")

    # 5. Check sufficient disk space
    available_space = shutil.disk_usage(config['output_dir']).free
    required_space = len(metadata) * 50 * 1024**3  # 50 GB per sample
    if available_space < required_space:
        errors.append(f"Insufficient disk space: {available_space / 1024**3:.1f} GB available, "
                     f"{required_space / 1024**3:.1f} GB required")

    if errors:
        raise ValueError(f"Input validation failed:\n" + "\n".join(f"  - {e}" for e in errors))

    print("✓ All input validation checks passed")
```

**Checkpoint:** Create validation report before proceeding

### Stage 1: Quality Control (FastQC)

**Purpose:** Assess raw read quality and identify potential issues

**Processing:**

```bash
# For each sample, run FastQC on both R1 and R2
fastqc \
    --outdir ${output_dir}/fastqc_raw \
    --threads 4 \
    ${sample}_R1.fastq.gz \
    ${sample}_R2.fastq.gz
```

**Validation Checks:**

```python
def validate_fastqc_results(fastqc_data_file):
    """Parse FastQC results and check quality metrics"""
    warnings = []
    failures = []

    with open(fastqc_data_file) as f:
        for line in f:
            if line.startswith('>>'):
                module = line.split('\t')[0][2:]
                status = line.split('\t')[1]

                if status == 'FAIL':
                    failures.append(module)
                elif status == 'WARN':
                    warnings.append(module)

            # Check per base sequence quality
            if line.startswith('>>Per base sequence quality'):
                # Parse quality scores
                scores = []
                for qline in f:
                    if qline.startswith('>>'):
                        break
                    if not qline.startswith('#'):
                        score = float(qline.split('\t')[1])
                        scores.append(score)

                median_quality = statistics.median(scores)
                if median_quality < 28:
                    failures.append(f"Low median quality score: {median_quality:.1f}")

    # Report results
    if failures:
        raise ValueError(f"FastQC quality check failed: {failures}")
    if warnings:
        print(f"FastQC warnings (non-critical): {warnings}")

    print("✓ FastQC quality checks passed")
```

**Quality Metrics to Track:**
- Per base sequence quality (should be >28)
- Per sequence quality scores (median >30)
- Adapter content (should be minimal)
- Overrepresented sequences (check for contamination)
- GC content (should match expected for species)

**Error Handling:**
- If quality is very poor (median <20), flag sample for review
- If adapter contamination is high (>10%), proceed to trimming
- Log all quality metrics for downstream analysis

### Stage 2: Read Trimming (Trimmomatic)

**Purpose:** Remove adapters and low-quality bases

**Processing:**

```bash
trimmomatic PE \
    -threads 8 \
    -phred33 \
    ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    ${sample}_R1_trimmed.fastq.gz ${sample}_R1_unpaired.fastq.gz \
    ${sample}_R2_trimmed.fastq.gz ${sample}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

**Validation Checks:**

```python
def validate_trimming_results(sample_id, raw_r1, trimmed_r1, log_file):
    """Validate trimming preserved sufficient reads"""

    # Parse trimming statistics from log
    with open(log_file) as f:
        log_content = f.read()

    # Extract input and output read counts
    input_reads = int(re.search(r'Input Read Pairs: (\d+)', log_content).group(1))
    surviving_reads = int(re.search(r'Both Surviving: (\d+)', log_content).group(1))

    survival_rate = surviving_reads / input_reads

    # Check survival rate is reasonable
    if survival_rate < 0.5:
        raise ValueError(f"Low read survival rate for {sample_id}: {survival_rate:.1%} "
                        f"({surviving_reads}/{input_reads} reads)")

    if survival_rate < 0.7:
        print(f"⚠ Warning: Moderate read loss for {sample_id}: {survival_rate:.1%}")

    # Check trimmed files are not empty
    if os.path.getsize(trimmed_r1) < 1000:
        raise ValueError(f"Trimmed FASTQ is too small: {trimmed_r1}")

    print(f"✓ Trimming validation passed: {survival_rate:.1%} reads retained")

    return {
        'sample_id': sample_id,
        'input_reads': input_reads,
        'surviving_reads': surviving_reads,
        'survival_rate': survival_rate
    }
```

**Checkpoint:** Save trimming statistics for all samples

### Stage 3: Alignment (STAR)

**Purpose:** Align reads to reference genome

**Index Preparation (one-time):**

```bash
STAR \
    --runMode genomeGenerate \
    --genomeDir ${genome_index_dir} \
    --genomeFastaFiles ${genome_fa} \
    --sjdbGTFfile ${annotation_gtf} \
    --sjdbOverhang 99 \
    --runThreadN 16
```

**Alignment:**

```bash
STAR \
    --genomeDir ${genome_index_dir} \
    --readFilesIn ${sample}_R1_trimmed.fastq.gz ${sample}_R2_trimmed.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ${output_dir}/${sample}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --runThreadN 8 \
    --quantMode GeneCounts
```

**Validation Checks:**

```python
def validate_alignment(sample_id, log_file, bam_file):
    """Validate alignment quality metrics"""

    # Parse STAR log file
    with open(log_file) as f:
        log_lines = f.readlines()

    metrics = {}
    for line in log_lines:
        if 'Number of input reads' in line:
            metrics['input_reads'] = int(line.split('|')[1].strip())
        elif 'Uniquely mapped reads number' in line:
            metrics['uniquely_mapped'] = int(line.split('|')[1].strip())
        elif 'Uniquely mapped reads %' in line:
            metrics['uniquely_mapped_pct'] = float(line.split('|')[1].strip().rstrip('%'))
        elif 'Number of reads mapped to multiple loci' in line:
            metrics['multi_mapped'] = int(line.split('|')[1].strip())
        elif '% of reads unmapped: too short' in line:
            metrics['unmapped_too_short_pct'] = float(line.split('|')[1].strip().rstrip('%'))

    # Validate alignment rate
    if metrics['uniquely_mapped_pct'] < 50:
        raise ValueError(f"Low alignment rate for {sample_id}: "
                        f"{metrics['uniquely_mapped_pct']:.1f}% "
                        "(expected >70%)")

    if metrics['uniquely_mapped_pct'] < 70:
        print(f"⚠ Warning: Suboptimal alignment rate for {sample_id}: "
              f"{metrics['uniquely_mapped_pct']:.1f}%")

    # Check BAM file exists and is not empty
    if not os.path.exists(bam_file):
        raise ValueError(f"BAM file not created: {bam_file}")

    bam_size = os.path.getsize(bam_file)
    if bam_size < 1_000_000:  # Less than 1 MB is suspicious
        raise ValueError(f"BAM file too small: {bam_size} bytes")

    # Index BAM file for downstream use
    subprocess.run(['samtools', 'index', bam_file], check=True)

    print(f"✓ Alignment validation passed: {metrics['uniquely_mapped_pct']:.1f}% aligned")

    return metrics
```

**Error Handling:**
- If alignment rate <50%, stop and investigate (wrong reference, poor quality)
- If alignment rate 50-70%, flag but continue (may be acceptable for some datasets)
- Retry once if STAR crashes (may be transient memory issue)
- Log all alignment metrics for QC report

**Checkpoint:** Save alignment metrics before quantification

### Stage 4: Quantification (featureCounts)

**Purpose:** Count reads per gene

**Processing:**

```bash
featureCounts \
    -p \
    -B \
    -C \
    -T 8 \
    -a ${annotation_gtf} \
    -o ${output_dir}/counts.txt \
    ${output_dir}/*_Aligned.sortedByCoord.out.bam
```

**Validation Checks:**

```python
def validate_counts(counts_file, expected_samples):
    """Validate count matrix"""

    # Read counts file
    counts = pd.read_csv(counts_file, sep='\t', comment='#')

    # Check all samples are present
    sample_columns = [col for col in counts.columns if col.endswith('.bam')]
    if len(sample_columns) != expected_samples:
        raise ValueError(f"Expected {expected_samples} samples, found {len(sample_columns)}")

    # Check number of genes
    n_genes = len(counts)
    if n_genes < 10000:
        print(f"⚠ Warning: Low number of genes detected: {n_genes} (expected >20,000 for human)")

    # Check for all-zero genes
    numeric_cols = counts.select_dtypes(include=[np.number]).columns
    all_zero_genes = (counts[numeric_cols] == 0).all(axis=1).sum()
    if all_zero_genes / n_genes > 0.5:
        print(f"⚠ Warning: {all_zero_genes}/{n_genes} genes have zero counts across all samples")

    # Calculate per-sample statistics
    sample_stats = []
    for col in sample_columns:
        sample_counts = counts[col]
        stats = {
            'sample': col.replace('_Aligned.sortedByCoord.out.bam', ''),
            'total_counts': sample_counts.sum(),
            'detected_genes': (sample_counts > 0).sum(),
            'median_count': sample_counts.median(),
            'mean_count': sample_counts.mean()
        }

        # Validate total counts
        if stats['total_counts'] < 1_000_000:
            raise ValueError(f"Low total counts for {stats['sample']}: "
                           f"{stats['total_counts']:,} (expected >5M for RNA-seq)")

        # Validate detected genes
        if stats['detected_genes'] < 10_000:
            print(f"⚠ Warning: Low detected genes for {stats['sample']}: "
                  f"{stats['detected_genes']:,}")

        sample_stats.append(stats)

    # Check coefficient of variation across samples
    total_counts = [s['total_counts'] for s in sample_stats]
    cv = np.std(total_counts) / np.mean(total_counts)
    if cv > 0.5:
        print(f"⚠ Warning: High variation in library sizes (CV={cv:.2f})")

    print("✓ Count matrix validation passed")
    return pd.DataFrame(sample_stats)
```

**Checkpoint:** Save count matrix and sample statistics

### Stage 5: Multi-Sample QC (MultiQC)

**Purpose:** Aggregate QC metrics across all samples

**Processing:**

```bash
multiqc \
    --outdir ${output_dir}/multiqc \
    --force \
    --filename multiqc_report.html \
    ${output_dir}
```

**Validation:**

```python
def validate_multiqc_report(report_dir, expected_samples):
    """Check MultiQC report was generated successfully"""

    report_html = os.path.join(report_dir, 'multiqc_report.html')
    if not os.path.exists(report_html):
        raise ValueError("MultiQC report not generated")

    # Parse MultiQC data
    data_json = os.path.join(report_dir, 'multiqc_data', 'multiqc_data.json')
    if os.path.exists(data_json):
        with open(data_json) as f:
            data = json.load(f)

        # Check all samples are in report
        samples_in_report = len(data.get('report_saved_raw_data', {}).get('multiqc_fastqc', {}))
        if samples_in_report < expected_samples:
            print(f"⚠ Warning: Only {samples_in_report}/{expected_samples} samples in QC report")

    print("✓ MultiQC report generated successfully")
```

### Stage 6: Final Assembly and Validation

**Purpose:** Create final output files and comprehensive validation report

**Processing:**

```python
def create_final_outputs(config, sample_stats, alignment_metrics, trim_stats):
    """Assemble final outputs and validation report"""

    # 1. Create cleaned count matrix
    counts = pd.read_csv(f"{config['output_dir']}/counts.txt", sep='\t', comment='#')

    # Rename columns to sample IDs
    counts.columns = [col.replace('_Aligned.sortedByCoord.out.bam', '')
                     if col.endswith('.bam') else col
                     for col in counts.columns]

    # Save clean count matrix
    count_matrix = counts.set_index('Geneid').drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
    count_matrix.to_csv(f"{config['output_dir']}/counts_matrix.csv")

    # 2. Create comprehensive QC metrics file
    qc_metrics = {
        'pipeline_version': config['pipeline_version'],
        'reference_genome': config['genome_build'],
        'annotation_version': config['annotation_version'],
        'processing_date': datetime.now().isoformat(),
        'total_samples': len(sample_stats),
        'samples': []
    }

    for idx, sample_row in sample_stats.iterrows():
        sample_id = sample_row['sample']
        sample_metrics = {
            'sample_id': sample_id,
            'raw_reads': trim_stats[sample_id]['input_reads'],
            'trimmed_reads': trim_stats[sample_id]['surviving_reads'],
            'trim_survival_rate': trim_stats[sample_id]['survival_rate'],
            'aligned_reads': alignment_metrics[sample_id]['uniquely_mapped'],
            'alignment_rate': alignment_metrics[sample_id]['uniquely_mapped_pct'],
            'total_counts': int(sample_row['total_counts']),
            'detected_genes': int(sample_row['detected_genes']),
            'qc_pass': True
        }

        # Determine if sample passes QC thresholds
        if (sample_metrics['alignment_rate'] < 70 or
            sample_metrics['total_counts'] < 5_000_000 or
            sample_metrics['detected_genes'] < 10_000):
            sample_metrics['qc_pass'] = False
            sample_metrics['qc_flags'] = []

            if sample_metrics['alignment_rate'] < 70:
                sample_metrics['qc_flags'].append('low_alignment_rate')
            if sample_metrics['total_counts'] < 5_000_000:
                sample_metrics['qc_flags'].append('low_total_counts')
            if sample_metrics['detected_genes'] < 10_000:
                sample_metrics['qc_flags'].append('low_detected_genes')

        qc_metrics['samples'].append(sample_metrics)

    # Save QC metrics
    with open(f"{config['output_dir']}/qc_metrics.json", 'w') as f:
        json.dump(qc_metrics, f, indent=2)

    # 3. Create sample statistics CSV
    stats_df = pd.DataFrame(qc_metrics['samples'])
    stats_df.to_csv(f"{config['output_dir']}/sample_stats.csv", index=False)

    # 4. Print summary
    n_passed = sum(s['qc_pass'] for s in qc_metrics['samples'])
    n_failed = len(qc_metrics['samples']) - n_passed

    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60)
    print(f"Total samples processed: {len(qc_metrics['samples'])}")
    print(f"Samples passed QC: {n_passed}")
    print(f"Samples failed QC: {n_failed}")

    if n_failed > 0:
        print("\nFailed samples:")
        for sample in qc_metrics['samples']:
            if not sample['qc_pass']:
                print(f"  - {sample['sample_id']}: {', '.join(sample['qc_flags'])}")

    print(f"\nOutput files:")
    print(f"  - Count matrix: {config['output_dir']}/counts_matrix.csv")
    print(f"  - QC metrics: {config['output_dir']}/qc_metrics.json")
    print(f"  - Sample stats: {config['output_dir']}/sample_stats.csv")
    print(f"  - MultiQC report: {config['output_dir']}/multiqc/multiqc_report.html")
    print("="*60)
```

## Error Handling and Recovery

### Retry Strategy

```python
class PipelineError(Exception):
    """Base pipeline error"""
    pass

class TransientError(PipelineError):
    """Recoverable error - should retry"""
    pass

class PermanentError(PipelineError):
    """Non-recoverable error - should not retry"""
    pass

def run_stage_with_retry(stage_func, sample_id, max_retries=3):
    """Run pipeline stage with retry logic"""

    for attempt in range(max_retries):
        try:
            result = stage_func(sample_id)
            return result

        except TransientError as e:
            if attempt == max_retries - 1:
                raise PermanentError(f"Stage failed after {max_retries} attempts: {e}")

            wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
            print(f"Attempt {attempt + 1} failed, retrying in {wait_time}s: {e}")
            time.sleep(wait_time)

        except PermanentError as e:
            # Don't retry permanent errors
            raise

        except Exception as e:
            # Unknown errors - log and don't retry
            logging.error(f"Unexpected error in stage: {e}", exc_info=True)
            raise PermanentError(f"Unexpected error: {e}")
```

### Checkpointing

```python
def run_pipeline_with_checkpoint(samples, config):
    """Run pipeline with checkpointing for resumability"""

    checkpoint_file = f"{config['output_dir']}/pipeline_checkpoint.json"

    # Load checkpoint if exists
    completed_stages = {}
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file) as f:
            completed_stages = json.load(f)
        print(f"Resuming from checkpoint: {len(completed_stages)} samples partially completed")

    stages = ['fastqc', 'trim', 'align', 'count']

    for sample_id in samples:
        sample_completed = completed_stages.get(sample_id, [])

        for stage in stages:
            if stage in sample_completed:
                print(f"Skipping {stage} for {sample_id} (already completed)")
                continue

            # Run stage
            try:
                run_stage(stage, sample_id, config)

                # Update checkpoint
                sample_completed.append(stage)
                completed_stages[sample_id] = sample_completed

                with open(checkpoint_file, 'w') as f:
                    json.dump(completed_stages, f, indent=2)

            except Exception as e:
                logging.error(f"Stage {stage} failed for {sample_id}: {e}")
                raise

    # Clean up checkpoint file on success
    os.remove(checkpoint_file)
```

## Performance Optimization

### Parallel Sample Processing

```python
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_samples_parallel(samples, config, n_workers=4):
    """Process multiple samples in parallel"""

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all samples
        future_to_sample = {
            executor.submit(process_single_sample, sample, config): sample
            for sample in samples
        }

        # Process results as they complete
        results = {}
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                result = future.result()
                results[sample] = result
                print(f"✓ Completed {sample}")
            except Exception as e:
                print(f"✗ Failed {sample}: {e}")
                results[sample] = {'error': str(e)}

        return results
```

## Monitoring Dashboard

```python
def generate_progress_dashboard(samples, completed_stages):
    """Generate real-time progress dashboard"""

    stages = ['fastqc', 'trim', 'align', 'count']

    # Calculate progress
    total_stages = len(samples) * len(stages)
    completed = sum(len(completed_stages.get(s, [])) for s in samples)
    progress_pct = (completed / total_stages) * 100

    # Print dashboard
    print("\n" + "="*60)
    print(f"PIPELINE PROGRESS: {completed}/{total_stages} stages ({progress_pct:.1f}%)")
    print("="*60)

    for sample in samples:
        sample_stages = completed_stages.get(sample, [])
        stage_status = [('✓' if s in sample_stages else '○') for s in stages]
        print(f"{sample}: {' '.join(stage_status)} [{' '.join(stages)}]")

    print("="*60 + "\n")
```

## Summary

This example demonstrates a production-ready RNA-seq pipeline with:

- Comprehensive input validation before expensive processing
- Quality checks at each stage with clear thresholds
- Robust error handling with retry logic
- Checkpointing for resumability after failures
- Parallel sample processing for efficiency
- Detailed QC metrics and reporting
- Final validation ensuring output correctness

The pipeline can be adapted for other bioinformatics workflows by changing the specific tools while maintaining the validation and error handling framework.
