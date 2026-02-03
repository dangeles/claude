# Pipeline Debugging Example

## Scenario: Failed RNA-seq Pipeline

You receive an alert that your RNA-seq pipeline failed overnight. Multiple samples failed at different stages, and you need to diagnose and fix the issues quickly.

## Initial Investigation

### Step 1: Check Pipeline Status

```bash
# Check pipeline log for errors
tail -n 100 pipeline.log

# Output shows:
ERROR: Stage 'align' failed for sample_042
ERROR: Stage 'align' failed for sample_043
ERROR: Stage 'count' failed for sample_037
WARNING: Low alignment rate for sample_039: 45.2%
```

### Step 2: Review Checkpoint Status

```python
import json

# Load checkpoint to see what completed
with open('pipeline_checkpoint.json') as f:
    checkpoint = json.load(f)

# Analyze completion status
stages = ['fastqc', 'trim', 'align', 'count']
for sample_id, completed in checkpoint.items():
    missing = [s for s in stages if s not in completed]
    if missing:
        print(f"{sample_id}: Failed at {missing[0]} (completed: {completed})")

# Output:
# sample_037: Failed at count (completed: ['fastqc', 'trim', 'align'])
# sample_039: Failed at count (completed: ['fastqc', 'trim', 'align'])
# sample_042: Failed at align (completed: ['fastqc', 'trim'])
# sample_043: Failed at align (completed: ['fastqc', 'trim'])
```

## Problem 1: Alignment Failures (sample_042, sample_043)

### Diagnosis Process

**Check alignment log files:**

```bash
# Look at STAR alignment log
cat output/sample_042_Log.final.out

# Key findings:
# Number of input reads |   45231
# Uniquely mapped reads number |   0
# % of reads unmapped: other |   98.23%
```

**Problem identified:** Almost no reads mapping, suggests wrong reference genome.

**Verify reference genome:**

```python
def diagnose_alignment_failure(sample_id, config):
    """Diagnose why alignment failed"""

    log_file = f"{config['output_dir']}/{sample_id}_Log.final.out"

    with open(log_file) as f:
        log_content = f.read()

    # Parse key metrics
    metrics = {}
    for line in log_content.split('\n'):
        if 'Number of input reads' in line:
            metrics['input_reads'] = int(line.split('|')[1].strip())
        elif 'Uniquely mapped reads %' in line:
            metrics['mapped_pct'] = float(line.split('|')[1].strip().rstrip('%'))
        elif '% of reads unmapped: too short' in line:
            metrics['too_short_pct'] = float(line.split('|')[1].strip().rstrip('%'))

    # Diagnose based on metrics
    if metrics['mapped_pct'] < 10:
        if metrics['too_short_pct'] > 50:
            return "DIAGNOSIS: Reads too short after trimming - check trimming parameters"
        else:
            return "DIAGNOSIS: Wrong reference genome or severe contamination"

    elif metrics['mapped_pct'] < 50:
        return "DIAGNOSIS: Poor quality reads or genome mismatch"

    return "DIAGNOSIS: Unknown alignment issue"

# Run diagnosis
print(diagnose_alignment_failure('sample_042', config))
# Output: DIAGNOSIS: Wrong reference genome or severe contamination
```

**Check for species mismatch:**

```python
def check_species_contamination(fastq_file, known_species_refs):
    """Quick alignment to multiple species to identify contamination"""

    results = {}

    # Try aligning subset of reads to different species
    for species, ref_genome in known_species_refs.items():
        # Align first 10,000 reads
        subprocess.run([
            'bowtie2',
            '-x', ref_genome,
            '-U', fastq_file,
            '-u', '10000',  # Only first 10k reads
            '--threads', '4',
            '--no-unal',
            '-S', f'/tmp/{species}_test.sam'
        ], capture_output=True)

        # Count aligned reads
        aligned = subprocess.run(
            ['samtools', 'view', '-c', '-F', '4', f'/tmp/{species}_test.sam'],
            capture_output=True, text=True
        )
        results[species] = int(aligned.stdout.strip())

    # Find best match
    best_species = max(results, key=results.get)
    best_count = results[best_species]

    print(f"\nContamination check results:")
    for species, count in sorted(results.items(), key=lambda x: x[1], reverse=True):
        print(f"  {species}: {count}/10000 reads ({count/100:.1f}%)")

    return best_species, best_count/100

# Run contamination check
species_refs = {
    'human': 'refs/hg38/genome',
    'mouse': 'refs/mm10/genome',
    'rat': 'refs/rn6/genome',
    'ecoli': 'refs/ecoli/genome'
}

detected_species, pct = check_species_contamination(
    'output/sample_042_R1_trimmed.fastq.gz',
    species_refs
)

# Output:
# Contamination check results:
#   mouse: 8234/10000 reads (82.3%)
#   human: 234/10000 reads (2.3%)
#   rat: 156/10000 reads (1.6%)
#   ecoli: 43/10000 reads (0.4%)
#
# Detected species: mouse (82.3% alignment)
```

**Root cause:** Sample metadata indicated human samples, but samples 042 and 043 are actually mouse samples.

### Solution

```python
def fix_species_mismatch(sample_ids, correct_species, config):
    """Re-align samples with correct reference genome"""

    # Update sample metadata
    metadata = pd.read_csv(config['metadata_file'])
    metadata.loc[metadata['sample_id'].isin(sample_ids), 'species'] = correct_species

    # Save corrected metadata
    metadata.to_csv(config['metadata_file'], index=False)

    print(f"Updated {len(sample_ids)} samples to species: {correct_species}")

    # Clear failed alignment results
    for sample_id in sample_ids:
        failed_files = [
            f"{config['output_dir']}/{sample_id}_Aligned.sortedByCoord.out.bam",
            f"{config['output_dir']}/{sample_id}_Log.final.out",
        ]
        for f in failed_files:
            if os.path.exists(f):
                os.remove(f)

    # Update checkpoint to re-run alignment
    with open('pipeline_checkpoint.json') as f:
        checkpoint = json.load(f)

    for sample_id in sample_ids:
        # Remove 'align' and 'count' from completed stages
        checkpoint[sample_id] = [s for s in checkpoint[sample_id]
                                if s not in ['align', 'count']]

    with open('pipeline_checkpoint.json', 'w') as f:
        json.dump(checkpoint, f, indent=2)

    print(f"✓ Ready to re-run alignment for {len(sample_ids)} samples")

# Apply fix
fix_species_mismatch(['sample_042', 'sample_043'], 'mouse', config)

# Re-run pipeline for these samples
run_pipeline_with_checkpoint(['sample_042', 'sample_043'], config)
```

## Problem 2: Quantification Failures (sample_037, sample_039)

### Diagnosis Process

**Check for file existence:**

```bash
# Check if BAM files exist
ls -lh output/sample_037_Aligned.sortedByCoord.out.bam
# File exists: 2.3 GB

# Check if BAM is valid
samtools quickcheck output/sample_037_Aligned.sortedByCoord.out.bam
echo $?
# Output: 1 (file is corrupted)
```

**Problem identified:** BAM files exist but are corrupted.

**Check for disk space issues:**

```python
def check_disk_space_during_failure(log_file):
    """Parse logs for disk space issues"""

    with open(log_file) as f:
        log_content = f.read()

    # Look for common disk space error messages
    error_patterns = [
        'No space left on device',
        'Disk quota exceeded',
        'Cannot allocate memory',
        'write error'
    ]

    found_errors = []
    for pattern in error_patterns:
        if pattern.lower() in log_content.lower():
            found_errors.append(pattern)

    if found_errors:
        return f"DIAGNOSIS: Disk/memory issue - {', '.join(found_errors)}"

    return "DIAGNOSIS: Unknown quantification issue"

# Check logs
diagnosis = check_disk_space_during_failure('output/sample_037_Log.out')
print(diagnosis)
# Output: DIAGNOSIS: Disk/memory issue - No space left on device
```

**Root cause:** Disk filled up during alignment, causing BAM files to be corrupted.

### Solution

```python
def fix_disk_space_and_rerun(sample_ids, config):
    """Clean up space and re-run failed samples"""

    # 1. Clean up corrupted files
    print("Cleaning up corrupted files...")
    for sample_id in sample_ids:
        corrupted_files = [
            f"{config['output_dir']}/{sample_id}_Aligned.sortedByCoord.out.bam",
            f"{config['output_dir']}/{sample_id}_Aligned.sortedByCoord.out.bam.bai"
        ]
        for f in corrupted_files:
            if os.path.exists(f):
                size = os.path.getsize(f)
                os.remove(f)
                print(f"  Removed {f} ({size / 1024**3:.2f} GB)")

    # 2. Clean up intermediate files to free space
    print("\nCleaning up intermediate files...")
    patterns_to_clean = [
        '*_Unmapped.out.mate*',
        '*_SJ.out.tab',
        '*_Log.progress.out'
    ]

    for pattern in patterns_to_clean:
        files = glob.glob(f"{config['output_dir']}/{pattern}")
        for f in files:
            size = os.path.getsize(f)
            os.remove(f)
            print(f"  Removed {f} ({size / 1024**2:.1f} MB)")

    # 3. Check available space
    available = shutil.disk_usage(config['output_dir']).free
    required = len(sample_ids) * 50 * 1024**3  # 50 GB per sample

    print(f"\nDisk space check:")
    print(f"  Available: {available / 1024**3:.1f} GB")
    print(f"  Required: {required / 1024**3:.1f} GB")

    if available < required:
        raise ValueError(f"Still insufficient disk space: {available / 1024**3:.1f} GB available, "
                        f"{required / 1024**3:.1f} GB required")

    # 4. Update checkpoint to re-run from alignment
    with open('pipeline_checkpoint.json') as f:
        checkpoint = json.load(f)

    for sample_id in sample_ids:
        # Keep only 'fastqc' and 'trim', remove 'align' and 'count'
        checkpoint[sample_id] = [s for s in checkpoint[sample_id]
                                if s in ['fastqc', 'trim']]

    with open('pipeline_checkpoint.json', 'w') as f:
        json.dump(checkpoint, f, indent=2)

    print(f"\n✓ Ready to re-run from alignment for {len(sample_ids)} samples")

# Apply fix
fix_disk_space_and_rerun(['sample_037', 'sample_039'], config)

# Re-run pipeline
run_pipeline_with_checkpoint(['sample_037', 'sample_039'], config)
```

## Problem 3: Low Alignment Rate (sample_039)

### Diagnosis Process

**Detailed quality check:**

```python
def diagnose_low_alignment(sample_id, config):
    """Detailed diagnosis of low alignment rate"""

    print(f"Diagnosing low alignment for {sample_id}...")

    # 1. Check raw read quality
    fastqc_data = f"{config['output_dir']}/fastqc_raw/{sample_id}_R1_fastqc/fastqc_data.txt"

    with open(fastqc_data) as f:
        content = f.read()
        if 'FAIL' in content:
            print("  ⚠ Raw read quality FAILED FastQC")

    # 2. Check trimming statistics
    trim_log = f"{config['output_dir']}/{sample_id}_trimmomatic.log"
    with open(trim_log) as f:
        log = f.read()
        survival_match = re.search(r'Both Surviving: (\d+) \(([\d.]+)%\)', log)
        if survival_match:
            survival_pct = float(survival_match.group(2))
            print(f"  Trimming survival rate: {survival_pct:.1f}%")
            if survival_pct < 70:
                print("  ⚠ High read loss during trimming")

    # 3. Check alignment details
    align_log = f"{config['output_dir']}/{sample_id}_Log.final.out"
    with open(align_log) as f:
        log = f.read()

        # Parse detailed unmapping reasons
        unmapped_reasons = {}
        for line in log.split('\n'):
            if '% of reads unmapped: too short' in line:
                unmapped_reasons['too_short'] = float(line.split('|')[1].strip().rstrip('%'))
            elif '% of reads unmapped: too many mismatches' in line:
                unmapped_reasons['mismatches'] = float(line.split('|')[1].strip().rstrip('%'))
            elif '% of reads unmapped: other' in line:
                unmapped_reasons['other'] = float(line.split('|')[1].strip().rstrip('%'))

        print(f"  Unmapped read breakdown:")
        for reason, pct in unmapped_reasons.items():
            print(f"    {reason}: {pct:.1f}%")

    # 4. Sample reads for contamination screen
    subprocess.run([
        'fastq-screen',
        '--aligner', 'bowtie2',
        '--conf', 'fastq_screen.conf',
        '--subset', '100000',
        '--outdir', f"{config['output_dir']}/fastq_screen",
        f"{config['output_dir']}/{sample_id}_R1_trimmed.fastq.gz"
    ])

    # Parse FastQ Screen results
    screen_file = f"{config['output_dir']}/fastq_screen/{sample_id}_R1_trimmed_screen.txt"
    print(f"\n  Contamination screening results:")
    with open(screen_file) as f:
        for line in f:
            if not line.startswith('#') and line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    species = parts[0]
                    pct_mapped = parts[1]
                    print(f"    {species}: {pct_mapped}% mapped")

# Run diagnosis
diagnose_low_alignment('sample_039', config)

# Output:
# Diagnosing low alignment for sample_039...
#   ⚠ Raw read quality FAILED FastQC
#   Trimming survival rate: 58.3%
#   ⚠ High read loss during trimming
#   Unmapped read breakdown:
#     too_short: 32.1%
#     mismatches: 12.4%
#     other: 10.3%
#
#   Contamination screening results:
#     Human: 45.2% mapped
#     rRNA: 28.3% mapped
#     Adapters: 15.2% mapped
#     PhiX: 0.3% mapped
```

**Root cause:** Sample has poor quality, high rRNA contamination, and residual adapters.

### Solution

This sample requires reprocessing from raw FASTQ with adjusted parameters:

```python
def reprocess_poor_quality_sample(sample_id, config):
    """Reprocess sample with more aggressive QC"""

    print(f"Reprocessing {sample_id} with enhanced QC...")

    # 1. More aggressive trimming
    trimmomatic_cmd = [
        'trimmomatic', 'PE',
        '-threads', '8',
        '-phred33',
        f"{config['fastq_dir']}/{sample_id}_R1.fastq.gz",
        f"{config['fastq_dir']}/{sample_id}_R2.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R1_trimmed.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R1_unpaired.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R2_trimmed.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R2_unpaired.fastq.gz",
        'ILLUMINACLIP:adapters.fa:2:30:10',
        'LEADING:5',  # More aggressive than default (3)
        'TRAILING:5',  # More aggressive than default (3)
        'SLIDINGWINDOW:4:20',  # Higher quality threshold (was 15)
        'MINLEN:50'  # Longer minimum length (was 36)
    ]

    subprocess.run(trimmomatic_cmd, check=True)

    # 2. Filter rRNA contamination with SortMeRNA
    sortmerna_cmd = [
        'sortmerna',
        '--ref', 'refs/rRNA_databases/silva-bac-16s-id90.fasta',
        '--ref', 'refs/rRNA_databases/silva-bac-23s-id98.fasta',
        '--ref', 'refs/rRNA_databases/silva-euk-18s-id95.fasta',
        '--ref', 'refs/rRNA_databases/silva-euk-28s-id98.fasta',
        '--reads', f"{config['output_dir']}/{sample_id}_R1_trimmed.fastq.gz",
        '--reads', f"{config['output_dir']}/{sample_id}_R2_trimmed.fastq.gz",
        '--aligned', f"{config['output_dir']}/{sample_id}_rRNA",
        '--other', f"{config['output_dir']}/{sample_id}_clean",
        '--fastx',
        '--paired_in',
        '--threads', '8'
    ]

    subprocess.run(sortmerna_cmd, check=True)

    # Rename clean files for downstream processing
    os.rename(
        f"{config['output_dir']}/{sample_id}_clean_fwd.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R1_trimmed.fastq.gz"
    )
    os.rename(
        f"{config['output_dir']}/{sample_id}_clean_rev.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R2_trimmed.fastq.gz"
    )

    # 3. Check if enough reads remain
    def count_fastq_reads(fastq_file):
        with gzip.open(fastq_file, 'rt') as f:
            return sum(1 for _ in f) // 4

    remaining_reads = count_fastq_reads(
        f"{config['output_dir']}/{sample_id}_R1_trimmed.fastq.gz"
    )

    print(f"  Reads after enhanced QC: {remaining_reads:,}")

    if remaining_reads < 1_000_000:
        print(f"  ⚠ WARNING: Very few reads remaining ({remaining_reads:,})")
        print(f"  Consider excluding {sample_id} from analysis")
        return False

    # 4. Re-run alignment with adjusted parameters
    star_cmd = [
        'STAR',
        '--genomeDir', config['genome_index'],
        '--readFilesIn',
        f"{config['output_dir']}/{sample_id}_R1_trimmed.fastq.gz",
        f"{config['output_dir']}/{sample_id}_R2_trimmed.fastq.gz",
        '--readFilesCommand', 'zcat',
        '--outFileNamePrefix', f"{config['output_dir']}/{sample_id}_",
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--outFilterMultimapNmax', '20',  # Allow more multimappers
        '--outFilterMismatchNmax', '999',  # Don't filter by mismatches
        '--outFilterMismatchNoverReadLmax', '0.04',  # Allow 4% mismatches
        '--alignIntronMin', '20',
        '--alignIntronMax', '1000000',
        '--runThreadN', '8'
    ]

    subprocess.run(star_cmd, check=True)

    print(f"✓ Reprocessing complete for {sample_id}")
    return True

# Reprocess problematic sample
success = reprocess_poor_quality_sample('sample_039', config)

if success:
    # Continue with quantification
    run_quantification(['sample_039'], config)
else:
    # Mark sample as excluded
    mark_sample_excluded('sample_039', 'insufficient_reads_after_qc')
```

## Systematic Debugging Workflow

### Create debugging toolkit:

```python
class PipelineDebugger:
    """Comprehensive pipeline debugging toolkit"""

    def __init__(self, config):
        self.config = config

    def diagnose_all_failures(self):
        """Diagnose all failed samples systematically"""

        # Load checkpoint
        with open('pipeline_checkpoint.json') as f:
            checkpoint = json.load(f)

        # Identify failed samples
        stages = ['fastqc', 'trim', 'align', 'count']
        failed_samples = {}

        for sample_id, completed in checkpoint.items():
            if len(completed) < len(stages):
                failed_stage = stages[len(completed)]
                if failed_stage not in failed_samples:
                    failed_samples[failed_stage] = []
                failed_samples[failed_stage].append(sample_id)

        # Generate diagnostic report
        print("="*60)
        print("PIPELINE FAILURE DIAGNOSTIC REPORT")
        print("="*60)

        for stage, samples in failed_samples.items():
            print(f"\n{stage.upper()} failures: {len(samples)} samples")
            print(f"  Samples: {', '.join(samples)}")

            # Stage-specific diagnosis
            if stage == 'align':
                self.diagnose_alignment_failures(samples)
            elif stage == 'count':
                self.diagnose_quantification_failures(samples)

    def diagnose_alignment_failures(self, sample_ids):
        """Diagnose alignment failures for multiple samples"""

        print("\n  Running alignment diagnostics...")

        for sample_id in sample_ids:
            log_file = f"{self.config['output_dir']}/{sample_id}_Log.final.out"

            if not os.path.exists(log_file):
                print(f"    {sample_id}: No log file (alignment did not start)")
                continue

            with open(log_file) as f:
                log = f.read()
                mapped_pct = re.search(r'Uniquely mapped reads %.*\|([\d.]+)%', log)

                if mapped_pct:
                    pct = float(mapped_pct.group(1))
                    if pct < 10:
                        print(f"    {sample_id}: Severe mapping failure ({pct:.1f}%) "
                              "- likely wrong reference or contamination")
                    elif pct < 50:
                        print(f"    {sample_id}: Poor mapping ({pct:.1f}%) "
                              "- check quality and reference")

    def diagnose_quantification_failures(self, sample_ids):
        """Diagnose quantification failures"""

        print("\n  Running quantification diagnostics...")

        for sample_id in sample_ids:
            bam_file = f"{self.config['output_dir']}/{sample_id}_Aligned.sortedByCoord.out.bam"

            if not os.path.exists(bam_file):
                print(f"    {sample_id}: BAM file missing (alignment failed)")
            else:
                # Check BAM validity
                result = subprocess.run(
                    ['samtools', 'quickcheck', bam_file],
                    capture_output=True
                )

                if result.returncode != 0:
                    print(f"    {sample_id}: BAM file corrupted - likely disk space issue")
                else:
                    print(f"    {sample_id}: BAM valid - check featureCounts logs")

    def generate_recovery_plan(self):
        """Generate step-by-step recovery plan"""

        print("\n" + "="*60)
        print("RECOVERY PLAN")
        print("="*60)

        # This would analyze all failures and generate specific commands to fix them
        print("\n1. Fix species mismatches (samples: 042, 043)")
        print("   python fix_species.py --samples sample_042 sample_043 --species mouse")

        print("\n2. Clean disk space and retry corrupted BAMs (samples: 037, 039)")
        print("   python clean_and_retry.py --samples sample_037 sample_039")

        print("\n3. Reprocess low-quality sample with enhanced QC (sample: 039)")
        print("   python enhanced_qc.py --sample sample_039")

        print("\n4. Resume pipeline from checkpoint")
        print("   python run_pipeline.py --resume")

# Run debugger
debugger = PipelineDebugger(config)
debugger.diagnose_all_failures()
debugger.generate_recovery_plan()
```

## Summary

This debugging example demonstrated:

1. **Systematic diagnosis** - Checking logs, checkpoints, and metrics to identify root causes
2. **Problem classification** - Identifying different failure types (wrong reference, disk space, poor quality)
3. **Targeted solutions** - Applying specific fixes for each problem type
4. **Recovery strategy** - Using checkpoints to resume processing without repeating successful stages
5. **Prevention** - Adding validation checks to catch similar issues earlier in future runs

Key debugging principles:
- Always check logs and metrics before making assumptions
- Use contamination screening to identify species mismatches
- Verify file integrity (not just existence)
- Monitor resource usage (disk space, memory)
- Implement checkpointing for efficient recovery
- Create systematic debugging workflows for reproducibility
