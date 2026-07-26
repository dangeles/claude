---
name: data-pipeline-manager
version: 1.0
last_updated: 2026-06-09
description: Use when designing or troubleshooting robust data pipelines — building stage-by-stage quality validation, error handling, retry/checkpointing, and monitoring for ETL or bioinformatics data flows (e.g., FASTQ→BAM→counts). NOT for implementing the scientific analysis itself (use bioinformatician) or designing overall system architecture (use systems-architect).
success_criteria:
  - Pipeline runs end-to-end without errors
  - Quality validation implemented at each stage
  - Error handling comprehensive and informative
  - Data format mismatches caught and reported
  - Monitoring and logging in place for troubleshooting
  - Pipeline failure modes documented with recovery strategies
  - Performance acceptable for expected data volumes
category: bioinformatics-workflow
integrates-with:
  - bioinformatician
  - systems-architect
---

# Data Pipeline Manager

## Purpose

Guidance for designing, implementing, and troubleshooting data processing pipelines, with
emphasis on quality validation, error handling, and reliability. Pipelines are a prevalent
source of failure in automated systems: poor data quality, format mismatches, missing
files, and transient errors cascade through stages and become hard to diagnose after the
fact.

This skill covers the layer between the bioinformatician skill, which implements the
specific analysis, and the systems-architect skill, which designs the overall system:
pipeline orchestration, validation, and reliability.

## When to Use

- **Design and implementation**: new pipelines from scratch; RNA-seq, ChIP-seq, or other
  genomics processing; ETL workflows for data integration; ML preprocessing; batch or
  streaming systems.
- **Debugging**: tracing errors through multi-stage workflows, diagnosing format
  mismatches between stages, finding bottlenecks, resolving dependency and ordering
  issues.
- **Quality assurance**: validation at each stage, format consistency between stages,
  integrity checks and checksums, quality metrics tracked over time.
- **Error recovery**: retry for transient failures, checkpoint and resume, partial
  failures in parallel processing, rollback, idempotent operations.
- **Optimization**: parallelizing independent stages, reducing inter-stage data transfer,
  caching expensive operations, scaling to larger volumes.

When a pipeline investigation spans several genuinely independent tracks — unrelated
stages, separate sample batches — fanning the work out can pay off; see
`../references/delegation-and-scope.md` for when the overhead is worth it.

## Pipeline Design

Start by pinning down the contract: input sources with their formats, schemas, and
expected volume; each transformation stage's purpose and dependencies; output formats,
destinations, and success criteria; resource estimates per stage; and the time, cost, and
resource constraints you are working inside. The dependency map tells you which stages are
sequential and which can run in parallel, which in turn selects the pattern below.

Build in observability from the start rather than retrofitting it: metrics per stage
(records processed, time taken, errors), a logging strategy (what, where, at what level),
alerting thresholds, and tracing for distributed components.

## Pipeline Patterns

**Sequential** — each stage consumes the previous stage's output. RNA-seq (FASTQ → QC →
alignment → count matrix → differential expression), classic ETL, image processing.
Simple to reason about, with natural checkpointing at each stage; the slowest stage is the
bottleneck and any failure stops the run. Give each stage a defined input/output contract,
validate between stages, and keep intermediate outputs for debugging and resume.

**Parallel** — independent branches with no dependencies: many samples, several tools over
the same input, partitioned data, fan-out of work to workers. High throughput and isolated
failures, at the cost of resource management and result aggregation. Partition data
deliberately, handle partial failures, and balance load to avoid stragglers.

**Conditional** — branch on data content or intermediate results: route by data type or
quality, skip expensive steps that aren't needed. Efficient and adaptive, but every path
needs testing, a consistent output format, and a log line recording which path each item
took.

**Hybrid** — a mix, for multi-omics integration, complex bioinformatics workflows, or
enterprise ETL across several sources. Most flexible and hardest to debug; document the
DAG explicitly and modularize components.

For anything beyond a simple sequence, use a workflow manager (Snakemake, Nextflow,
Airflow) rather than hand-rolled orchestration.

## Input Validation

Validate inputs before processing rather than discovering problems mid-run. Structural
checks: existence and accessibility, format against its specification (FASTQ, BAM, CSV,
JSON), completeness (checksums, file size), schema (columns and data types), non-empty
content, text encoding, read/write permissions on every path used.

Content checks go further: numeric values within expected ranges, required fields present
and non-null, referential integrity across identifiers, internal consistency between
related fields (start < end coordinates), domain quality metrics (quality scores, read
lengths), and sample identifiers matching the expected manifest.

Fail fast on validation errors. Report the specific file, line, and field; say what is
wrong and how to fix it; stop before consuming expensive compute; log the error so trends
are visible across runs.

Concrete validators — FASTQ, CSV, BAM/SAM, DataFrame schema, gene-count matrix, genome
coordinates, disk space, and cross-stage sample and dimension consistency — are in
`references/validation-patterns.md`.

## Transform

Each stage runs with a clear input/output contract and its own error handling. Practices
that make transformations survivable:

- **Idempotency**: operations safe to retry without side effects.
- **Atomic writes**: write to a temp path, then rename.
- **Batching and streaming**: process in memory-sized batches; stream rather than loading
  whole datasets; chunk or window large files.
- **Error context**: capture the input data, parameters, and environment alongside the
  exception.
- **Graceful degradation**: keep processing remaining records when an individual record
  fails.
- **Checkpointing**: persist intermediate results so a failed run resumes rather than
  restarts.
- **Progress and resource logging**: log progress every N records; track memory, CPU, and
  disk to catch resource exhaustion before it becomes an opaque crash.

For large-scale data, add backpressure so a fast stage cannot overwhelm a slow downstream
one, use efficient intermediate formats (Parquet, Arrow), and consider a distributed
framework (Spark, Dask) when a single machine stops being enough.

## Output Validation

Structure: correct format, expected schema, expected dimensions for tabular or array data,
all expected files present, file sizes plausible (not empty, not suspiciously large).

Content: values within biologically or logically plausible ranges; no unexpected NaN, Inf,
or null; distributions reasonable rather than all-zero or degenerate; logical relationships
between output fields; domain-specific sanity checks (gene counts non-negative); key
metrics compared against a baseline or previous run.

Report quality metrics with the outputs: record counts (input vs output, filtered, passed),
success and error rates by stage, data quality scores and flags, processing time and
resource usage, coverage and completeness, and domain metrics such as alignment and
duplication rates.

## Error Handling

Classify the failure first, since the class determines whether retrying can possibly help:
transient (network timeouts, temporary unavailability, rate limits), permanent (format
errors, missing required data, logic errors), configuration (wrong parameters, missing
credentials, invalid paths), resource (out of memory, disk full, quota exceeded), and data
(invalid values, corrupted files, schema mismatches).

**Retry** applies to transient failures: exponential backoff with jitter (1s, 2s, 4s, 8s),
3-5 retries as a typical ceiling, per-operation timeouts, a circuit breaker that stops
retrying once the error rate crosses a threshold, and idempotent operations so a retry
cannot double-process.

**Logging** for diagnosis: structured format (JSON) so logs can be aggregated; include
input data, parameters, environment, and stack traces; use appropriate severity levels;
carry correlation IDs across distributed components; keep credentials and PII out.

**Checkpointing** for resumability: save state after expensive operations, record which
stages completed and what they produced along with what remains, resume from the last good
checkpoint, clean up checkpoints on success, and handle concurrent checkpoint access in
distributed runs.

Implementations — exponential backoff, conditional and async retry, file- and
database-backed checkpoints, structured logging, circuit breaker, fallback and
partial-success patterns, automatic recovery — are in
`references/error-handling-strategies.md`.

## Monitoring and Alerting

Track while a pipeline runs: progress through each stage, processing rate
(records/bytes per second), input queue depth and backlog, error rate per stage and
overall, CPU/memory/disk/network usage, and end-to-end plus per-stage latency.

Alert on failures immediately; on performance degrading below threshold; on quality
metrics leaving expected ranges; on approaching resource limits; and on SLAs at risk.
Aggregate and prioritize so alerts stay meaningful.

A useful dashboard shows overall health, stage-by-stage progress, error rates and types
over time, resource utilization trends, data quality metrics, and comparison against
historical performance.

## Bioinformatics-Specific Considerations

### FASTQ to Count Matrix Pipeline

**Standard workflow**:
1. **Quality Control**: FastQC on raw FASTQ files
2. **Preprocessing**: Trimming adapters and low-quality bases (Trimmomatic, cutadapt)
3. **Alignment**: Align to reference genome (STAR, HISAT2, BWA)
4. **Quality Control**: Check alignment metrics (samtools, picard)
5. **Quantification**: Count reads per gene (featureCounts, HTSeq, RSEM)
6. **Count Matrix**: Assemble counts into matrix for downstream analysis

**Critical validation points**:
- Raw read quality scores (median Q30+)
- Adapter contamination levels
- Alignment rate (>70% for RNA-seq)
- Duplicate rate (<30%)
- Gene detection rate (>10,000 genes for human)
- Count distribution (not zero-inflated; mean count >100, coefficient of variation <10)

### Genome Build Consistency

Mixing genome builds (hg19 vs hg38, mm9 vs mm10) causes misalignments and silently
incorrect results. Document the build in metadata and file headers, validate build
consistency across all samples, check the reference FASTA matches the annotation GTF,
verify chromosome naming conventions (chr1 vs 1), pass the build as an explicit pipeline
parameter, and include it in output file names.

### Sample Metadata Tracking

Essential fields: sample identifier (unique across the study), biological
condition/phenotype, technical batch, sequencing run, library preparation protocol,
genome build and annotation version, processing date and pipeline version.

Validate that sample IDs match between FASTQ files and metadata, that no sample ID is
duplicated, that required fields are present, that controlled-vocabulary terms are valid,
and that metadata agrees with the data (species matches the reference).

### Reference Data Management

Version-control reference genomes and annotations, keep naming consistent, store
checksums, document source and version, maintain separate references per genome build, and
check reference integrity before a run.

## Integration with Other Skills

**bioinformatician** implements the analysis inside the framework this skill provides: it
selects tools and debugs analysis-specific issues and algorithm performance, while this
skill supplies the pipeline architecture, validation, error handling, and orchestration
around it.

**systems-architect** designs the surrounding system and infrastructure — scaling,
infrastructure reliability, monitoring platforms — while this skill designs the pipeline
components, parallelization, reliability, and pipeline-specific monitoring that sit on top.

## References

- `references/validation-patterns.md` — validators for filesystem, formats, schemas, data
  quality, bioinformatics-specific and cross-stage checks
- `references/error-handling-strategies.md` — retry, checkpointing, logging, circuit
  breakers, graceful degradation, recovery
- `../references/delegation-and-scope.md` — when to fan work out to subagents

## Examples

- `examples/rnaseq-pipeline-example.md` — full FASTQ → counts pipeline with validation
- `examples/pipeline-debugging-example.md` — tracing a failed multi-stage run
