# VNtyper Input Modes

## Current BioScript Milestone

The current BioScript port targets prealigned BAM input first.

The command planner covers:

- selecting the MUC1 broad BAM region,
- slicing the BAM with `samtools view`,
- indexing the sliced BAM,
- converting the slice to paired FASTQ with `samtools fastq`,
- calculating VNTR coverage with `samtools depth`,
- calling Kestrel over the extracted reads,
- sorting/indexing the Kestrel VCF with `bcftools`.

## FASTQ Input

FASTQ input is deferred. Upstream VNtyper can run fastp and BWA before the
Kestrel path, but BioScript does not need `bioscript.fastp` or `bioscript.bwa`
for the first BAM milestone.

When FASTQ support is reopened:

- add `bioscript.fastp` command builders for QC/trimming,
- add `bioscript.bwa` command builders for paired-end alignment,
- add FASTQ integration fixtures and expected BAM/Kestrel outputs,
- decide whether the first public FASTQ API accepts raw FASTQs or requires a
  preconfigured reference index.

Until then, users should provide prealigned BAM/BAI inputs.
