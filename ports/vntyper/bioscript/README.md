# VNtyper BioScript Port

This directory will contain the BioScript implementation of the VNtyper
pipeline. Keep code here focused on VNtyper behavior; reusable compatibility
layers should live in `rust/bioscript-libs` and be exposed through
`from bioscript import ...` modules.

## Target Interface

The user-facing BioScript program paths are:

```text
ports/vntyper/bioscript/vntyper.bs
ports/vntyper/bioscript/vntyper-bam-native.bs
ports/vntyper/bioscript/vntyper-fastq.bs
```

`vntyper.bs.py` remains an executable sketch and Python oracle support module;
the runtime-facing BAM flow now lives in `vntyper.bs`.

The port should expose two entry points.

### BAM Input

```python
run_vntyper(
    bam=input_file,
    bai=input_bai,
    bam_region=bam_region,
    vntr_region=vntr_region,
    reference_fasta=reference_fasta,
    kmer_size=20,
    minimum_difference=5,
    max_haplotypes=2,
    max_saved_states=2,
    reference_build="hg19",
    output_dir=output_dir,
    participant_id=participant_id,
)
```

Expected native flow:

```text
BAM -> bioscript.samtools.view_region_native
    -> bioscript.samtools.fastq_native
    -> bioscript.samtools.depth_native
    -> bioscript.kestrel.run_native
    -> bioscript.bcftools.sort_native/index_native
    -> VNtyper TSV/JSON/HTML report logic
```

### FASTQ Input

```python
run_vntyper_fastq(
    r1=fastq_1,
    r2=fastq_2,
    reference_fasta=reference_fasta,
    kmer_size=20,
    minimum_difference=5,
    max_haplotypes=2,
    max_saved_states=2,
    reference_build="hg19",
    output_dir=output_dir,
    participant_id=participant_id,
)
```

Expected native flow:

```text
FASTQ pair -> bioscript.kestrel.run_native
           -> bioscript.bcftools.sort/index
           -> bioscript.vcf.read_vntyper_kestrel
           -> TSV execution summary
```

`vntyper-fastq.bs` currently exercises this native BioScript runtime path on
tiny deterministic fixtures. Core Kestrel call-table conversion is now in
`bioscript-libs`; the FASTQ runtime slice also builds the VNtyper JSON report
payload through the VCF facade and writes `kestrel_result.tsv` from the same
rows. Full HTML report generation still lives in the Python scaffold.

`vntyper.bs` and `vntyper-bam-native.bs` exercise the BAM native runtime path
on caller-provided regions and reference paths. They run Samtools, Kestrel,
BCFtools, and VCF facades through the BioScript runtime and materialize
BAM-slice, FASTQ, VCF, TSV, and report JSON artifacts.

The tiny runtime tests pass small Kestrel values so short synthetic references
can produce calls. Current native parity gates use `kmer_size=20`,
`minimum_difference=5`, `max_haplotypes=2`, and `max_saved_states=2` as bounded
native-engine settings. Upstream Java VNtyper config uses `max_align_states=40`
and `max_hap_states=40`; matching those limits in `kestrel-rs` is still part of
the Kestrel parity/performance gap.

## Local Test Gates

Small VNtyper-port tests:

```sh
PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'
```

Opt-in large BAM parity:

```sh
BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_native_bam_pipeline_gate
```

Opt-in large FASTQ parity:

```sh
BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_native_fastq_pipeline_gate
```

This gate runs native Kestrel and native BCFtools against representative
positive and negative FASTQ fixtures, then compares the generated classification
and report shape to expected VNtyper outputs.
