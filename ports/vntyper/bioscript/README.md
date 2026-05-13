# VNtyper BioScript Port

This directory will contain the BioScript implementation of the VNtyper
pipeline. Keep code here focused on VNtyper behavior; reusable compatibility
layers should live in `rust/bioscript-libs` and be exposed through
`from bioscript import ...` modules.

## Target Interface

The user-facing BioScript program path is:

```text
ports/vntyper/bioscript/vntyper.bs
```

`vntyper.bs.py` remains an executable sketch until the runtime can execute the
same flow as real BioScript syntax.

The port should expose two entry points.

### BAM Input

```python
run_vntyper(
    bam=input_file,
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
    reference_build="hg19",
    output_dir=output_dir,
    participant_id=participant_id,
)
```

Expected native flow:

```text
FASTQ pair -> bioscript.kestrel.run_native
           -> bioscript.bcftools.sort_native/index_native
           -> VNtyper TSV/JSON/HTML report logic
```

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

FASTQ native parity is still a target gate and should be added next to the BAM
gate once Kestrel FASTQ fixture expectations are locked down.
