# `bioscript.kestrel` API Sketch

VNtyper uses Kestrel as a mapping-free MUC1-VNTR caller. The first BioScript
surface should be Python-shaped and structured, even if the first backend still
executes the Java Kestrel release.

Import form:

```python
from bioscript import kestrel
```

## Initial API

```python
result = kestrel.run(
    fastq_1="sample_R1.fastq.gz",
    fastq_2="sample_R2.fastq.gz",
    reference_vntr="All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
    output_vcf="output.vcf",
    output_sam="output.sam",
    temp_dir="work",
    sample_name="sample",
    kmer_size=20,
    max_align_states=40,
    max_hap_states=40,
    memory="12g",
    log_level="INFO",
)
```

Return shape:

```python
{
    "vcf": "output.vcf",
    "sam": "output.sam",
    "kmer_size": 20,
    "sample_name": "sample",
    "records": kestrel.read_vcf("output.vcf"),
}
```

## Command Builder

`kestrel.build_command(...)` should exist for tests, but it should return a
structured argv list, not a shell string:

```python
[
    "java",
    "-Xmx12g",
    "-jar",
    "kestrel.jar",
    "-k",
    "20",
    "--maxalignstates",
    "40",
    "--maxhapstates",
    "40",
    "-r",
    reference_vntr,
    "-o",
    output_vcf,
    "-sSAMPLE",
    fastq_1,
    fastq_2,
    "--hapfmt",
    "sam",
    "-p",
    output_sam,
    "--logstderr",
    "--logstdout",
    "--loglevel",
    "INFO",
    "--temploc",
    temp_dir,
]
```

This mirrors the exact Kestrel options VNtyper currently constructs in
`vntyper/scripts/kestrel_genotyping.py`.

## Backend Plan

1. `java` backend:
   Run a configured Kestrel JAR/release with safe argv construction. This is
   the first parity target.
2. `rust` backend:
   Port only the Kestrel internals VNtyper needs. Candidate Java packages:
   `counter`, `activeregion`, `align`, `variant`, and `writer.vcf`.
3. `auto` backend:
   Use Rust when feature-complete for the requested options, otherwise fall
   back to the Java adapter if allowed by runtime policy.

## VNtyper-Specific Defaults

```python
{
    "kmer_size": 20,
    "max_align_states": 40,
    "max_hap_states": 40,
    "memory": "12g",
    "additional_settings": "",
}
```

