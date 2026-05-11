# Upstream VNtyper Inventory

This inventory maps the upstream VNtyper implementation to the BioScript port.
It is the behavioral reference for the first BioScript milestone:

```text
BAM input -> MUC1 read extraction -> Kestrel VCF -> frameshift/depth
classification -> TSV/JSON report
```

## Source Paths Read

- `ports/vntyper/vntyper/vntyper/cli.py`
  Defines CLI arguments, input selection, reference assembly options, output
  paths, report generation, and module toggles.
- `ports/vntyper/vntyper/vntyper/scripts/pipeline.py`
  Orchestrates validation, output directories, input-type routing, BAM/CRAM or
  FASTQ preprocessing, Kestrel execution, summary files, reports, and optional
  modules.
- `ports/vntyper/vntyper/vntyper/scripts/fastq_bam_processing.py`
  Handles FASTQ QC, BAM/CRAM region slicing, unmapped-read retention, FASTQ
  extraction, coverage calculation, BAM header parsing, and assembly detection.
- `ports/vntyper/vntyper/vntyper/scripts/kestrel_genotyping.py`
  Builds the Kestrel Java command, runs Kestrel for configured k-mer sizes,
  converts Kestrel SAM to BAM, filters VCFs to indels, compresses with
  bcftools when available, splits insertion/deletion VCFs, processes k-mer
  results, flags variants, selects the best variant, and writes
  `kestrel_result.tsv`.
- `ports/vntyper/vntyper/vntyper/scripts/file_processing.py`
  Filters VCF rows to indels and splits indel VCFs into insertion/deletion
  VCFs.
- `ports/vntyper/vntyper/vntyper/scripts/variant_parsing.py`
  Reads VCF records into table rows and applies final ALT-based filtering.
- `ports/vntyper/vntyper/vntyper/scripts/motif_processing.py`
  Loads the MUC1 motif reference, preprocesses insertion/deletion rows, applies
  motif correction and annotation, and deduplicates frame-shift candidates.
- `ports/vntyper/vntyper/vntyper/scripts/scoring.py`
  Splits the Kestrel sample field into alternate/active-region depths,
  computes frame score, derives direction and frame-shift amount, and marks
  valid insertion/deletion frame-shift patterns.
- `ports/vntyper/vntyper/vntyper/scripts/confidence_assignment.py`
  Computes depth score and assigns `Negative`, `Low_Precision`,
  `High_Precision`, or `High_Precision*` from Kestrel config thresholds.
- `ports/vntyper/vntyper/vntyper/scripts/flagging.py`
  Applies configured row-level flag rules and duplicate detection before final
  variant selection.
- `ports/vntyper/vntyper/vntyper/scripts/region_utils.py`
  Resolves assembly aliases, detects chromosome naming from BAM headers, and
  builds MUC1 region strings.
- `ports/vntyper/vntyper/vntyper/scripts/reference_registry.py`
  Defines canonical assembly names, coordinate systems, reference sources, MUC1
  coordinate ranges, chromosome naming, and registry validation.
- `ports/vntyper/vntyper/vntyper/scripts/chromosome_utils.py`
  Detects assembly and chromosome naming from contigs and validates chromosome
  names for UCSC, NCBI, and Ensembl styles.
- `ports/vntyper/vntyper/vntyper/scripts/generate_report.py`
  Builds screening summaries, loads fastp/log/summary data, renders HTML, and
  optionally adds IGV content.
- `ports/vntyper/vntyper/vntyper/scripts/kestrel_config.json`
  Provides Kestrel, frame-score, depth-confidence, ALT-filtering, motif, and
  flagging thresholds.
- `ports/vntyper/vntyper/vntyper/scripts/report_config.json`
  Provides Kestrel/adVNTR screening summary decision rules.

## Minimal Pipeline Surface

The first BioScript port should keep the optional module surface out of the
critical path and implement this narrow path first:

1. Validate one input mode: BAM first, FASTQ later.
2. Resolve MUC1 broad BAM region and VNTR coverage region for the selected
   assembly and chromosome naming convention.
3. Build safe external-tool argv for `samtools view`, `samtools index`,
   `samtools fastq`, and `samtools depth`.
4. Build safe external-tool argv for Kestrel with VNtyper defaults:
   k-mer `20`, Java memory `12g`, max align states `40`, max hap states `40`,
   SAM haplotype output, stdout/stderr logging, and temporary directory.
5. Parse Kestrel VCF rows, filter to indels, split insertion/deletion records,
   and normalize sample-depth fields.
6. Compute frame score, direction, frame-shift amount, valid frame-shift flag,
   depth score, confidence, ALT filters, motif annotations, row flags, and
   final best-variant selection.
7. Emit deterministic `kestrel_result.tsv` and structured JSON before HTML.

## Current BioScript Coverage

Already implemented:

- `ports/vntyper/bioscript/vntyper_regions.py` for assembly aliases,
  coordinate lookup, chromosome naming, naming-convention detection, and region
  string construction.
- `ports/vntyper/bioscript/vntyper_commands.py` for deterministic BAM-path
  command planning across region slicing, indexing, FASTQ extraction, coverage,
  Kestrel, and bcftools post-processing.
- `bioscript.samtools` command builders for `view_region`, `fastq`, `depth`,
  and `index`.
- `bioscript.bcftools` command builders for `sort`, `index`, `view_filter`,
  and `norm`.
- `bioscript.kestrel.build_command` matching the VNtyper Kestrel defaults.
- `bioscript.vcf.read_kestrel` for Kestrel VCF rows.
- `ports/vntyper/bioscript/vntyper_port.py` for Kestrel VCF parsing,
  frame/depth/confidence post-processing, and report JSON from fixture rows.

Still missing for parity:

- BAM-header-aware chromosome naming detection.
- Full `process_bam_to_fastq` command plan including unmapped-read retention.
- Kestrel post-processing parity for motif annotation, duplicate flagging, and
  final best-variant selection.
- Coverage QC parsing from `samtools depth`.
- Deterministic TSV parity against upstream `kestrel_result.tsv` and
  `kestrel_pre_result.tsv`.
- HTML report parity.

## Upstream Outputs To Match

The core parity checks should compare:

- `kestrel/output.vcf`
- `kestrel/output_indel.vcf`
- `kestrel/output_insertion.vcf`
- `kestrel/output_deletion.vcf`
- `kestrel/kestrel_pre_result.tsv`
- `kestrel/kestrel_result.tsv`
- pipeline summary JSON
- coverage summary TSV
- selected report JSON fields used by the first HTML report

The copied large data currently provides BAM/BAI and FASTQ inputs, but not the
expected VCF/TSV/JSON outputs, so those still need to be generated from
upstream VNtyper or added as fixtures.
