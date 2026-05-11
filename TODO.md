# BioScript VNtyper Port TODO

This tracks the VNtyper port from <https://github.com/hassansaei/VNtyper> into
BioScript. The goal is not to rewrite every upstream dependency immediately.
The goal is to make VNtyper run in BioScript while extending BioScript,
`bioscript-libs`, and native/external tool wrappers only where the VNtyper
surface requires it.

## Directory Layout

- [x] `ports/vntyper/vntyper`
      Upstream VNtyper source, added as a git submodule for source reading and
      parity checks.
- [x] `ports/vntyper/kestrel`
      Upstream Kestrel Java source, added as a git submodule for source reading
      and eventual native porting.
- [x] `ports/vntyper/bioscript`
      BioScript implementation of the VNtyper pipeline and small ported helper
      modules.
- [x] `ports/vntyper/tests`
      BioScript-owned tests that compare the port against upstream behavior.
- [x] `ports/vntyper/test-data`
      Ignored local data drop zone for BAM/BAI, FASTQ, VCF, and expected output
      files copied in from elsewhere.

## Port Strategy

- [x] Treat upstream VNtyper as the behavioral reference.
- [x] Start with the smallest faithful path:
      BAM input -> MUC1 read extraction -> Kestrel VCF -> frameshift/depth
      classification -> TSV/JSON report.
- [x] Keep Kestrel as an external tool adapter first.
- [x] Keep samtools/bcftools/fastp/bwa as external tool adapters first, then
      replace the narrowest needed pieces with Rust wrappers when useful.
- [x] Keep optional modules separate:
      `adVNTR`, `SHARK`, cohort summaries, HTML reports, and mutation counter.
- [x] Prefer deterministic TSV/JSON parity tests before full HTML/report
      parity.

## Upstream Inventory

- [x] Read upstream CLI flow:
      `ports/vntyper/vntyper/vntyper/cli.py`.
- [x] Read upstream pipeline orchestration:
      `ports/vntyper/vntyper/vntyper/scripts/pipeline.py`.
- [x] Read Kestrel flow:
      `ports/vntyper/vntyper/vntyper/scripts/kestrel_genotyping.py`.
- [x] Read alignment/read extraction flow:
      `ports/vntyper/vntyper/vntyper/scripts/fastq_bam_processing.py`.
- [x] Read VCF and indel processing:
      `ports/vntyper/vntyper/vntyper/scripts/file_processing.py`,
      `variant_parsing.py`, and `motif_processing.py`.
- [x] Read scoring/confidence/filtering:
      `scoring.py`, `confidence_assignment.py`, `flagging.py`, and
      `kestrel_config.json`.
- [x] Read region/reference registry helpers:
      `region_utils.py`, `reference_registry.py`, and `chromosome_utils.py`.
- [x] Identify exact upstream outputs needed for parity:
      `kestrel_result.tsv`, `kestrel_pre_result.tsv`, filtered VCFs,
      pipeline summary JSON, and selected report fields.

## Test Data

- [x] Copy local VNtyper test data into `ports/vntyper/test-data`.
- [x] Inventory copied data:
      117 files, about 1.2 GiB, including hg19/hg38 subset BAM/BAI files,
      paired FASTQs, and remapped BWA BAM/BAI files across GRCh37/GRCh38,
      hg19/hg38, and Ensembl naming variants.
- [ ] Add or generate expected Kestrel VCF/TSV outputs for large integration
      data; copied data currently contains alignment/FASTQ inputs but no
      `.vcf`, `.tsv`, or result `.json` files. Tiny expected TSV/JSON fixtures
      exist for unit tests.
- [x] Mirror upstream `tests/test_data_config.json` filenames and MD5s in a
      BioScript-side manifest.
- [x] Add a data validator that checks required files.
- [x] Wire the data validator into integration tests so they skip with a clear
      message when large data is absent.
- [x] Keep large copied data out of git.
- [x] Add tiny synthetic VCF fixtures for unit tests that do not need BAM or
      Kestrel.

## BioScript Port Files

- [x] Add `ports/vntyper/bioscript/vntyper.bs.py` or equivalent top-level
      BioScript pipeline entry point.
- [x] Add BioScript modules for:
      region selection, command planning, Kestrel VCF parsing, frameshift
      classification, confidence assignment, and report row generation.
- [x] Add first BioScript-side post-processing module for Kestrel VCF parsing,
      frameshift classification, confidence assignment, and report JSON.
- [x] Keep BioScript code close to upstream naming where that helps parity.
- [x] Use `from bioscript import ...` imports for supported libraries and tool
      wrappers.
- [x] Avoid class-heavy ports until Monty class support is ready; use functions
      and plain dict/list records for the first pass.

## `bioscript-libs` Work

- [x] Add a `bioscript-libs::tools` or module-specific external tool wrapper
      layer with safe command construction.
- [x] Add `bioscript.samtools` wrapper surface for the VNtyper subset:
      `view`, `fastq`, `depth`, `index`, and possibly `faidx`.
- [x] Add `bioscript.bcftools` wrapper surface for optional VCF sort/compress
      fallback behavior.
- [x] Add `bioscript.kestrel` wrapper surface for invoking the vendored or
      configured Kestrel JAR.
- [x] Design `bioscript.kestrel` as a Python-shaped API rather than a direct
      Java clone. Initial surface:
      `kestrel.run(...)`, `kestrel.build_command(...)`, and
      `kestrel.read_vcf(...)`.
- [ ] Port the Kestrel Java internals only after the external-tool-backed
      wrapper passes VNtyper parity. Candidate internal packages:
      `counter`, `activeregion`, `align`, `variant`, and `writer.vcf`.
- [ ] Add `bioscript.fastp` wrapper surface only if FASTQ QC is in the first
      milestone.
- [ ] Add `bioscript.bwa` wrapper surface only if FASTQ input alignment is in
      the first milestone.
- [x] Add lightweight `bioscript.vcf` parsing helpers for Kestrel VCF rows.
- [x] Add TSV/CSV/table helpers if the port would otherwise need a pandas-like
      surface.

## Runtime / Security

- [x] Decide the external command policy for BioScript:
      allowlist commands, fixed argv builders, workspace-confined inputs, and
      controlled output paths.
- [x] Add runtime bindings for new modules imported via
      `from bioscript import samtools, kestrel, vcf, bcftools`.
- [x] Add runtime method bindings for `samtools` and `kestrel` command-builder
      calls.
- [x] Ensure `bioscript.kestrel` accepts structured arguments only; no arbitrary
      command strings.
- [x] Add tests that unsupported shell strings, remote paths, and write modes
      fail closed.
- [ ] Record tool execution in runtime trace/timing output.

## Python Compatibility Package

- [x] Add Python-side `bioscript.kestrel` command builders matching the Rust
      structured argv surface.
- [x] Add Python-side `bioscript.samtools` command builders matching the Rust
      structured argv surface.
- [x] Add Python tests for VNtyper tool command builders.

## Test Plan

- [ ] Port upstream unit tests first:
      confidence assignment, scoring, flagging, variant parsing, motif
      filtering, region utilities, chromosome utilities, and reference registry.
- [ ] Add parity tests that run the upstream Python function and BioScript port
      on the same tiny fixture and compare TSV/JSON values.
- [ ] Add integration tests against `ports/vntyper/test-data` once copied:
      one positive BAM, one negative BAM, and one FASTQ pair if available.
- [ ] Run upstream VNtyper tests from the submodule as a reference check when
      Python dependencies and external tools are installed.
- [x] Run BioScript tests without external tools by using fixed Kestrel VCF
      fixtures.
- [ ] Run full pipeline tests only when Kestrel/samtools/test data are present.

## Reporting / UI Parity

- [ ] Treat upstream `generate_report.py`, `report_template.html`, and
      `report_config.json` as the reporting reference.
- [x] Emit a structured BioScript report JSON before generating HTML.
- [ ] Include run metadata:
      report date, VNtyper version, input files, alignment pipeline, detected
      assembly/contig, and BAM header warnings.
- [ ] Include VNTR coverage QC:
      mean, median, stdev, min, max, region length, uncovered bases, percent
      uncovered, and pass/warning status.
- [ ] Include fastp QC when available:
      sequencing setup, duplication rate, Q20 rate, Q30 rate, passed-filter read
      rate, and threshold pass/warning status.
- [ ] Include screening summary logic from `report_config.json`:
      Kestrel result, optional adVNTR result, quality pass/fail, and validation
      recommendations.
- [ ] Include cross-match summary when adVNTR results are present.
- [x] Include Kestrel identified variants table:
      motif, variant, position, REF, ALT, motif sequence, variant depth,
      active-region depth, depth score, confidence, and flag.
- [ ] Include adVNTR identified variants table when available:
      VID, variant, supporting reads, mean coverage, p-value, RU, POS, REF,
      ALT, and flag.
- [ ] Preserve interactive HTML features after JSON parity:
      searchable/sortable tables, show/hide flagged rows, colored confidence
      values, flag icons/tooltips, detailed coverage toggle, and collapsible
      pipeline log.
- [ ] Add IGV visualization after core report parity:
      embedded IGV.js, variant selector table, and BAM/VCF track sessions.
- [ ] Make the first BioScript HTML report useful without IGV or adVNTR:
      final screening summary, Kestrel table, VNTR coverage QC, metadata, and
      pipeline log.

## Milestones

- [x] M1: Upstream source vendored and BioScript port skeleton committed.
- [x] M2: Kestrel VCF post-processing works in BioScript from fixture VCFs.
- [ ] M3: Confidence/depth/frame classification parity with upstream unit
      tests.
- [ ] M4: BAM path works using external samtools and Kestrel wrappers.
- [ ] M5: Native Rust Kestrel feasibility spike:
      reproduce Kestrel VCF output for one tiny fixture or document why the JVM
      adapter remains the practical first target.
- [ ] M6: Structured report JSON parity for the minimal BAM/Kestrel path.
- [ ] M7: HTML report parity for core summary, Kestrel table, coverage QC, and
      logs.
- [ ] M8: FASTQ path works using external fastp/bwa or documented prealigned
      inputs.
- [ ] M9: Optional adVNTR/SHARK/cohort/report modules triaged.
- [ ] M10: IGV visualization parity.
- [ ] M11: Replace selected external-tool behavior with Rust/noodles wrappers
      where the benefit is clear.

## Open Decisions

- [x] Whether the first public BioScript API should be command-like:
      `vntyper.run(config)` or step-oriented:
      `vntyper.extract_reads`, `vntyper.call_kestrel`, `vntyper.classify`.
- [x] Whether Kestrel is stored under `ports/vntyper/test-data/tools`, resolved
      from `PATH`, or configured via an environment/runtime option.
- [x] Whether pandas-like table operations should become `bioscript.table` or
      remain VNtyper-local helper functions.
- [x] Whether VNtyper references should be copied into BioScript-owned fixtures
      or read from the upstream submodule reference directory.
