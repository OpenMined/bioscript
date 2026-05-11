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
- [x] Add or generate expected Kestrel VCF/TSV outputs for large integration
      data; copied data currently contains alignment/FASTQ inputs but no
      `.vcf`, `.tsv`, or result `.json` files. Tiny expected TSV/JSON fixtures
      exist for unit tests. A dry-run generator now exists at
      `ports/vntyper/tests/generate_expected_outputs.py`; it records sample
      labels, planned commands, and the ignored expected-output layout. Without
      `--dry-run`, it uses the external pipeline runner to materialize VCF, TSV,
      and JSON outputs once local samtools/bcftools/Kestrel prerequisites and
      validated sample labels are available. `--fastq-only` can bootstrap
      Kestrel VCF/TSV/report outputs from existing copied FASTQ pairs without
      samtools/bcftools. Generated local ignored FASTQ-backed outputs now exist
      under `ports/vntyper/test-data/expected/{positive,negative}`.
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
- [x] Add a BioScript-owned Kestrel build helper for environments without Ant:
      `ports/vntyper/tests/build_kestrel_jar.py` compiles the vendored Java
      sources with Java 8 compatibility and packages an ignored local
      `ports/vntyper/test-data/tools/kestrel/kestrel.jar` for integration
      tests.
- [ ] Port the Kestrel Java internals only after the external-tool-backed
      wrapper passes VNtyper parity. Candidate internal packages:
      `counter`, `activeregion`, `align`, `variant`, and `writer.vcf`.
      The first native surface now exists in
      `rust/bioscript-libs/src/kestrel/native.rs`: a Rust Kestrel VCF writer
      model that mirrors the Java `writer.vcf` headers, FORMAT fields,
      multi-sample genotype/depth fields, validation, and record ordering.
      Native SNP/insertion/deletion VCF normalization now mirrors the Java
      `variant` package rules for VCF POS/REF/ALT anchoring against a reference
      region. The native Rust code has also been split into focused
      `native::{kmer, variant, vcf}` modules, and `native::KmerCountMap`
      provides the first `counter` package equivalent for counting canonical
      sequence k-mers from in-memory strings, FASTQ, and FASTQ.gz inputs. The
      native `ActiveRegion` and `RegionStats` types now mirror the Java
      `activeregion` data model for anchor k-mers and percentile count
      summaries. A first native `detect_active_regions` candidate scanner now
      computes reference k-mer counts and Java-shaped difference thresholds,
      then emits anchored and right-open depth-drop regions for downstream
      haplotype work. It now also exposes Java-shaped `anchor_both_ends`
      behavior, defaults reverse-kmer counting and both-end anchoring to the
      Java detector defaults, and emits left-open candidates for near-left-end
      active regions when unanchored ends are explicitly allowed. Java's
      exponential recovery-threshold shape is now implemented with `decay_min`
      and `decay_alpha` controls and exposed through the native/Python wrapper
      path. The right-scan peak detection heuristic is also partially ported:
      `peak_scan_length` controls stable-recovery scanning and short recovery
      spikes inside a low-count valley no longer prematurely terminate the
      active region. The native detector also exposes a Java-shaped
      `scan_limit_factor` control plus an explicit `max_gap_size` input for the
      Java `maxGapSize + scanLimitFactor * k` shape; BioScript now ports the
      Java default `AlignmentWeight.getMaxExclusiveGapSize(k)` calculation and
      uses it as the native wrapper default when callers do not provide an
      explicit gap component. The Java alignment-weight vector parser shape is
      also ported for default/partial vectors, surrounding bounds, sign
      normalization, and Java integer literal formats. Both left and right
      scans discard candidates that exceed that limit. Java's
      default `recoverRightAnchor` behavior is now partially ported as
      `recover_right_anchor`: when the normal recovery threshold is never
      reached inside the scan limit, the native detector searches for a later
      abrupt count increase and uses that k-mer as a recovered right anchor.
      The first left-scan peak suppression
      rule is also ported: short isolated count increases can be skipped rather
      than being emitted as left-end active regions, and left-open candidate
      scans now respect the same scan-limit length used by right scans. The
      left-scan recovery check now also follows Java's discard shape when
      counts recover before the scan reaches the left end, which prevents those
      internal recoveries from being emitted as left-end active regions.
      Java's `callAmbiguousRegions` switch is now exposed as
      `call_ambiguous_regions` through the native detector and Python wrapper,
      with default-on behavior and optional rejection of active regions whose
      reference span contains ambiguous bases. The
      native `align_haplotype` and `call_alignment_variants` helpers provide a
      first deterministic reference-vs-haplotype edit surface that emits
      SNP/insertion/deletion calls using the same native VCF normalization path.
      The upstream compiled Kestrel JUnit reference-reader fixture set has also
      been ported into Rust tests: native reference parsing now covers FASTA,
      FASTQ, mixed case, legal IUPAC/gap characters, and Kestrel's deterministic
      ambiguous-base-to-ACGT k-mer normalization for k sizes 1, 2, 21, 32, and
      64.
      The separate upstream `paudano/kescases` publication pipeline is now
      vendored as `ports/vntyper/kescases` for the next parity layer; it
      contains Kestrel CLI/Snakemake workflows, bundled Kestrel jars, reference
      FASTA data, and comparison pipelines rather than ordinary unit-test
      sources.
      `call_explicit_haplotypes_to_vcf` now ties explicit haplotype evidence to
      the native aligner, variant caller, and VCF writer for an end-to-end
      non-assembling caller path. The first graph-backed Rust haplotype
      assembler now walks counted k-mer paths between active-region anchors and
      feeds assembled haplotypes into the native VCF caller. The native
      `call_sequences_to_vcf` path now ties read sequence counting, active-region
      detection, graph haplotype assembly, alignment, variant calling, and VCF
      writing together for small synthetic fixtures, and
      `bioscript.kestrel.call_sequences_native` exposes that path through the
      Python wrapper/PyO3 layer. `call_fastq_paths_to_vcf` and
      `bioscript.kestrel.call_fastq_native` extend the same native caller to
      FASTQ inputs produced by the BioScript samtools extraction path. A
      multi-reference native VCF path now counts FASTQ reads once, emits all
      reference contig headers, and scans each reference region for variants,
      with Python/PyO3 wrapper access through
      `bioscript.kestrel.call_fastq_references_native`. Python-side
      `bioscript.kestrel.load_reference_regions` reads multi-record FASTA files
      into `(name, sequence, md5)` triples for that native path, matching the
      shape of VNtyper motif dictionaries. The Java parity gate now includes a
      multi-reference FASTQ fixture that emits all contig headers and calls the
      matching reference record, which is the next required shape for full
      VNtyper motif-reference parity. The
      haplotype assembler now tracks repeated k-mers and trims saved states by
      path depth using exposed `max_repeat_count` and `max_saved_states`
      controls. A first opt-in Java parity gate now exists at
      `rust/bioscript-libs/tests/kestrel_java_parity.rs`; when
      `BIOSCRIPT_RUN_KESTREL_JAVA_PARITY=1` and a Kestrel jar are available, it
      compares native FASTQ-to-VCF output with Java Kestrel on tiny
      perfect-reference no-variant, MUC1 SNP, nonrepetitive SNP, adjacent
      nonrepetitive SNPs, k=20 nonrepetitive SNP/deletion/insertion fixtures,
      mixed reference/alternate SNP and deletion depth, a mixed insertion
      no-call, sparse split-read, and multi-reference fixtures.
      The native assembler now tracks observed adjacent k-mer transitions from
      each read/FASTQ record and refuses to bridge k-mers that were never
      adjacent in an input read, which fixes the Java-confirmed sparse
      reference-consistent case (`AAAACCC`, `CCCTGGG`, `GGGTTTT`) against
      `AAAACCCCGGGGTTTT`. It also assigns VCF DP from the total assembled
      active-region haplotype depth, matching Java's mixed reference/alternate
      depth shape. The remaining work is the full Java active-region detector
      heuristics and broader parity against Java Kestrel outputs on larger
      synthetic and VNtyper fixtures.
- [x] Add `bioscript.fastp` wrapper surface only if FASTQ QC is in the first
      milestone.
- [x] Add `bioscript.bwa` wrapper surface only if FASTQ input alignment is in
      the first milestone.
- [x] Add lightweight `bioscript.vcf` parsing helpers for Kestrel VCF rows.
- [x] Add TSV/CSV/table helpers if the port would otherwise need a pandas-like
      surface.
- [x] Add a first native noodles replacement slice:
      `pysam.AlignmentFile.fetch` now supports indexed BAM inputs through
      `bioscript-formats::alignment::query_bam_records`, while CRAM continues
      through the existing noodles CRAM path.

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
- [x] Record tool execution in runtime trace/timing output.

## Python Compatibility Package

- [x] Add Python-side `bioscript.kestrel` command builders matching the Rust
      structured argv surface.
- [x] Add Python-side `bioscript.samtools` command builders matching the Rust
      structured argv surface.
- [x] Add Python tests for VNtyper tool command builders.

## Test Plan

- [x] Port upstream unit tests first:
      confidence assignment, scoring, flagging, variant parsing, motif
      filtering, region utilities, chromosome utilities, and reference registry.
- [x] Add parity tests that run the upstream Python function and BioScript port
      on the same tiny fixture and compare TSV/JSON values.
- [x] Add integration tests against `ports/vntyper/test-data` once copied:
      one positive BAM, one negative BAM, and one FASTQ pair if available.
      Current coverage plans commands for two representative BAMs and one FASTQ
      pair, and a fake-runner test covers the BAM path running slice, index,
      FASTQ extraction, depth, Kestrel, bcftools, and TSV/JSON materialization.
      A second fake-runner path now covers native BioScript samtools slice,
      FASTQ extraction, and depth followed by Kestrel without requiring
      bcftools. A gated real-data native BAM pipeline test now exists and skips
      until explicitly enabled with `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1` and
      `bioscript._native`, Java/Kestrel, BAM/BAI inputs, and expected outputs
      are all available. A separate `samtools` oracle gate now exists at
      `ports/vntyper/tests/test_samtools_fastq_oracle.py`; it is opt-in with
      `BIOSCRIPT_RUN_SAMTOOLS_ORACLE=1` and compares native FASTQ extraction
      counts against `samtools view -P`, name-sort, and `samtools fastq`.
      The local environment is Arch Linux, and `sudo pacman -Sy --needed
      --noconfirm samtools bcftools` cannot run non-interactively here because
      sudo requires a terminal password. To unblock comparison gates, local
      ignored builds of `htslib`, `samtools`, and `bcftools` 1.23.1 were built
      under `ports/vntyper/test-data/tools/local`; the manifest discovers those
      binaries when system installs are absent.
      FASTQ-backed Kestrel expected outputs are gated by
      `test_fastq_expected_outputs.py`; native BAM-backed positive and negative
      representative samples are gated by `test_native_bam_pipeline_gate.py`.
      The native BAM gate verifies sample classification, report shape,
      screening summary, nonempty Kestrel rows, variant-table linkage, and VNTR
      coverage fields against the generated expected report set. The external
      BAM gate in `test_full_pipeline_gate.py` is opt-in with
      `BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1` and runs the samtools/bcftools plus
      Kestrel path for the representative positive and negative BAMs.
- [x] Run upstream VNtyper tests from the submodule as a reference check when
      Python dependencies and external tools are installed.
- [x] Run BioScript tests without external tools by using fixed Kestrel VCF
      fixtures.
- [x] Run full pipeline tests only when Kestrel/samtools/test data are present.

## Reporting / UI Parity

- [x] Treat upstream `generate_report.py`, `report_template.html`, and
      `report_config.json` as the reporting reference.
- [x] Emit a structured BioScript report JSON before generating HTML.
- [x] Include run metadata:
      report date, VNtyper version, input files, alignment pipeline, detected
      assembly/contig, and BAM header warnings.
- [x] Include VNTR coverage QC:
      mean, median, stdev, min, max, region length, uncovered bases, percent
      uncovered, and pass/warning status.
- [x] Include fastp QC when available:
      sequencing setup, duplication rate, Q20 rate, Q30 rate, passed-filter read
      rate, and threshold pass/warning status.
- [x] Include screening summary logic from `report_config.json`:
      Kestrel result, optional adVNTR result, quality pass/fail, and validation
      recommendations.
- [x] Include cross-match summary when adVNTR results are present.
- [x] Include Kestrel identified variants table:
      motif, variant, position, REF, ALT, motif sequence, variant depth,
      active-region depth, depth score, confidence, and flag.
- [x] Include adVNTR identified variants table when available:
      VID, variant, supporting reads, mean coverage, p-value, RU, POS, REF,
      ALT, and flag.
- [x] Preserve interactive HTML features after JSON parity:
      searchable/sortable tables, show/hide flagged rows, colored confidence
      values, flag icons/tooltips, detailed coverage toggle, and collapsible
      pipeline log.
- [x] Add IGV visualization after core report parity:
      embedded IGV.js, variant selector table, and BAM/VCF track sessions.
- [x] Make the first BioScript HTML report useful without IGV or adVNTR:
      final screening summary, Kestrel table, VNTR coverage QC, metadata, and
      pipeline log.

## Milestones

- [x] M1: Upstream source vendored and BioScript port skeleton committed.
- [x] M2: Kestrel VCF post-processing works in BioScript from fixture VCFs.
- [x] M3: Confidence/depth/frame classification parity with upstream unit
      tests.
- [x] M4: BAM path works using external samtools and Kestrel wrappers.
      The execution layer now exists in
      `ports/vntyper/bioscript/vntyper_external_pipeline.py` and is covered
      with an injected fake runner. Local ignored `htslib`, `samtools`, and
      `bcftools` 1.23.1 builds provide comparison tools when system packages are
      unavailable. The opt-in external BAM gate runs the real-tool path against
      representative positive and negative BAMs, requires nonempty Kestrel rows,
      and compares classification/report shape with generated expected reports.
      Native BioScript BAM FASTQ extraction now writes complete primary R1/R2
      pairs only and matches the copied representative FASTQ fixture counts
      for `example_6449_hg19_subset.bam` (`82523/82523`) and
      `example_66bf_hg19_subset.bam` (`19877/19877`). The native BAM/Kestrel
      gate now passes locally when explicitly enabled with
      `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1` and a temporarily copied
      `bioscript._native` extension.
- [x] M5: Native Rust Kestrel feasibility spike:
      reproduce Kestrel VCF output for one tiny fixture or document why the JVM
      adapter remains the practical first target.
- [x] M6: Structured report JSON parity for the minimal BAM/Kestrel path.
      Fake-runner coverage now captures `samtools depth -a` output and feeds
      mean/median/stdev/min/max/uncovered-base fields into the structured JSON;
      FASTQ-backed Kestrel reports are now generated locally, and the runner
      can use native BioScript samtools wrappers before Kestrel. The opt-in
      native BAM gate now validates copied positive and negative BAM samples
      against generated expected reports, including report schema, Kestrel
      classification, screening summary, variant-table linkage, and populated
      VNTR coverage metrics. The external `samtools`/`bcftools` gate also runs
      locally through the ignored user-space tool build.
- [x] M7: HTML report parity for core summary, Kestrel table, coverage QC, and
      logs.
- [x] M8: FASTQ path works using external fastp/bwa or documented prealigned
      inputs.
- [x] M9: Optional adVNTR/SHARK/cohort/report modules triaged.
- [x] M10: IGV visualization parity.
- [x] M11: Replace selected external-tool behavior with Rust/noodles wrappers
      where the benefit is clear.
      Selected replacements now cover indexed BAM region fetch, indexed BAM
      depth summary, BAM region slicing, and BAM-region-to-paired-FASTQ
      extraction through noodles, with `bioscript-python` native samtools
      wrappers for CPython tests. Remaining future candidates are VCF
      sorting/indexing.

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
