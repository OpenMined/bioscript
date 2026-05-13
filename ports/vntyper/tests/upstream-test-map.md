# Upstream VNtyper Test Map

Reference source: `ports/vntyper/vntyper/tests`.

This map decides where each upstream VNtyper test area belongs in the BioScript
port. The goal is not to run upstream pytest verbatim; it is to preserve the
same behavior with tests at the right layer: BioScript runtime, `bioscript-libs`
facade, Rust engine crate, or VNtyper-port logic.

## Integration And Orchestration

| Upstream file | BioScript mapping | Status |
| --- | --- | --- |
| `test_orchestration.py` | Port to BioScript/VNtyper large-data gates. BAM, FASTQ, and optional adVNTR runners should map to BioScript runner functions or runtime program execution. | Partial: BAM native gate exists and passes classification parity; `vntyper-fastq.bs` runs native Kestrel/BCFtools through the runtime on tiny fixtures; FASTQ large-data classification parity is blocked by `kestrel-rs`; adVNTR remains deferred. |
| `integration/test_pipeline_integration.py` | Port to opt-in large-data parity tests under `ports/vntyper/tests`. | Partial: external/native BAM gates exist; FASTQ runtime execution exists for tiny fixtures, but large-data parity is blocked by Kestrel; full TSV/report output parity remains open. |
| `docker/test_docker_pipeline.py` | Out of scope for BioScript core; replace with native binary/runtime smoke tests if BioScript gets a container image. | Deferred. |
| `parametrization.py` | Keep equivalent manifest-driven case selection in `ports/vntyper/tests/data_manifest.py`. | Covered for current positive/negative BAM and FASTQ representative cases by `data_manifest.py` and skip-message tests; upstream download/checksum behavior is intentionally out of scope for normal BioScript tests. |
| `test_data_utils.py` | Keep only local manifest validation and skip messages. BioScript should not auto-download large data during normal tests. | Covered by `test_data_manifest.py`; checksum/download behavior is out of scope. |

## Unit Behavior

| Upstream file | BioScript mapping | Status |
| --- | --- | --- |
| `unit/test_alignment_processing.py` | `bioscript-libs` Samtools facade tests plus VNtyper command-plan tests. Exact FASTQ parity belongs in `samtools-rs`. | Covered for VNtyper-required behavior: native facade tests cover tiny BAM/index handling, and the opt-in samtools oracle gate verifies `view -P | sort -n | fastq -1/-2/-0/-s` against real samtools for representative fixtures. |
| `unit/test_bcftools_optional.py` | `bioscript-libs` BCFtools facade tests and Python wrapper tests. | Covered for VNtyper-required behavior: native sort/index and VCF materialization are tested. Optional native `view -i/-e` expression execution is deferred because the current VNtyper port filters records in port logic rather than through bcftools expressions. |
| `unit/test_chromosome_utils.py` | Port to `ports/vntyper/tests/test_vntyper_regions.py` or config tests. | Covered for VNtyper-required naming conventions by `test_vntyper_regions.py`; upstream pytest subset gate includes the upstream file when dependencies are installed. |
| `unit/test_confidence_assignment.py` | Port to VNtyper post-processing tests. | Covered for current thresholds and boundary behavior by `test_ported_upstream_units.py`, `test_vntyper_port.py`, and `test_upstream_scoring_parity.py`; upstream pytest subset gate includes the upstream file when dependencies are installed. |
| `unit/test_flagging.py` | Port to VNtyper post-processing/report tests. | Covered for rule evaluation, duplicate flags, and report visibility by `test_ported_upstream_units.py` and `test_vntyper_report.py`; keep expanding if new upstream flag rules are added. |
| `unit/test_grch_support.py` | Port to region/config tests and BAM/FASTQ parity cases for hg19/hg38. | Partial: hg19/hg38 coordinate/config behavior is covered by `test_vntyper_regions.py` and `test_vntyper_config.py`; representative large-data gates currently exercise hg19 fixtures only. |
| `unit/test_haplo_count_and_selection.py` | Port to VNtyper post-processing tests; engine-specific haplotype behavior belongs in `kestrel-rs`. | Partial: VNtyper best-call selection is covered by port tests; Kestrel haplotype parity is represented by the opt-in `kestrel-rs` FASTQ parity gate and currently fails against Java expected VCF counts. |
| `unit/test_install_references.py` | Mostly out of scope; BioScript uses vendored/reference paths rather than installing upstream reference bundles at runtime. | Deferred. |
| `unit/test_motif_filtering_issue_136.py` | Port directly to VNtyper post-processing tests. | Covered for current right/left motif filtering and issue-style conserved motif exclusions by `test_ported_upstream_units.py`. |
| `unit/test_reference_registry.py` | Port to VNtyper config tests. | Covered for current explicit reference paths and report schema config by `test_vntyper_config.py`; upstream install/download behavior is deferred with `unit/test_install_references.py`. |
| `unit/test_region_utils.py` | Port to `test_vntyper_regions.py` and config tests. | Covered for assembly aliases, coordinate strings, contig naming conventions, and invalid coordinates by `test_vntyper_regions.py`; upstream pytest subset gate includes the upstream file when dependencies are installed. |
| `unit/test_scoring.py` | Port directly to VNtyper post-processing tests and upstream scoring parity tests. | Covered for frame scoring, frameshift extraction, confidence assignment, depth score, and upstream subset parity by `test_ported_upstream_units.py`, `test_vntyper_port.py`, and `test_upstream_scoring_parity.py`; upstream pytest subset gate includes the upstream file when dependencies are installed. |
| `unit/test_utils.py` | Split by behavior: path/config behavior to VNtyper tests, command behavior to facade tests, unrelated CLI helpers out of scope. | Partial: sample-name/path validation and manifest skip behavior are covered by `test_vntyper_commands.py` and `test_data_manifest.py`; remaining unrelated CLI helper behavior should stay out of BioScript core unless the final runtime CLI needs it. |
| `unit/test_variant_parsing.py` | Port directly to VNtyper VCF parsing/post-processing tests; Rust VCF parsing tests should be added if logic moves to `bioscript-libs`. | Covered for VNtyper-required VCF parsing, ALT filtering, named sample columns, expected TSV rows, and expected report summary by `test_vntyper_port.py`, `test_ported_upstream_units.py`, and `test_upstream_scoring_parity.py`. Core Rust call-table conversion is covered by `rust/bioscript-libs/tests/vntyper_vcf.rs`; upstream pytest subset gate includes the upstream file when dependencies are installed. |

## Benchmark Tests

| Upstream file | BioScript mapping | Status |
| --- | --- | --- |
| `benchmark/*.py` | Out of scope for correctness. Add separate performance tracking only after parity is complete. | Deferred. |

## Required New BioScript Tests

- Runtime tests executing BioScript VNtyper programs: covered by
  `rust/bioscript-runtime/tests/vntyper_program.rs`. `vntyper.bs` is still a
  BAM command-plan execution test. `vntyper-fastq.bs` now runs native
  Kestrel/BCFtools/VNtyper Kestrel call-table parsing on tiny generated
  FASTQ/reference fixtures, writes `kestrel_result.tsv` plus a TSV summary, and
  materializes report JSON through the VCF facade. Full TSV/JSON/HTML parity
  remains open.
- Rust `bioscript-libs` test for native Samtools/Kestrel/BCFtools orchestration
  on tiny fixtures: covered by `rust/bioscript-libs/tests/vntyper_facades.rs`.
- Opt-in BAM large-data parity for positive and negative fixtures: covered by
  `ports/vntyper/tests/test_native_bam_pipeline_gate.py`; classification parity
  passes.
- Opt-in FASTQ large-data parity for positive and negative fixtures: covered by
  `ports/vntyper/tests/test_native_fastq_pipeline_gate.py`; the gate runs but
  currently fails because native Kestrel output differs from Java expected data.
- JSON/TSV normalized comparisons with explicit ignored fields for paths,
  timestamps, and tool-version metadata: open. Current BAM generated TSV row
  counts differ from expected fixtures even when report summary classification
  matches.
- HTML report structure comparison from generated report JSON: covered by
  `ports/vntyper/tests/test_vntyper_report.py`.

## Known Engine-Owned Gaps

- `kestrel-rs`: Java Kestrel parity for VNtyper positive/negative FASTQ
  expected VCF records. Reduced into
  `vendor/rust/kestrel-rs/crates/kestrel/tests/vntyper_fastq_parity.rs`;
  opt-in failures currently show fewer Rust records than Java expected records.
- `bcftools-rs`: native `view -i/-e` expression execution only if a future
  BioScript VNtyper path needs it.
