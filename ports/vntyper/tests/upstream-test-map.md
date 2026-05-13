# Upstream VNtyper Test Map

Reference source: `ports/vntyper/vntyper/tests`.

This map decides where each upstream VNtyper test area belongs in the BioScript
port. The goal is not to run upstream pytest verbatim; it is to preserve the
same behavior with tests at the right layer: BioScript runtime, `bioscript-libs`
facade, Rust engine crate, or VNtyper-port logic.

## Integration And Orchestration

| Upstream file | BioScript mapping | Status |
| --- | --- | --- |
| `test_orchestration.py` | Port to BioScript/VNtyper large-data gates. BAM, FASTQ, and optional adVNTR runners should map to BioScript runner functions or runtime program execution. | Partial: BAM native gate exists; FASTQ native parity and adVNTR remain open. |
| `integration/test_pipeline_integration.py` | Port to opt-in large-data parity tests under `ports/vntyper/tests`. | Partial: external/native BAM gates exist; FASTQ and full upstream output checks remain open. |
| `docker/test_docker_pipeline.py` | Out of scope for BioScript core; replace with native binary/runtime smoke tests if BioScript gets a container image. | Deferred. |
| `parametrization.py` | Keep equivalent manifest-driven case selection in `ports/vntyper/tests/data_manifest.py`. | Partial. |
| `test_data_utils.py` | Keep only local manifest validation and skip messages. BioScript should not auto-download large data during normal tests. | Covered by `test_data_manifest.py`; checksum/download behavior is out of scope. |

## Unit Behavior

| Upstream file | BioScript mapping | Status |
| --- | --- | --- |
| `unit/test_alignment_processing.py` | `bioscript-libs` Samtools facade tests plus VNtyper command-plan tests. Exact FASTQ parity belongs in `samtools-rs`. | Partial. |
| `unit/test_bcftools_optional.py` | `bioscript-libs` BCFtools facade tests and Python wrapper tests. | Partial; native sort/index covered, optional filter expression execution deferred unless needed. |
| `unit/test_chromosome_utils.py` | Port to `ports/vntyper/tests/test_vntyper_regions.py` or config tests. | Partial. |
| `unit/test_confidence_assignment.py` | Port to VNtyper post-processing tests. | Partial. |
| `unit/test_flagging.py` | Port to VNtyper post-processing/report tests. | Partial. |
| `unit/test_grch_support.py` | Port to region/config tests and BAM/FASTQ parity cases for hg19/hg38. | Partial. |
| `unit/test_haplo_count_and_selection.py` | Port to VNtyper post-processing tests; engine-specific haplotype behavior belongs in `kestrel-rs`. | Partial. |
| `unit/test_install_references.py` | Mostly out of scope; BioScript uses vendored/reference paths rather than installing upstream reference bundles at runtime. | Deferred. |
| `unit/test_motif_filtering_issue_136.py` | Port directly to VNtyper post-processing tests. | Partial. |
| `unit/test_reference_registry.py` | Port to VNtyper config tests. | Partial. |
| `unit/test_region_utils.py` | Port to `test_vntyper_regions.py` and config tests. | Partial. |
| `unit/test_scoring.py` | Port directly to VNtyper post-processing tests and upstream scoring parity tests. | Partial. |
| `unit/test_utils.py` | Split by behavior: path/config behavior to VNtyper tests, command behavior to facade tests, unrelated CLI helpers out of scope. | Open. |
| `unit/test_variant_parsing.py` | Port directly to VNtyper VCF parsing/post-processing tests; Rust VCF parsing tests should be added if logic moves to `bioscript-libs`. | Partial. |

## Benchmark Tests

| Upstream file | BioScript mapping | Status |
| --- | --- | --- |
| `benchmark/*.py` | Out of scope for correctness. Add separate performance tracking only after parity is complete. | Deferred. |

## Required New BioScript Tests

- Runtime test executing the final `ports/vntyper/bioscript/vntyper.bs` program
  on tiny checked-in fixtures.
- Rust `bioscript-libs` test for native Samtools/Kestrel/BCFtools orchestration
  on tiny fixtures.
- Opt-in BAM large-data parity for positive and negative fixtures.
- Opt-in FASTQ large-data parity for positive and negative fixtures.
- JSON/TSV normalized comparisons with explicit ignored fields for paths,
  timestamps, and tool-version metadata.
- HTML report structure comparison from generated report JSON.

## Known Engine-Owned Gaps

- `samtools-rs`: exact parity for `view -P | sort -n | fastq -1/-2/-0/-s`.
- `kestrel-rs`: Java Kestrel parity for VNtyper positive/negative FASTQ
  expected VCF records.
- `bcftools-rs`: native `view -i/-e` expression execution only if a future
  BioScript VNtyper path needs it.
