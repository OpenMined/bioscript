# BioScript Library Support TODO

Goal: make BioScript support recognizable bioinformatics library/tool surfaces
through thin `bioscript-libs` facades backed by vendored Rust engine crates.
Build the reusable primitives first, wire Samtools next, and then make the
VNtyper BioScript port a small amount of pipeline code plus data/config that
uses those built-in primitives.

## Direction

- [x] Use explicit BioScript imports:
      `from bioscript import samtools, bcftools, kestrel, pysam, pyfaidx`.
- [ ] Treat BioScript library support as the product:
      common pipeline code should read like standard bioinformatics workflows,
      not like private BioScript internals.
- [ ] Build in layers:
      engine crates -> BioScript facades -> facade tests -> VNtyper port.
- [x] Keep BioScript-owned code as compatibility/adaptation code, not full
      algorithm ports.
- [x] Put heavy native implementations in reusable Rust engine crates under
      `vendor/rust`.
- [x] Keep upstream Python API references under `vendor/python`.
- [ ] Refactor existing BioScript methods to call these higher-level facades
      instead of private lower-level helpers where the public bioinformatics
      name is clearer.

## Vendor Layout

- [x] Move Python reference submodules:
      `vendor/python/pysam`
      `vendor/python/pyfaidx`
- [x] Add Kestrel Rust engine:
      `vendor/rust/kestrel-rs`
- [x] Add HTS Rust engine:
      `vendor/rust/htslib-rs`
- [x] Add BCFtools Rust engine:
      `vendor/rust/bcftools-rs`
- [ ] Add Samtools Rust engine once ready:
      `vendor/rust/samtools-rs`
- [ ] Keep vendored engine crate tests inside their own repos/workspaces.
- [ ] Keep BioScript tests focused on adapter behavior and pipeline integration.

## Rust Crate Wiring

- [x] Wire `rust/bioscript-libs` to local `kestrel-rs` path dependencies:
      `kestrel` and `kanalyze`.
- [ ] Wire `rust/bioscript-libs` to local `htslib-rs`.
      Top-level `vendor/rust/htslib-rs` is present, but direct Cargo wiring is
      deferred until the duplicated nested `htslib-rs` dependency inside
      `bcftools-rs` is unified.
- [x] Wire `rust/bioscript-libs` to local `bcftools-rs`.
- [ ] Wire `rust/bioscript-libs` to local `samtools-rs` when available.
- [ ] Add `[patch]` entries only where nested crate dependencies would
      otherwise pull remote git/crates.io versions instead of local submodules.
- [x] Document the dependency graph:
      BioScript -> `bioscript-libs` facade -> vendored Rust engine crate.

## Crate Publishing

- [ ] Keep local path dependencies while `kestrel-rs`, `htslib-rs`,
      `bcftools-rs`, and `samtools-rs` APIs are still changing quickly.
- [ ] Publish those engine crates once their public APIs and test suites are
      stable enough for external consumers.
- [ ] After publishing, replace stable path dependencies with versioned crates
      where that simplifies the Cargo graph.
- [ ] Keep submodules available for upstream test fixtures, source comparison,
      and local patching even after published crates are used by default.

## Milestones

- [x] M1: Kestrel Rust engine is vendored and callable through BioScript.
- [ ] M2: HTS and BCFtools Rust engines are vendored and wired by path.
      Both engines are vendored. BCFtools is wired into `bioscript-libs`;
      top-level HTS direct wiring is still pending dependency unification.
- [ ] M3: Samtools Rust engine is vendored and wired by path.
- [ ] M4: BioScript facades expose a minimal, recognizable built-in toolkit:
      `samtools`, `bcftools`, `kestrel`, `pysam`, `pyfaidx`, and VCF/table
      helpers.
- [ ] M5: Existing BioScript lower-level helper paths are refactored to use the
      public facades where possible.
- [ ] M6: VNtyper is reimplemented as a small BioScript pipeline that mostly
      coordinates built-in primitives and carries only VNtyper-specific
      constants, motif data, filtering rules, and report logic.

## Kestrel Facade

- [x] Remove old in-tree custom Rust Kestrel algorithm modules from
      `rust/bioscript-libs/src/kestrel/native/`.
- [x] Replace them with `rust/bioscript-libs/src/kestrel/native.rs`, a thin
      adapter around `vendor/rust/kestrel-rs`.
- [x] Preserve the Python-facing API names used by VNtyper:
      `call_sequences_native`, `call_fastq_native`,
      `call_fastq_references_native`.
- [x] Add adapter support for `.fastq.gz` inputs by normalizing them before
      calling `kestrel-rs`.
- [x] Remove the stale BioScript Java-parity test that targeted the deleted
      in-tree Kestrel internals.
- [x] Add small deterministic adapter tests proving `kestrel-rs` emits an
      expected SNP VCF through the BioScript facade.
- [x] Decide whether BioScript should expose a more direct `kestrel.run(...)`
      path that writes output files, or keep the current string-returning VCF
      helpers for Python/VNtyper integration.
      Decision: keep string-returning low-level helpers and expose
      `kestrel.run_native(...)` as the file-writing convenience path.
- [ ] Move any remaining Kestrel algorithm parity expectations into
      `vendor/rust/kestrel-rs`.

## Samtools Facade

- [x] Existing BioScript command-builder surface:
      `samtools.view_region`, `samtools.fastq`, `samtools.depth`.
- [x] Existing native prototype supports BAM slicing, FASTQ extraction, and
      depth summary through BioScript-owned primitives.
- [ ] Replace native prototype internals with calls into `samtools-rs` once the
      crate is available.
- [ ] Prioritize Samtools after vendoring HTS/BCFtools because VNtyper's BAM
      path should become:
      `samtools.view` -> `samtools.index/sort` if needed ->
      `samtools.fastq` -> `samtools.depth`.
- [ ] Keep the public BioScript API shaped like familiar samtools operations:
      `view`, `fastq`, `sort`, `index`, `depth`, `faidx`.
- [ ] Add adapter tests for:
      region parsing, indexed BAM input, `.bam/.bai` discovery, paired FASTQ
      output counts, depth summary fields, and error mapping.
- [ ] Keep oracle tests against real samtools opt-in only.

## BCFtools Facade

- [x] Existing BioScript command-builder surface:
      `bcftools.sort`, `bcftools.view_filter`.
- [x] Add `vendor/rust/bcftools-rs`.
- [x] Inspect the `bcftools-rs` public API and choose the thinnest adapter
      surface for VNtyper.
- [x] Replace command-only behavior with native calls where the Rust crate
      supports them.
      Initial native methods: `view_header_native`, `view_native`, and
      `index_native`, backed by `bcftools_rs::commands::{view,index}`.
      Native sort now calls `bcftools_rs::commands::sort` for the VNtyper
      `sort -o output.vcf.gz -W -O z` path.
- [ ] Initial target operations:
      `view`, `sort`, `norm`, compression/index helpers if needed.
- [ ] Add adapter tests for VCF input/output, compressed output, filter
      expressions used by VNtyper, and useful error messages.
      Initial coverage verifies `bcftools-rs` header extraction, VCF output,
      BGZF-compressed output, native sort, CSI/TBI indexing, Python wrapper
      delegation, and the real PyO3 native extension when installed. Filter
      expression coverage remains pending until `bcftools-rs view` supports
      `-i/-e`.

## HTS / Pysam / Pyfaidx Facades

- [x] Keep `pysam` and `pyfaidx` as recognizable compatibility namespaces.
- [x] `pyfaidx.Fasta` has a small Rust/Python-compatible FASTA slice surface.
- [x] `pysam.AlignmentFile.fetch` has initial BAM/CRAM read support.
- [ ] Refactor lower-level alignment code to flow through `pysam` or
      `samtools` facades where that makes scripts more recognizable.
- [ ] Use `htslib-rs` as the shared backend for BAM/CRAM/VCF/FASTA primitives
      once vendored.
- [ ] Add parity tests from focused upstream `pysam` and `pyfaidx` cases, not
      the full upstream test suites.

## Python Package

- [x] Keep top-level `python/bioscript` matching BioScript import names.
- [x] Keep optional delegation to real Python libraries where useful.
- [x] Expose native functions through `rust/bioscript-python`.
- [ ] Add Python tests that call the real native extension for each engine
      facade with tiny fixtures.
- [x] Keep mocked-extension tests for argument normalization and missing-native
      behavior.
- [x] Make Python-only fallback behavior explicit per module:
      real Python library, pure Python fallback, or native-required.

## Runtime / Monty Integration

- [x] Support `from bioscript import x` import rewriting for current modules.
- [x] Bind initial module objects and method calls in `bioscript-runtime`.
- [ ] Add runtime method bindings for native samtools/bcftools operations once
      facades are stable.
      BCFtools native bindings now cover `view_header_native`, `view_native`,
      `sort_native`, and `index_native`; Samtools native bindings are still pending the
      `samtools-rs` backend.
- [ ] Keep runtime responsible for language/object adaptation only.
- [ ] Keep file/path/security policy centralized and reused across facades.

## VNtyper Proof Port

- [x] Keep upstream VNtyper source vendored at `ports/vntyper/vntyper`.
- [x] Keep local large test data ignored under `ports/vntyper/test-data`.
- [x] Keep BioScript VNtyper port under `ports/vntyper/bioscript`.
- [x] Keep BioScript-owned VNtyper tests under `ports/vntyper/tests`.
- [x] Current tests cover command planning, Kestrel VCF parsing, scoring,
      report JSON/HTML shape, and fake-runner pipeline behavior.
- [x] Current adapter smoke tests prove BioScript can call `kestrel-rs`.
- [ ] Reframe the final VNtyper port as its own BioScript code, not as a copy
      of every upstream dependency. The VNtyper-specific layer should contain:
      MUC1 regions, motif/reference data, Kestrel parameter choices,
      frameshift/depth classification, report rows, and CLI/pipeline glue.
- [ ] Keep generic work out of the VNtyper port. Generic work belongs in
      BioScript facades:
      BAM/CRAM slicing, FASTQ extraction, depth, VCF parsing/filtering,
      Kestrel calling, FASTA lookup, TSV/JSON helpers.
- [ ] Refactor VNtyper pipeline code to prefer:
      `samtools.*`, `bcftools.*`, `kestrel.*`, `pysam.*`, and `pyfaidx.*`
      over private helper names.
      Native Kestrel execution now goes through `kestrel.run_native(...)`
      instead of VNtyper manually loading references and writing VCF text.
      The FASTQ-only path can now optionally run native Kestrel followed by
      native BCFtools sort/index without Java or external bcftools.
- [x] Define the minimal VNtyper BioScript interface, for example:
      `run_vntyper(bam=..., reference_build="hg19", output_dir=...)` and
      `run_vntyper_fastq(r1=..., r2=..., reference_build="hg19", output_dir=...)`.
- [ ] Keep VNtyper data/config small and explicit:
      MUC1 coordinates, motif FASTA path, confidence thresholds, report schema,
      and optional validation toggles.
- [ ] Once `samtools-rs` and `bcftools-rs` are wired, rerun the BAM path using
      only BioScript native facades.
- [ ] Compare native-facade VNtyper output against expected positive/negative
      fixtures for:
      FASTQ path, BAM path, report JSON, and HTML report.
- [ ] Keep large real-data parity tests opt-in with clear skip messages.

## Test Policy

- [x] Engine crates own engine correctness:
      e.g. `vendor/rust/kestrel-rs` owns Kestrel Java/algorithm parity.
- [x] BioScript owns facade correctness:
      argument normalization, path handling, output shape, error mapping, and
      integration with BioScript/Python/VNtyper.
- [ ] Add tiny fixture tests for every facade method before wiring it into
      VNtyper.
- [ ] Add opt-in oracle tests against real CLI tools where useful.
- [ ] Add one end-to-end VNtyper native-facade test after each major backend is
      swapped in.

## Near-Term Order

- [x] Commit the Kestrel vendor/facade swap.
- [x] Add `vendor/rust/htslib-rs`.
- [x] Add `vendor/rust/bcftools-rs`.
- [x] Inspect `bcftools-rs` and `htslib-rs` APIs.
- [x] Implement the first `bcftools` native adapter method.
- [x] Add adapter tests for that method.
- [ ] Add `vendor/rust/samtools-rs` when ready.
- [ ] Implement the Samtools native facade methods needed for VNtyper.
- [ ] Add Samtools adapter tests using tiny BAM/FASTQ/depth fixtures.
- [ ] Refactor existing BioScript helper methods to call public facades.
- [ ] Build the minimal VNtyper BioScript pipeline on top of those facades.

## Verification Commands

```sh
cd rust
cargo test -p bioscript-libs -p bioscript-python -p bioscript-runtime
cargo test --manifest-path ../vendor/rust/kestrel-rs/Cargo.toml
```

```sh
PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'
PYTHONPATH=python:ports/vntyper/bioscript python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'
```
