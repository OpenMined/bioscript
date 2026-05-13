# BioScript Library Support TODO

Goal: make BioScript support recognizable bioinformatics library/tool surfaces
through thin `bioscript-libs` facades backed by vendored Rust engine crates.
Build the reusable primitives first, wire Samtools next, and then make the
VNtyper BioScript port a small amount of pipeline code plus data/config that
uses those built-in primitives.

## Direction

- [x] Use explicit BioScript imports:
      `from bioscript import samtools, bcftools, kestrel, pysam, pyfaidx`.
- [x] Treat BioScript library support as the product:
      common pipeline code should read like standard bioinformatics workflows,
      not like private BioScript internals.
- [x] Build in layers:
      engine crates -> BioScript facades -> facade tests -> VNtyper port.
      Current layering is engine crates under `vendor/rust`, public facades in
      `rust/bioscript-libs` plus `python/bioscript`, adapter/runtime tests, and
      VNtyper pipeline code under `ports/vntyper/bioscript`.
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
- [x] Add Samtools Rust engine:
      `vendor/rust/samtools-rs` from
      `git@github.com:madhavajay/samtools-rs.git`.
      The repo contains the VNtyper-needed `view`, `fastq`, `depth`, `index`,
      and related API surface.
- [x] Keep vendored engine crate tests inside their own repos/workspaces.
      `kestrel-rs`, `samtools-rs`, `bcftools-rs`, and `htslib-rs` keep their
      engine tests under their own vendored workspaces; BioScript only points
      at the submodule revisions and calls their public APIs.
- [x] Keep BioScript tests focused on adapter behavior and pipeline integration.
      BioScript-owned tests now cover argument normalization, runtime/Python
      wrappers, tiny fixture adapters, and VNtyper integration gates rather
      than re-testing whole engines.

## Rust Crate Wiring

- [x] Wire `rust/bioscript-libs` to local `kestrel-rs` path dependencies:
      `kestrel` and `kanalyze`.
- [x] Wire `rust/bioscript-libs` to local `htslib-rs`.
      The top-level submodule and the nested BCFtools HTS backend are advanced
      to `2f63d19` on `bioscript-samtools-template-fastq`, which includes the
      Samtools-native support and template-expanded BAM region writer needed by
      `samtools-rs`.
- [x] Wire `rust/bioscript-libs` to local `bcftools-rs`.
- [x] Wire `rust/bioscript-libs` to local `samtools-rs`.
      `bioscript-libs` depends on
      `vendor/rust/samtools-rs/crates/samtools-rs`, and the vendored
      `samtools-rs` workspace is patched on
      `bioscript-use-shared-htslib` to share the BCFtools HTS backend path so
      Cargo has one unambiguous `htslib-rs` package.
- [x] Add `[patch]` entries only where nested crate dependencies would
      otherwise pull remote git/crates.io versions instead of local submodules.
      No new engine-crate patches were needed: `bioscript-libs` uses path
      dependencies and the vendored `samtools-rs` workspace points at the
      shared nested `bcftools-rs/htslib-rs` path. Existing workspace patches
      remain limited to the local noodles/lexical overrides.
- [x] Document the dependency graph:
      BioScript -> `bioscript-libs` facade -> vendored Rust engine crate.

## Crate Publishing

- [x] Keep local path dependencies while `kestrel-rs`, `htslib-rs`,
      `bcftools-rs`, and `samtools-rs` APIs are still changing quickly.
- [ ] Publish those engine crates once their public APIs and test suites are
      stable enough for external consumers.
- [ ] After publishing, replace stable path dependencies with versioned crates
      where that simplifies the Cargo graph.
- [x] Keep submodules available for upstream test fixtures, source comparison,
      and local patching even after published crates are used by default.

## Milestones

- [x] M1: Kestrel Rust engine is vendored and callable through BioScript.
- [x] M2: HTS and BCFtools Rust engines are vendored and wired by path.
- [x] M3: Samtools Rust engine is vendored and wired by path.
- [x] M4: BioScript facades expose a minimal, recognizable built-in toolkit:
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
- [x] Move any remaining Kestrel algorithm parity expectations into
      `vendor/rust/kestrel-rs`.
      Java/Rust parity and algorithm behavior tests live in the Kestrel engine
      workspace, including `crates/kestrel/tests/cli_parity.rs` and the
      Java-compatible unit tests. BioScript keeps only facade smoke coverage.

## Samtools Facade

- [x] Existing BioScript command-builder surface:
      `samtools.view_region`, `samtools.fastq`, `samtools.depth`.
- [x] Existing native prototype supports BAM slicing, FASTQ extraction, and
      depth summary through BioScript-owned primitives.
- [x] Replace native prototype internals with calls into `samtools-rs`.
      `view_region_native`, `fastq_native`, and `depth_native` now call
      `samtools_rs::native` and adapt the results back to BioScript's existing
      return shapes.
- [x] Prioritize Samtools now that `samtools-rs` is available because VNtyper's BAM
      path should become:
      `samtools.view` -> `samtools.index/sort` if needed ->
      `samtools.fastq` -> `samtools.depth`.
- [x] Keep the public BioScript API shaped like familiar samtools operations:
      `view`, `fastq`, `sort`, `index`, `depth`, `faidx`.
      Command-builder facades for those names are exposed in Rust, Python, and
      the runtime; VNtyper-specific template extraction stays in the native
      `fastq_native` adapter.
- [x] Add adapter tests for:
      region parsing, indexed BAM input, `.bam/.bai` discovery, paired FASTQ
      output counts, depth summary fields, and error mapping.
      Covered by `samtools_native_adapter_handles_tiny_indexed_bam`, which
      creates a tiny SAM/BAM fixture in a temp dir and exercises the BioScript
      Samtools facade end to end.
- [x] Keep oracle tests against real samtools opt-in only.
      `test_samtools_fastq_oracle.py` is gated by
      `BIOSCRIPT_RUN_SAMTOOLS_ORACLE=1` and external samtools availability.

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
- [x] Initial target operations:
      `view`, `sort`, `norm`, compression/index helpers if needed.
      Command-builder facades now cover `view`, `sort`, `norm`,
      `view_filter`, and `index`; native helpers cover `view`, `sort`, and
      indexing where `bcftools-rs` already supports them.
- [ ] Add adapter tests for VCF input/output, compressed output, filter
      expressions used by VNtyper, and useful error messages.
      Initial coverage verifies `bcftools-rs` header extraction, VCF output,
      BGZF-compressed output, native sort, CSI/TBI indexing, Python wrapper
      delegation, malformed-input error propagation, and the real PyO3 native
      extension when installed. Filter expression coverage at the command-builder
      layer exists; native filter expression coverage remains pending until
      `bcftools-rs view` supports `-i/-e`.

## HTS / Pysam / Pyfaidx Facades

- [x] Keep `pysam` and `pyfaidx` as recognizable compatibility namespaces.
- [x] `pyfaidx.Fasta` has a small Rust/Python-compatible FASTA slice surface.
- [x] `pysam.AlignmentFile.fetch` has initial BAM/CRAM read support.
- [ ] Refactor lower-level alignment code to flow through `pysam` or
      `samtools` facades where that makes scripts more recognizable.
- [x] Use `htslib-rs` as the shared backend for BAM/CRAM/VCF/FASTA primitives
      once vendored.
      FASTA access in `bioscript-libs` `pyfaidx` now builds and queries
      through `htslib_rs::faidx_compat`; Samtools/BCFtools already enter via
      their vendored engine crates. The pysam-style BAM/CRAM fetch path now
      routes through `htslib_rs::alignment_compat` indexed query helpers and
      converts HTS records into the BioScript `AlignedSegment` surface.
- [ ] Add parity tests from focused upstream `pysam` and `pyfaidx` cases, not
      the full upstream test suites.

## Python Package

- [x] Keep top-level `python/bioscript` matching BioScript import names.
- [x] Keep optional delegation to real Python libraries where useful.
- [x] Expose native functions through `rust/bioscript-python`.
- [x] Add Python tests that call the real native extension for each engine
      facade with tiny fixtures.
      `python/tests/test_tools.py` now exercises real `_native` calls for
      Kestrel, Samtools, and BCFtools. `pyfaidx` now has a Rust-backend Python
      wrapper around `pyfaidx_fetch_native` with mocked-extension coverage and
      `bioscript-python` compile coverage; `pysam` remains documented as a
      pending Python native facade.
- [x] Keep mocked-extension tests for argument normalization and missing-native
      behavior.
- [x] Make Python-only fallback behavior explicit per module:
      real Python library, pure Python fallback, or native-required.

## Runtime / Monty Integration

- [x] Support `from bioscript import x` import rewriting for current modules.
- [x] Bind initial module objects and method calls in `bioscript-runtime`.
- [x] Add runtime method bindings for native samtools/bcftools operations once
      facades are stable.
      BCFtools native bindings now cover `view_header_native`, `view_native`,
      `sort_native`, and `index_native`; Samtools native bindings now cover
      `view_region_native`, `fastq_native`, and `depth_native` through the
      BioScript facade, which is backed by `samtools-rs`.
- [x] Keep runtime responsible for language/object adaptation only.
      Runtime methods now adapt Monty objects, paths, and return shapes while
      delegating tool behavior to `bioscript-libs` facades.
- [x] Keep file/path/security policy centralized and reused across facades.
      Native Samtools and BCFtools runtime bindings use the same
      `resolve_existing_user_path` / `resolve_user_write_path` sandbox checks
      as other host-facing methods, with security tests covering materialized
      outputs.

## VNtyper Proof Port

- [x] Keep upstream VNtyper source vendored at `ports/vntyper/vntyper`.
- [x] Keep local large test data ignored under `ports/vntyper/test-data`.
- [x] Keep BioScript VNtyper port under `ports/vntyper/bioscript`.
- [x] Keep BioScript-owned VNtyper tests under `ports/vntyper/tests`.
- [x] Current tests cover command planning, Kestrel VCF parsing, scoring,
      report JSON/HTML shape, and fake-runner pipeline behavior.
- [x] Current adapter smoke tests prove BioScript can call `kestrel-rs`.
- [x] Reframe the final VNtyper port as its own BioScript code, not as a copy
      of every upstream dependency. The VNtyper-specific layer should contain:
      MUC1 regions, motif/reference data, Kestrel parameter choices,
      frameshift/depth classification, report rows, and CLI/pipeline glue.
- [x] Keep generic work out of the VNtyper port. Generic work belongs in
      BioScript facades:
      BAM/CRAM slicing, FASTQ extraction, depth, VCF parsing/filtering,
      Kestrel calling, FASTA lookup, TSV/JSON helpers.
- [x] Refactor VNtyper pipeline code to prefer:
      `samtools.*`, `bcftools.*`, `kestrel.*`, `pysam.*`, and `pyfaidx.*`
      over private helper names.
      `ports/vntyper/bioscript/vntyper_commands.py` builds the BAM plan
      through `bioscript.samtools`, `bioscript.bcftools`, and
      `bioscript.kestrel`; `vntyper_external_pipeline.py` uses the same public
      facade modules for native Samtools, Kestrel, and BCFtools execution.
      Native Kestrel execution now goes through `kestrel.run_native(...)`
      instead of VNtyper manually loading references and writing VCF text.
      The FASTQ-only path can now optionally run native Kestrel followed by
      native BCFtools sort/index without Java or external bcftools.
      The BAM path also has a native BCFtools sort/index switch, so native or
      external Kestrel output can be materialized as sorted/indexed VCF through
      the same `bcftools.sort_native(...)` facade.
- [x] Define the minimal VNtyper BioScript interface, for example:
      `run_vntyper(bam=..., reference_build="hg19", output_dir=...)` and
      `run_vntyper_fastq(r1=..., r2=..., reference_build="hg19", output_dir=...)`.
- [x] Keep VNtyper data/config small and explicit:
      MUC1 coordinates, motif FASTA path, confidence thresholds, report schema,
      and optional validation toggles.
      `ports/vntyper/bioscript/vntyper_config.py` centralizes the MUC1
      GRCh37/GRCh38 regions, motif FASTA path, Kestrel thresholds, report
      schema keys, native Kestrel bounds, and disabled-by-default adVNTR toggle.
      `ports/vntyper/tests/test_vntyper_config.py` guards that the explicit
      config matches the generated report surface.
- [x] Now that `samtools-rs` and `bcftools-rs` are wired, rerun the BAM path using
      only BioScript native facades.
      Verified the opt-in native-Samtools BAM gate with Java Kestrel, native
      Kestrel, and the all-native native-Samtools/native-Kestrel/native-BCFtools
      path for the representative positive and negative fixtures. The all-native
      gate now asserts matching Kestrel classification, matching screening
      summary, and creation of the native BCFtools sorted VCF plus CSI index.
- [ ] Compare native-facade VNtyper output against expected positive/negative
      fixtures for:
      FASTQ path, BAM path, report JSON, and HTML report.
      BAM report JSON/classification parity is covered by the opt-in all-native
      gate. FASTQ native parity and HTML report comparisons remain open.
- [x] Keep large real-data parity tests opt-in with clear skip messages.
      Large VNtyper data gates live behind explicit environment switches such
      as `BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1`,
      `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1`, and
      `BIOSCRIPT_RUN_SAMTOOLS_ORACLE=1`; missing data, tools, expected
      outputs, and native extensions raise `unittest.SkipTest` with concrete
      prerequisite messages.

## Test Policy

- [x] Engine crates own engine correctness:
      e.g. `vendor/rust/kestrel-rs` owns Kestrel Java/algorithm parity.
- [x] BioScript owns facade correctness:
      argument normalization, path handling, output shape, error mapping, and
      integration with BioScript/Python/VNtyper.
- [x] Add tiny fixture tests for every facade method before wiring it into
      VNtyper.
      Coverage now spans Samtools, BCFtools, Kestrel, pysam, pyfaidx, VCF/table
      helpers, Python wrapper delegation, and runtime imports/materialization.
- [x] Add opt-in oracle tests against real CLI tools where useful.
      Real-tool gates are opt-in, including the Samtools FASTQ oracle and
      VNtyper external/native BAM gates.
- [x] Add one end-to-end VNtyper native-facade test after each major backend is
      swapped in.
      `test_native_bam_pipeline_gate.py` exercises the native Samtools facade
      with the VNtyper BAM path, then native Kestrel, then the all-native
      native-Samtools/native-Kestrel/native-BCFtools path for representative
      positive and negative fixtures.

## Near-Term Order

- [x] Commit the Kestrel vendor/facade swap.
- [x] Add `vendor/rust/htslib-rs`.
- [x] Add `vendor/rust/bcftools-rs`.
- [x] Inspect `bcftools-rs` and `htslib-rs` APIs.
- [x] Implement the first `bcftools` native adapter method.
- [x] Add adapter tests for that method.
- [x] Add `vendor/rust/samtools-rs` from
      `git@github.com:madhavajay/samtools-rs.git`.
      The stale local config/worktree state was reused with the SSH remote.
- [x] Implement the Samtools native facade methods needed for VNtyper.
      `view_region_native`, `fastq_native`, and `depth_native` are backed by
      `samtools-rs`; native `index/sort` can be exposed later if VNtyper needs
      them after BAM slicing.
- [x] Add Samtools adapter tests using tiny BAM/FASTQ/depth fixtures.
      `samtools_native_adapter_handles_tiny_indexed_bam` writes a tiny SAM
      fixture, converts it to BAM, indexes it, and checks native view, FASTQ,
      depth, and error behavior through the BioScript facade. `samtools-rs`
      owns broader command/native-wrapper engine tests.
      Opt-in oracle testing against real `samtools fastq` is close but not
      exact yet: the native path currently emits +20 read1 records on the
      positive fixture and +3 on the negative fixture versus real samtools.
      Keep this open until `samtools-rs` fully matches `view -P | sort -n |
      fastq -1/-2/-0/-s` behavior.
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
