# BioScript Native Library + VNtyper Port TODO

Goal: ship a BioScript version that includes the vendored native bioinformatics
libraries, preserves all existing BioScript behavior, and adds a VNtyper test
program ported to BioScript that passes parity tests comparable to upstream
VNtyper.

This is not just a facade spike. The finish line is:

- Existing BioScript scripts, runtime tests, Python wrapper tests, and Rust
  crate tests still pass.
- `vendor/rust` engines are wired through `bioscript-libs` and are the default
  native implementation path for the supported tool surfaces.
- A VNtyper BioScript program exists as the user-facing port, with the Python
  scaffold retained only as test/oracle support if still useful.
- VNtyper parity tests cover representative positive and negative samples,
  FASTQ and BAM entry points, report JSON, TSV calls, and HTML report structure.
- Any remaining gap against upstream VNtyper is documented with a concrete owner:
  BioScript runtime, `bioscript-libs`, `samtools-rs`, `bcftools-rs`,
  `kestrel-rs`, or VNtyper-port logic.

## Current Baseline

- [x] Vendored Rust engines exist under `vendor/rust`:
      `kestrel-rs`, `htslib-rs`, `bcftools-rs`, and `samtools-rs`.
- [x] Python reference libraries are kept under `vendor/python` where needed.
- [x] `rust/bioscript-libs` exposes recognizable facades for:
      `samtools`, `bcftools`, `kestrel`, `pysam`, `pyfaidx`, and VCF helpers.
- [x] `python/bioscript` exposes matching import names for Python-side tests and
      wrapper use.
- [x] `ports/vntyper/vntyper` contains the upstream VNtyper source as the
      reference implementation.
- [x] `ports/vntyper/test-data` contains ignored representative BAM/FASTQ data
      and expected output material.
- [x] `ports/vntyper/bioscript` contains the current Python-style VNtyper port
      scaffold and report logic.

## Non-Negotiable Gates

- [x] Establish one command that runs the old BioScript test suite.
      Suggested gate:
      `cd rust && CC=cc AR=ar cargo test --workspace`
      plus Python tests:
      `PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'`.
      Verified 2026-05-14:
      `CC=cc AR=ar cargo test --workspace` from `rust/` passes after restoring
      wasm compatibility with the current `VariantSpec` shape and report
      analysis visibility. The gate includes APOL1 real-file tests and the
      first-party Rust source-size guard.
      `PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'`
      passes: 31 tests, 2 skipped.
- [x] Establish one command that runs all BioScript facade tests against the
      vendored native engines.
      Suggested gate:
      `cd rust && CC=cc AR=ar cargo test -p bioscript-libs -p bioscript-python -p bioscript-runtime`.
      Verified 2026-05-14: passes. Coverage includes `bioscript-libs`,
      `bioscript-python`, and `bioscript-runtime` facade/runtime tests.
- [x] Establish one command that runs the VNtyper port tests that do not require
      large data or external tools.
      Suggested gate:
      `PYTHONPATH=python:ports/vntyper/bioscript python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'`.
      Verified 2026-05-14: 70 tests, 7 skipped. Skips are opt-in large-data or
      external-tool gates.
- [x] Establish opt-in commands for large-data parity gates:
      `BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1`,
      `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1`, and any new FASTQ/native parity gate.
      Documented in `docs/lib-support.md`. Added
      `BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1` via
      `ports/vntyper/tests/test_native_fastq_pipeline_gate.py`.
- [x] Add a short `docs/lib-support.md` or equivalent section documenting these
      gates so future work cannot silently regress the old BioScript behavior.
      See `docs/lib-support.md` "Verification Gates".

## Native Library Integration

- [x] Confirm `bioscript-libs` depends on vendored `kestrel-rs`, `htslib-rs`,
      `bcftools-rs`, and `samtools-rs` by local path or submodule revision.
      Confirmed in `rust/bioscript-libs/Cargo.toml`:
      `bcftools-rs`, `htslib-rs`, `kanalyze`, `kestrel`, and `samtools-rs`
      are all local paths under `vendor/rust`.
- [x] Add a dependency graph note in `docs/`:
      BioScript syntax/runtime -> `bioscript-libs` facade -> vendored engine.
      See `docs/lib-support.md` "Current Dependency Graph".
- [ ] Make native facades the default path for BioScript runtime calls where a
      native implementation exists.
- [ ] Keep command-builder fallbacks for dry-run/planning, but mark them as
      planning surfaces rather than the primary implementation.
- [x] Audit Python wrappers and runtime methods so supported names match:
      `from bioscript import samtools, bcftools, kestrel, pysam, pyfaidx`.
      Confirmed by `python/bioscript/__init__.py`, module wrapper tests, and
      runtime import tests for the supported names.
- [x] Add a test that imports each supported module from BioScript runtime syntax
      and verifies at least one method dispatch reaches the Rust facade.
      Existing runtime tests cover library imports, command builders, native
      Samtools/BCFtools materialization, Kestrel/VCF helpers, Pyfaidx aliasing,
      and Pysam fetch through runtime dispatch.
- [x] Add a test that imports each supported module from `python/bioscript` and
      verifies native extension delegation or a documented fallback.
      Existing Python tests cover backend policy, pure Python fallbacks, and
      native delegation for the supported wrappers.

## Existing BioScript Compatibility

- [x] Run all existing Rust tests before changing VNtyper behavior and save the
      command/output summary in this TODO.
      Verified 2026-05-14: `CC=cc AR=ar cargo test --workspace` passes from
      `rust/`.
- [x] Run all existing Python tests before changing VNtyper behavior and save the
      command/output summary in this TODO.
      Verified 2026-05-14:
      `PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'`
      passes: 31 tests, 2 skipped.
- [x] Run existing `bioscripts/` examples or their current tests if available.
      The Rust workspace gate includes CLI and APOL1 real-file tests:
      `tests/apol1_real_files.rs` and `tests/cli.rs` pass.
- [x] Keep APOL1/load-genotypes behavior unchanged unless a dedicated parity
      test proves the refactor is equivalent.
      No APOL1/load-genotypes refactor was made in this pass; existing APOL1
      tests pass under the Rust workspace gate.
- [ ] Add regression tests before replacing any old helper with a facade-backed
      implementation.
- [x] Check first-party production Rust source files under
      `rust/bioscript-*/src/**/*.rs` stay at or below 500 lines after edits.
      Verified by `bioscript-core/tests/source_size.rs` in the Rust workspace
      gate.

## VNtyper Program Shape

- [x] Decide the final user-facing program path.
      Proposed path: `ports/vntyper/bioscript/vntyper.bio` or
      `ports/vntyper/bioscript/vntyper.bs`.
      Decision: use `ports/vntyper/bioscript/vntyper.bs` for the final
      BioScript program. Documented in `ports/vntyper/bioscript/README.md`.
- [x] Keep `ports/vntyper/bioscript/vntyper.bs.py` only as an executable sketch
      until the real BioScript/Monty program can run.
      Documented in `ports/vntyper/bioscript/README.md`.
- [x] Define the public BioScript interface for VNtyper:
      input BAM or FASTQ pair, reference build, output directory, participant ID,
      optional report flags.
      Documented BAM and FASTQ entry points in
      `ports/vntyper/bioscript/README.md`.
- [ ] Port the current Python scaffold into actual BioScript syntax supported by
      the runtime.
      Initial command-planning program exists at
      `ports/vntyper/bioscript/vntyper.bs` and runs through the CLI. The native
      execution/post-processing pipeline still needs to move from the Python
      scaffold into runnable BioScript/runtime-supported calls.
- [ ] If Monty syntax is missing required features, add the smallest runtime or
      syntax support needed and cover it with runtime tests.
- [x] Keep VNtyper-specific constants in one config surface:
      MUC1 regions, reference FASTA path, Kestrel parameters, confidence
      thresholds, report fields, and optional adVNTR flags.
      `ports/vntyper/bioscript/vntyper_config.py` centralizes the current
      VNtyper-specific regions, reference paths, Kestrel parameters,
      thresholds, report keys, and optional-module toggles.
- [x] Keep the BioScript VNtyper program small: it should coordinate facades and
      call VNtyper-specific functions, not reimplement samtools/bcftools/kestrel
      internals.
      `vntyper.bs` and `vntyper-fastq.bs` are command-plan coordinator scripts;
      reusable tool behavior remains in `bioscript-libs` facades and vendored
      Rust engines.

## VNtyper Native Execution Path

- [x] BAM path:
      `samtools.view_region_native` -> `samtools.fastq_native` ->
      `samtools.depth_native` -> `kestrel.run_native` ->
      `bcftools.sort_native/index_native` -> VNtyper post-processing/report.
      Verified by the opt-in all-native BAM gate for representative positive
      and negative fixtures.
- [ ] FASTQ path:
      input FASTQ pair -> `kestrel.run_native` ->
      `bcftools.sort_native/index_native` -> VNtyper post-processing/report.
- [x] Ensure the BAM path can run without Java Kestrel, external samtools, or
      external bcftools when native gates are enabled.
      `require_all_native_bam_pipeline_prerequisites()` no longer requires
      Java or a Kestrel jar, and the all-native BAM parity test passed on
      2026-05-14 with `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1`.
- [x] Ensure the FASTQ path can run without Java Kestrel or external bcftools
      when native gates are enabled.
      Verified 2026-05-14 that the native FASTQ gate executes through native
      Kestrel and native BCFtools without Java/external tools. Parity is not
      yet correct: the negative fixture currently reports `High_Precision`
      instead of expected `negative`.
- [x] Add one CLI/runtime command that runs the BioScript VNtyper program against
      a BAM fixture.
      `vntyper_bioscript_program_runs_via_cli_and_writes_command_plan` runs
      `ports/vntyper/bioscript/vntyper.bs` with the representative positive BAM
      fixture and verifies the generated command plan.
- [x] Add one CLI/runtime command that runs the BioScript VNtyper program against
      a FASTQ fixture pair.
      Added `ports/vntyper/bioscript/vntyper-fastq.bs` and runtime coverage in
      `rust/bioscript-runtime/tests/vntyper_program.rs`.

## VNtyper Parity Tests

- [x] Inventory upstream VNtyper tests under
      `ports/vntyper/vntyper/tests` and map each relevant test to one of:
      port directly, replace with Rust facade test, replace with BioScript
      runtime test, or intentionally out of scope.
      See `ports/vntyper/tests/upstream-test-map.md`.
- [x] Create `ports/vntyper/tests/upstream-test-map.md` with that mapping.
- [x] Add unit tests for VNtyper-specific post-processing:
      VCF parsing, frameshift classification, depth score, confidence class,
      motif filtering, final best-call selection, TSV output, report JSON.
      Existing tests cover this in `test_vntyper_port.py`,
      `test_ported_upstream_units.py`, `test_upstream_scoring_parity.py`, and
      `test_vntyper_report.py`.
- [x] Add Rust tests where the behavior belongs in `bioscript-libs` rather than
      Python scaffolding.
      Candidate areas: VCF parsing, report-neutral call table generation,
      facade error mapping, and native command result shapes.
      Added `rust/bioscript-libs/tests/vntyper_facades.rs` for the native
      Samtools/Kestrel/BCFtools facade path on tiny generated fixtures. Existing
      `api.rs` tests cover VCF parsing and facade error mapping.
- [x] Add BioScript runtime tests that execute the VNtyper BioScript program on
      tiny deterministic fixtures.
      Added `rust/bioscript-runtime/tests/vntyper_program.rs`, which executes
      `ports/vntyper/bioscript/vntyper.bs` through `BioscriptRuntime` and
      verifies the generated command plan.
- [x] Add large-data opt-in parity tests for positive and negative BAM fixtures.
      Covered by `test_native_bam_pipeline_gate.py` and the existing external
      BAM gate.
- [x] Add large-data opt-in parity tests for positive and negative FASTQ
      fixtures.
      Added `test_native_fastq_pipeline_gate.py`, gated by
      `BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1`.
- [ ] Compare generated `kestrel_result.tsv` to expected fixture output.
- [ ] Compare generated `report.json` to expected fixture output, with explicit
      allowances for paths, timestamps, and tool-version metadata.
- [x] Compare generated HTML report structure against expected report content:
      summary, coverage QC, variant table, flags, pipeline log, and optional IGV
      configuration.
      `test_vntyper_report.py` covers generated report structure from fixture
      JSON/report rows, including summary, coverage QC, variant table controls,
      flags, pipeline log, and optional IGV configuration. Byte-for-byte
      upstream HTML parity is not available as an upstream fixture target.
- [ ] Make every large-data parity skip message list exactly which file, tool,
      environment variable, or native extension is missing.

## Engine Parity Gaps To Close Or Escalate

- [ ] `samtools-rs`: verify FASTQ extraction matches the VNtyper command chain
      `view -P | sort -n | fastq -1/-2/-0/-s` for representative fixtures.
- [ ] `samtools-rs`: if counts differ from real samtools, reduce to a small
      fixture and fix in the engine crate or document an intentional difference.
- [ ] `kestrel-rs`: run VNtyper FASTQ positive/negative fixtures and compare
      VCF records against Java Kestrel expected outputs.
      Attempted 2026-05-14 via
      `BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1 PYTHONPATH=python:ports/vntyper/bioscript python -m unittest ports.vntyper.tests.test_native_fastq_pipeline_gate`.
      The gate failed on the negative fixture: native Kestrel classification was
      `High_Precision`, expected VNtyper classification was `negative`.
- [ ] `kestrel-rs`: any Java parity gaps should be reduced into
      `vendor/rust/kestrel-rs` tests, not hidden in BioScript tests.
- [x] `bcftools-rs`: confirm the VNtyper-required sort/compress/index path is
      complete for all generated VCFs.
      Confirmed for tiny Kestrel-generated VCFs in
      `rust/bioscript-libs/tests/vntyper_facades.rs` and existing BCFtools
      adapter tests. Large-data generated VCF coverage remains part of the
      opt-in VNtyper parity gates.
- [ ] `bcftools-rs`: only implement native `view -i/-e` filtering if the
      BioScript VNtyper port actually needs it.
- [ ] `htslib-rs`: confirm shared BAM/CRAM/FASTA/VCF primitives are used through
      facades, not duplicated in BioScript-specific code.

## Rust Test Targets To Add

- [x] `rust/bioscript-libs/tests/vntyper_facades.rs`
      for native Samtools/Kestrel/BCFtools orchestration on tiny fixtures.
- [x] `rust/bioscript-libs/tests/vntyper_vcf.rs`
      for VNtyper-relevant VCF parsing and call-table conversion if moved to
      Rust.
- [x] `rust/bioscript-runtime/tests/vntyper_program.rs`
      for executing the BioScript VNtyper test program through the runtime.
- [ ] Keep large real-data tests opt-in and out of normal `cargo test` unless
      they use tiny checked-in fixtures.

## Python/Test Harness Work

- [x] Keep `ports/vntyper/tests/data_manifest.py` as the single source for
      large fixture paths and expected output paths.
      Existing large-data gates and manifest tests route through this helper.
- [x] Add FASTQ native prerequisites to the manifest, parallel to the existing
      BAM native prerequisites.
      Added `require_native_fastq_pipeline_prerequisites()` and
      `REPRESENTATIVE_FASTQ_CASES` in `ports/vntyper/tests/data_manifest.py`.
- [ ] Add or regenerate expected outputs for any checked-in representative
      FASTQ native fixtures.
- [x] Keep `ports/vntyper/test-data` ignored except for README/manifest files.
      Current git status shows no tracked test-data payload changes.
- [x] Remove generated `__pycache__` files from the repo if any are tracked.
      Verified with `git ls-files 'ports/vntyper/**/__pycache__/*'
      'python/**/__pycache__/*'`: no tracked generated cache files.
- [ ] Keep Python scaffold tests until equivalent Rust/BioScript runtime tests
      cover the behavior.

## Documentation

- [x] Document the supported BioScript imports and their backend engines.
      See `docs/lib-support.md`.
- [x] Document the VNtyper BioScript interface with one BAM example and one
      FASTQ example.
      See `ports/vntyper/bioscript/README.md`.
- [x] Document how to run small tests, full local tests, and opt-in large-data
      parity tests.
      See `docs/lib-support.md` and `ports/vntyper/bioscript/README.md`.
- [x] Document known gaps separately from TODO checkboxes once a gap is accepted
      as engine-owned or out of scope.
      See `ports/vntyper/tests/upstream-test-map.md`.

## Completion Criteria

- [x] Old BioScript Rust test gate passes.
      Verified 2026-05-14 with `CC=cc AR=ar cargo test --workspace`.
- [x] Old BioScript Python test gate passes.
      Verified 2026-05-14 with
      `PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'`.
- [x] Native facade Rust/Python tests pass.
      Verified 2026-05-14 with
      `CC=cc AR=ar cargo test -p bioscript-libs -p bioscript-python -p bioscript-runtime`
      and Python wrapper tests.
- [x] VNtyper small fixture tests pass without external Java/samtools/bcftools.
      Verified 2026-05-14 with
      `PYTHONPATH=python:ports/vntyper/bioscript python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'`.
- [x] VNtyper BAM positive/negative native parity gate passes.
      Verified 2026-05-14:
      `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 PYTHONPATH=python:ports/vntyper/bioscript python -m unittest ports.vntyper.tests.test_native_bam_pipeline_gate.VntyperNativeBamPipelineGateTests.test_native_bam_pipeline_with_native_kestrel_and_bcftools_matches_expected_classification`
      passed in 91.426s.
- [ ] VNtyper FASTQ positive/negative native parity gate passes.
      Current status 2026-05-14: gate runs but fails negative-fixture parity
      (`High_Precision` vs expected `negative`).
- [ ] VNtyper report JSON and TSV outputs match expected fixtures with explicit
      normalized fields.
- [ ] VNtyper HTML report structure test passes.
- [ ] Upstream VNtyper test map is complete and every relevant upstream behavior
      has a ported test, Rust facade test, runtime test, or documented exclusion.
- [ ] `TODO.md` contains no ambiguous "done enough" items; each completed item
      points to a file, test, command, or documented decision.
