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

- [ ] Establish one command that runs the old BioScript test suite.
      Suggested gate:
      `cd rust && CC=cc AR=ar cargo test --workspace`
      plus Python tests:
      `PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'`.
- [ ] Establish one command that runs all BioScript facade tests against the
      vendored native engines.
      Suggested gate:
      `cd rust && CC=cc AR=ar cargo test -p bioscript-libs -p bioscript-python -p bioscript-runtime`.
- [ ] Establish one command that runs the VNtyper port tests that do not require
      large data or external tools.
      Suggested gate:
      `PYTHONPATH=python:ports/vntyper/bioscript python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'`.
- [ ] Establish opt-in commands for large-data parity gates:
      `BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1`,
      `BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1`, and any new FASTQ/native parity gate.
- [ ] Add a short `docs/lib-support.md` or equivalent section documenting these
      gates so future work cannot silently regress the old BioScript behavior.

## Native Library Integration

- [ ] Confirm `bioscript-libs` depends on vendored `kestrel-rs`, `htslib-rs`,
      `bcftools-rs`, and `samtools-rs` by local path or submodule revision.
- [ ] Add a dependency graph note in `docs/`:
      BioScript syntax/runtime -> `bioscript-libs` facade -> vendored engine.
- [ ] Make native facades the default path for BioScript runtime calls where a
      native implementation exists.
- [ ] Keep command-builder fallbacks for dry-run/planning, but mark them as
      planning surfaces rather than the primary implementation.
- [ ] Audit Python wrappers and runtime methods so supported names match:
      `from bioscript import samtools, bcftools, kestrel, pysam, pyfaidx`.
- [ ] Add a test that imports each supported module from BioScript runtime syntax
      and verifies at least one method dispatch reaches the Rust facade.
- [ ] Add a test that imports each supported module from `python/bioscript` and
      verifies native extension delegation or a documented fallback.

## Existing BioScript Compatibility

- [ ] Run all existing Rust tests before changing VNtyper behavior and save the
      command/output summary in this TODO.
- [ ] Run all existing Python tests before changing VNtyper behavior and save the
      command/output summary in this TODO.
- [ ] Run existing `bioscripts/` examples or their current tests if available.
- [ ] Keep APOL1/load-genotypes behavior unchanged unless a dedicated parity
      test proves the refactor is equivalent.
- [ ] Add regression tests before replacing any old helper with a facade-backed
      implementation.
- [ ] Check first-party production Rust source files under
      `rust/bioscript-*/src/**/*.rs` stay at or below 500 lines after edits.

## VNtyper Program Shape

- [ ] Decide the final user-facing program path.
      Proposed path: `ports/vntyper/bioscript/vntyper.bio` or
      `ports/vntyper/bioscript/vntyper.bs`.
- [ ] Keep `ports/vntyper/bioscript/vntyper.bs.py` only as an executable sketch
      until the real BioScript/Monty program can run.
- [ ] Define the public BioScript interface for VNtyper:
      input BAM or FASTQ pair, reference build, output directory, participant ID,
      optional report flags.
- [ ] Port the current Python scaffold into actual BioScript syntax supported by
      the runtime.
- [ ] If Monty syntax is missing required features, add the smallest runtime or
      syntax support needed and cover it with runtime tests.
- [ ] Keep VNtyper-specific constants in one config surface:
      MUC1 regions, reference FASTA path, Kestrel parameters, confidence
      thresholds, report fields, and optional adVNTR flags.
- [ ] Keep the BioScript VNtyper program small: it should coordinate facades and
      call VNtyper-specific functions, not reimplement samtools/bcftools/kestrel
      internals.

## VNtyper Native Execution Path

- [ ] BAM path:
      `samtools.view_region_native` -> `samtools.fastq_native` ->
      `samtools.depth_native` -> `kestrel.run_native` ->
      `bcftools.sort_native/index_native` -> VNtyper post-processing/report.
- [ ] FASTQ path:
      input FASTQ pair -> `kestrel.run_native` ->
      `bcftools.sort_native/index_native` -> VNtyper post-processing/report.
- [ ] Ensure the BAM path can run without Java Kestrel, external samtools, or
      external bcftools when native gates are enabled.
- [ ] Ensure the FASTQ path can run without Java Kestrel or external bcftools
      when native gates are enabled.
- [ ] Add one CLI/runtime command that runs the BioScript VNtyper program against
      a BAM fixture.
- [ ] Add one CLI/runtime command that runs the BioScript VNtyper program against
      a FASTQ fixture pair.

## VNtyper Parity Tests

- [ ] Inventory upstream VNtyper tests under
      `ports/vntyper/vntyper/tests` and map each relevant test to one of:
      port directly, replace with Rust facade test, replace with BioScript
      runtime test, or intentionally out of scope.
- [ ] Create `ports/vntyper/tests/upstream-test-map.md` with that mapping.
- [ ] Add unit tests for VNtyper-specific post-processing:
      VCF parsing, frameshift classification, depth score, confidence class,
      motif filtering, final best-call selection, TSV output, report JSON.
- [ ] Add Rust tests where the behavior belongs in `bioscript-libs` rather than
      Python scaffolding.
      Candidate areas: VCF parsing, report-neutral call table generation,
      facade error mapping, and native command result shapes.
- [ ] Add BioScript runtime tests that execute the VNtyper BioScript program on
      tiny deterministic fixtures.
- [ ] Add large-data opt-in parity tests for positive and negative BAM fixtures.
- [ ] Add large-data opt-in parity tests for positive and negative FASTQ
      fixtures.
- [ ] Compare generated `kestrel_result.tsv` to expected fixture output.
- [ ] Compare generated `report.json` to expected fixture output, with explicit
      allowances for paths, timestamps, and tool-version metadata.
- [ ] Compare generated HTML report structure against expected report content:
      summary, coverage QC, variant table, flags, pipeline log, and optional IGV
      configuration.
- [ ] Make every large-data parity skip message list exactly which file, tool,
      environment variable, or native extension is missing.

## Engine Parity Gaps To Close Or Escalate

- [ ] `samtools-rs`: verify FASTQ extraction matches the VNtyper command chain
      `view -P | sort -n | fastq -1/-2/-0/-s` for representative fixtures.
- [ ] `samtools-rs`: if counts differ from real samtools, reduce to a small
      fixture and fix in the engine crate or document an intentional difference.
- [ ] `kestrel-rs`: run VNtyper FASTQ positive/negative fixtures and compare
      VCF records against Java Kestrel expected outputs.
- [ ] `kestrel-rs`: any Java parity gaps should be reduced into
      `vendor/rust/kestrel-rs` tests, not hidden in BioScript tests.
- [ ] `bcftools-rs`: confirm the VNtyper-required sort/compress/index path is
      complete for all generated VCFs.
- [ ] `bcftools-rs`: only implement native `view -i/-e` filtering if the
      BioScript VNtyper port actually needs it.
- [ ] `htslib-rs`: confirm shared BAM/CRAM/FASTA/VCF primitives are used through
      facades, not duplicated in BioScript-specific code.

## Rust Test Targets To Add

- [ ] `rust/bioscript-libs/tests/vntyper_facades.rs`
      for native Samtools/Kestrel/BCFtools orchestration on tiny fixtures.
- [ ] `rust/bioscript-libs/tests/vntyper_vcf.rs`
      for VNtyper-relevant VCF parsing and call-table conversion if moved to
      Rust.
- [ ] `rust/bioscript-runtime/tests/vntyper_program.rs`
      for executing the BioScript VNtyper test program through the runtime.
- [ ] Keep large real-data tests opt-in and out of normal `cargo test` unless
      they use tiny checked-in fixtures.

## Python/Test Harness Work

- [ ] Keep `ports/vntyper/tests/data_manifest.py` as the single source for
      large fixture paths and expected output paths.
- [ ] Add FASTQ native prerequisites to the manifest, parallel to the existing
      BAM native prerequisites.
- [ ] Add or regenerate expected outputs for any checked-in representative
      FASTQ native fixtures.
- [ ] Keep `ports/vntyper/test-data` ignored except for README/manifest files.
- [ ] Remove generated `__pycache__` files from the repo if any are tracked.
- [ ] Keep Python scaffold tests until equivalent Rust/BioScript runtime tests
      cover the behavior.

## Documentation

- [ ] Document the supported BioScript imports and their backend engines.
- [ ] Document the VNtyper BioScript interface with one BAM example and one
      FASTQ example.
- [ ] Document how to run small tests, full local tests, and opt-in large-data
      parity tests.
- [ ] Document known gaps separately from TODO checkboxes once a gap is accepted
      as engine-owned or out of scope.

## Completion Criteria

- [ ] Old BioScript Rust test gate passes.
- [ ] Old BioScript Python test gate passes.
- [ ] Native facade Rust/Python tests pass.
- [ ] VNtyper small fixture tests pass without external Java/samtools/bcftools.
- [ ] VNtyper BAM positive/negative native parity gate passes.
- [ ] VNtyper FASTQ positive/negative native parity gate passes.
- [ ] VNtyper report JSON and TSV outputs match expected fixtures with explicit
      normalized fields.
- [ ] VNtyper HTML report structure test passes.
- [ ] Upstream VNtyper test map is complete and every relevant upstream behavior
      has a ported test, Rust facade test, runtime test, or documented exclusion.
- [ ] `TODO.md` contains no ambiguous "done enough" items; each completed item
      points to a file, test, command, or documented decision.
