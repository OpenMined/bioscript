# BioScript Library Support

BioScript should support standard bioinformatics workflows through a small set of
Python-like library shims backed by Rust native code. The first target syntax is:

```python
from bioscript import pysam
```

This makes the compatibility boundary explicit. The imported `pysam` module is a
BioScript-supported subset that mimics the real Python `pysam` API where useful;
it is not a promise that every Python import or every upstream `pysam` feature is
available inside BioScript.

## Goals

- Let assay and pipeline authors write familiar bioinformatics code.
- Keep BioScript execution fast, deterministic, and sandboxable.
- Back common APIs with Rust crates such as `noodles`.
- Share the same shim surface between BioScript/Monty and normal Python tests.
- Use upstream source and tests to guide compatibility, without committing to
  full-library parity up front.

## Verification Gates

Use these gates when changing BioScript library support or the VNtyper port.
They are intentionally split so old BioScript compatibility, native facade
coverage, and VNtyper-port behavior can be diagnosed separately.

```sh
cd rust
CC=cc AR=ar cargo test --workspace
```

This is the old BioScript compatibility gate. It covers the Rust workspace,
CLI tests, APOL1 real-file tests, runtime tests, wasm tests, facade tests, and
the first-party Rust source-size guard.

```sh
PYTHONPATH=python python -m unittest discover -s python/tests -p 'test_*.py'
```

This is the Python wrapper gate. It verifies the `python/bioscript` package,
backend policies, pure-Python fallbacks, and native-extension delegation tests.

```sh
cd rust
CC=cc AR=ar cargo test -p bioscript-libs -p bioscript-python -p bioscript-runtime
```

This is the focused native facade gate. It verifies the Rust facades, PyO3
native extension crate, and BioScript runtime dispatch for supported library
imports.

```sh
PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'
```

This is the small VNtyper-port gate. It excludes large data unless the opt-in
environment variables are set, but still covers command planning,
post-processing, report rendering, upstream-unit ports, and skip behavior.

Large-data parity gates are opt-in and must keep clear skip messages:

```sh
BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_full_pipeline_gate

BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_native_bam_pipeline_gate

BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 \
BIOSCRIPT_RUN_NATIVE_BAM_OUTPUT_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_native_bam_pipeline_gate.VntyperNativeBamPipelineGateTests.test_native_bam_output_fingerprints_match_expected_outputs

BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_native_fastq_pipeline_gate

BIOSCRIPT_RUN_SAMTOOLS_ORACLE=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_samtools_fastq_oracle
```

The native BAM gate checks classification/report-shape parity. The separate
native BAM output gate compares normalized `kestrel_result.tsv` fingerprints
and normalized report summaries, so it should remain opt-in until Kestrel
record-level parity is fixed or a field-level allowance is explicitly accepted.

## Proposed Stack

```text
BioScript source
  from bioscript import pysam
        |
        v
bioscript-runtime import binding
        |
        v
bioscript-libs module registry
        |
        v
pysam-compatible shim
        |
        v
Rust backends: bioscript-formats, noodles, vendored Rust engine crates
```

The runtime should only own language binding and object adaptation. The
bioinformatics API behavior should live in `bioscript-libs` so it can be reused
from the CLI, wasm, FFI, and a Python package.

## First Supported Syntax

Initial BioScript support should be narrow:

```python
from bioscript import pysam
from bioscript import pyfaidx
```

Later, if Monty import support matures, this can expand to:

```python
from bioscript import pysam as ps
import pysam
```

The plain `import pysam` form should be treated as optional compatibility sugar,
because it may conflict with real Python packages and implies broader Python
module resolution than BioScript needs at first.

## Folder Plan

```text
rust/
  bioscript-libs/
    Cargo.toml
    src/
      lib.rs
      module_registry.rs
      errors.rs
      value.rs
      pysam/
        mod.rs
        alignment_file.rs
        aligned_segment.rs
        pileup.rs
      pyfaidx/
        mod.rs
        fasta.rs
      vcf/
        mod.rs
        variant_file.rs
        record.rs

  bioscript-runtime/
    src/runtime/imports.rs
    src/runtime/modules.rs

  bioscript-python/
    Cargo.toml
    src/lib.rs

python/
  pyproject.toml
  bioscript/
    __init__.py
    pysam.py
    pyfaidx.py
    vcf.py
    _native.py
  tests/
    test_pysam_subset.py
    test_pyfaidx_subset.py
    test_runtime_parity.py

vendor/
  python/
    pysam/
    pyfaidx/
  rust/
    kestrel-rs/
    bcftools-rs/
    htslib-rs/
  testdata/
    pysam/
    samtools/
    vcf/
    fasta/
```

## Rust Crate Responsibilities

`bioscript-libs` owns the compatibility APIs:

- module registry for supported shim modules
- Rust-native objects that model selected external APIs
- conversion-neutral data structures that the runtime and Python bindings can
  adapt into their own object models
- compatibility errors with clear unsupported-feature messages

`bioscript-runtime` owns Monty integration:

- parsing or intercepting `from bioscript import <module>`
- binding a supported shim module to the local BioScript name
- dispatching method calls on shim objects into `bioscript-libs`
- enforcing runtime path, resource, and sandbox rules

`bioscript-python` and `python/bioscript` expose the same import names in
CPython:

- expose native helpers through `bioscript._native` when available
- optionally compare against real Python libraries during tests
- let authors run the same scripts in normal Python before running them in
  BioScript

Each Python shim exposes `BACKEND_POLICY`, a `ModuleBackendPolicy` with `auto`,
`python`, and `rust` descriptions. This makes fallback behavior explicit:

- `bcftools` in BioScript runtime: `sort`, `index`, and `view` default to
  native Rust execution; `plan_sort`, `plan_index`, and `plan_view` are the
  command-planning surfaces for dry runs and audit logs.
- `samtools` in BioScript runtime: `view`, `depth`, `sort`, and `index`
  default to native Rust execution where the native signature matches the
  public method. `view_region` defaults to native execution when
  `include_unmapped=False`; use `plan_view_region` for the command-plan form
  or for `include_unmapped=True`. `fastq` defaults to whole-BAM native FASTQ
  conversion; VNtyper's region-sliced FASTQ path still uses
  `fastq_native(bam, region, ...)`. `plan_*` remains the command-planning
  surface.
- `kestrel`: `plan_command` is the command-planning surface. Native execution
  uses `kestrel.run_native(reference_fasta, fastq_paths, output_vcf, ...)` in
  the BioScript runtime, or lower-level Python helpers such as `run_native` /
  `call_fastq_references_native`.
- Python shim modules keep command-builder compatibility names for now; native
  helpers still require `bioscript._native`.
- `pyfaidx`: auto mode uses real `pyfaidx` when installed, otherwise a small
  pure Python FASTA fallback; rust mode is pending
- `pysam`: auto/python mode delegates to real `pysam` when installed; rust mode
  is pending

## Current Dependency Graph

The graph should stay narrow: BioScript owns language/runtime adaptation,
`bioscript-libs` owns compatibility facades, and vendored Rust engine crates own
native bioinformatics behavior.

```text
BioScript source
  -> bioscript-runtime import/method binding
  -> bioscript-libs facade module
  -> vendored Rust engine crate
  -> lower-level format/statistics crates as needed
```

Current wired paths:

```text
from bioscript import kestrel
  -> bioscript-runtime KestrelModule or python/bioscript/kestrel.py
  -> rust/bioscript-libs::kestrel
  -> vendor/rust/kestrel-rs/crates/kestrel
  -> vendor/rust/kestrel-rs/crates/kanalyze

from bioscript import bcftools
  -> bioscript-runtime BcftoolsModule or python/bioscript/bcftools.py
  -> rust/bioscript-libs::bcftools
  -> vendor/rust/bcftools-rs/crates/bcftools-rs
  -> vendor/rust/bcftools-rs/htslib-rs

from bioscript import pysam / samtools / pyfaidx
  -> bioscript-runtime module binding or python/bioscript module
  -> rust/bioscript-libs facade
  -> vendor/rust/samtools-rs for samtools operations
  -> vendor/rust/bcftools-rs/htslib-rs for shared HTS-backed primitives
  -> bioscript-formats and noodles where BioScript owns the domain helper
```

Python authors can call low-level `kestrel.call_*_native(...)` helpers when
they need VCF text, or `kestrel.run_native(reference_fasta, fastq_paths,
output_vcf)` when a pipeline wants Kestrel-like file output.
BioScript runtime code can use `vcf.read_vntyper_kestrel(path)` for the current
VNtyper Kestrel call-table rows used by the FASTQ runtime slice, then
`vcf.build_vntyper_report_json(sample_name, input_files, rows)` for the current
VNtyper JSON report payload. The runtime facade also accepts optional metadata
and coverage dictionaries as fourth and fifth arguments for BAM-style reports.

When `kestrel-rs`, `bcftools-rs`, `htslib-rs`, and `samtools-rs` stabilize,
the default Cargo dependencies can move from local paths to published crate
versions. Keep the submodules for source comparison, fixture access, and local
patching.

## Initial Library Targets

### `bioscript.pysam`

Start with the subset needed for alignment-backed assays:

```python
from bioscript import pysam

with pysam.AlignmentFile(input_file, "rc", reference_filename=reference_file) as bam:
    for read in bam.fetch("22", 36265859, 36266005):
        print(read.query_name, read.reference_start, read.reference_end)
```

Initial surface:

- `AlignmentFile(path, mode="r", reference_filename=None, index_filename=None)`
- `AlignmentFile.fetch(contig, start=None, stop=None)`
- context manager behavior in Python package, equivalent lifecycle in BioScript
- read fields: `query_name`, `reference_name`, `reference_start`,
  `reference_end`, `query_sequence`, `mapping_quality`, `cigarstring`,
  `is_unmapped`, `is_reverse`
- explicit unsupported errors for mutation, writing, remote files, complex tags,
  and full htslib behavior not yet implemented

Backends:

- CRAM and reference FASTA through `noodles` and the existing streaming CRAM
  path.
- BAM can be added after CRAM fetch parity is stable.

Support matrix: [`pysam-support.md`](pysam-support.md).

### `bioscript.pyfaidx`

Start with indexed FASTA lookup:

```python
from bioscript import pyfaidx

fasta = pyfaidx.Fasta(reference_file)
seq = fasta["22"][36265859:36266005]
```

Initial surface:

- `Fasta(path)`
- contig lookup by name
- Python-style slicing
- string conversion for fetched sequence windows

Support matrix: [`pyfaidx-support.md`](pyfaidx-support.md).

### `bioscript.vcf` or `bioscript.pysam.VariantFile`

Prefer `pysam.VariantFile` first if the goal is to minimize import surfaces.
Support:

- open VCF/VCF.GZ
- iterate records
- fetch by region when indexed
- expose `chrom`, `pos`, `id`, `ref`, `alts`, and sample genotype fields

The initial implementation decision is `bioscript.pysam.VariantFile` first for
general VCF compatibility, with `bioscript.vcf` reserved for BioScript-native
helpers. Current helpers include raw `read_kestrel(path)` parsing,
VNtyper-specific `read_vntyper_kestrel(path)` call-table conversion, and
`build_vntyper_report_json(sample_name, input_files, rows)` report generation.
The report helper accepts optional metadata and coverage maps in runtime calls
so BioScript BAM and FASTQ slices can preserve pipeline labels and depth
summaries without embedding report rules in scripts.

## Upstream Source And Tests

Vendored upstream repositories should be kept under `vendor/` as git
submodules when practical:

```text
vendor/python/pysam
vendor/python/pyfaidx
vendor/rust/kestrel-rs
vendor/rust/bcftools-rs
vendor/rust/htslib-rs
```

Reasons to clone upstream code:

- read the real API behavior while implementing shims
- port focused tests for the subset BioScript claims to support
- run selected upstream tests against real libraries where possible
- preserve fixtures and edge cases that are hard to rediscover

Do not run whole upstream suites as a compatibility gate initially. Instead,
copy or adapt targeted tests into BioScript-owned test files, with comments
linking back to upstream test names or files.

## Compatibility Policy

Each shim should document:

- supported constructors, methods, attributes, and argument combinations
- unsupported features with deliberate error messages
- parity tests against real Python libraries when available
- BioScript-specific restrictions caused by sandboxing or deterministic runtime
  requirements

Compatibility should expand by test case. A feature is supported when:

1. It is documented in this file or a module-specific support file.
2. It has Rust tests for `bioscript-libs`.
3. It has runtime tests for BioScript/Monty binding.
4. It has Python package tests when the Python wrapper exists.

Python parity testing is described in [`python-parity.md`](python-parity.md).

## Migration Path For Current Assays

Current assays use:

```python
G1_SITE = bioscript.variant(...)
genotypes = bioscript.load_genotypes(input_file)
site = genotypes.lookup_variant(G1_SITE)
```

Keep that API working while adding shim-based examples. The first migration
target should be an APOL1 proof that computes the same result through:

```python
from bioscript import pysam
```

This lets the project compare current high-level variant lookup behavior against
lower-level alignment-read iteration before replacing any production assay
surface.
