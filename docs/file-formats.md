# File Formats

This document defines the file-format inspection layer for Bioscript: how we detect supported genomics files, what metadata we infer from them, which signals are authoritative vs heuristic, and which fixtures we use to verify that behavior.

The current code already contains pieces of this logic in [`rust/bioscript-formats/src/genotype.rs`](/Users/madhavajay/dev/biovault-app/main/bioscript/rust/bioscript-formats/src/genotype.rs:124) and [`rust/bioscript-formats/src/prepare.rs`](/Users/madhavajay/dev/biovault-app/main/bioscript/rust/bioscript-formats/src/prepare.rs:34). The goal is to make that logic explicit and reusable through a dedicated inspection entrypoint rather than keeping detection embedded inside each loader.

## Goals

- Detect the broad file family from path, container, and sampled contents.
- Support SNP-array genotype exports, VCF/VCF.GZ, CRAM/BAM, and reference FASTA files.
- Return structured metadata that other parts of Bioscript can reuse.
- Separate strong evidence from best-effort guesses.
- Keep detection cheap enough to run before import, in the CLI, in tests, and eventually in UI flows.

## Public API

The inspection layer now lives in [`rust/bioscript-formats/src/inspect.rs`](/Users/madhavajay/dev/biovault-app/main/bioscript/rust/bioscript-formats/src/inspect.rs:1) and currently exposes one main entrypoint:

```rust
pub fn inspect_file(path: &Path, options: &InspectOptions) -> Result<FileInspection, RuntimeError>
```

Current result shape:

```rust
pub struct FileInspection {
    pub path: PathBuf,
    pub container: FileContainer,
    pub detected_kind: DetectedKind,
    pub confidence: DetectionConfidence,
    pub source: Option<SourceMetadata>,
    pub assembly: Option<Assembly>,
    pub phased: Option<bool>,
    pub selected_entry: Option<String>,
    pub has_index: Option<bool>,
    pub index_path: Option<PathBuf>,
    pub reference_matches: Option<bool>,
    pub evidence: Vec<String>,
    pub warnings: Vec<String>,
    pub duration_ms: u128,
}
```

Current source/platform shape:

```rust
pub struct SourceMetadata {
    pub vendor: Option<String>,
    pub platform_version: Option<String>,
    pub confidence: DetectionConfidence,
    pub evidence: Vec<String>,
}
```

Use `platform_version` as a string, not an integer. These labels are vendor-specific and are not guaranteed to form a simple numeric ordering.

The `InspectOptions` surface is intentionally small in the first pass:

```rust
pub struct InspectOptions {
    pub input_index: Option<PathBuf>,
    pub reference_file: Option<PathBuf>,
    pub reference_index: Option<PathBuf>,
}
```

The CLI exposes the same functionality through:

```text
bs inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>]
```

The current implementation prints a stable tab-separated report from `FileInspection::render_text()`, which is suitable for shell use and for downstream callers that just need a lightweight CLI probe.

Important rule: prefer `None` over a false claim.

Examples:
- `assembly: None` means "not enough evidence", not "not GRCh38".
- `phased: None` means "not applicable or not checked", not "unphased".
- `reference_matches: None` means "not validated", not "validated successfully".

## Supported Families

### Genotype text exports

These are consumer SNP-array exports such as:
- 23andMe
- AncestryDNA
- deCODEme
- Dynamic DNA / GSAv3
- FamilyTreeDNA
- Genes for Good
- Living DNA
- MyHeritage

Observed format families:
- Standard 4-column layout: `rsid chromosome position genotype`
- Split-allele layout: `rsid chromosome position allele1 allele2`
- Extended metric layout: genotype columns plus metrics such as `gs`, `baf`, `lrr`

The detector should not require fixed column numbers. It should scan candidate columns for:
- SNP id: `rs...` or `i...`
- chromosome
- integer position
- genotype or split alleles

This is the main robustness improvement needed over the current fixed/default fallback.

### VCF

Supported:
- `.vcf`
- `.vcf.gz`
- `.zip` containing `.vcf` or `.vcf.gz`

Metadata we want:
- sample count
- phased vs unphased
- contig naming style (`chr1` vs `1`)
- assembly if clearly encoded in header or reference metadata

### Raw alignments

Supported:
- `.cram`
- `.bam`

Metadata we want:
- alignment subtype (`cram` / `bam`)
- index presence and location (`.crai` / `.bai`)
- whether a supplied reference is required
- whether a supplied reference can be used successfully
- whether the reference appears to mismatch the CRAM's encoding reference

### Reference FASTA

Supported:
- `.fa`
- `.fasta`
- compressed support can be added later if we need it

Metadata we want:
- FASTA subtype
- `.fai` presence
- contig naming style
- assembly hint when file/path clearly encode it

## Detection Pipeline

Detection should run in layers from cheap to expensive.

### 1. Path-based hints

Use filename and extension to form an initial guess:
- `.zip` => archive container
- `.vcf`, `.vcf.gz` => VCF candidate
- `.cram`, `.bam` => alignment candidate
- `.fa`, `.fasta` => reference candidate
- filename or header tokens like `GRCh38`, `hg19`, `GSAv3`, `v5`, `23andMe`, and similar assembly/source labels are hints, not proof

### 2. Container inspection

If the file is a zip:
- list entries
- ignore directories and `__MACOSX/`
- prefer entries ending in `.vcf`, `.vcf.gz`, `.txt`, `.tsv`, `.csv`
- keep the selected entry name in the inspection result

This is already partially implemented in `select_zip_entry(...)` in [`genotype.rs`](/Users/madhavajay/dev/biovault-app/main/bioscript/rust/bioscript-formats/src/genotype.rs:1189).

### 3. Content sniffing

Sample a bounded number of lines and extract signals.

For genotype text files:
- delimiter detection: tab, comma, whitespace
- header alias mapping
- comment preamble scanning
- flexible row validation: rsid/id, chromosome, position, genotype
- classify as genotype only if a strong fraction of sampled rows look valid

For VCF:
- `##fileformat=VCF`
- `#CHROM`
- parse `GT` fields just far enough to determine sample count and phasing
- scan header metadata for source/vendor fingerprints when present, e.g. `sequencing.com`

For alignments:
- subtype comes from extension
- index presence comes from adjacent path checks
- reference compatibility requires an explicit validation step

### 4. Metadata inference

Once the file family is known, infer metadata with explicit precedence.

Recommended precedence:

1. File-internal authoritative metadata
2. File-internal human-readable header text
3. Container entry name
4. Outer filename
5. Source-specific fallback rules

## Heuristics By Family

### Genotype text

Signals to infer:
- delimiter
- column mapping
- vendor/source
- platform version
- assembly

Source/vendor detection should look at header text and column shapes. Useful fingerprints already observed in local datasets:
- 23andMe: header mentions `23andMe`
- Dynamic DNA: header mentions `Dynamic DNA`, `DDNA`, or `dynamicdnalabs`
- Extended GSAv3-style files expose `gs`, `baf`, `lrr`

Platform-version detection should stay conservative:
- filename tokens like `v2`, `v3`, `v4`, `v5`, `GSAv3`, `GSAv3-DTC`
- vendor-specific export names
- presence or absence of a vendor-specific marker panel, using a curated set of informative rsids as evidence
- maybe later: row-count ranges as a weak signal, but not as primary evidence

For 23andMe specifically, chip version should not be inferred from filename alone. A better approach is:
- use filename/header/export-name hints as weak evidence
- compare the file's marker set against curated version-specific rsid panels
- return the `platform_version` string with supporting evidence rather than a bare guess

For other vendors, use the same model:
- `vendor = "23andMe"`, `platform_version = Some("v5")`
- `vendor = "Dynamic DNA"`, `platform_version = Some("GSAv3-DTC")`
- `vendor = "AncestryDNA"`, `platform_version = None`

Assembly detection should use header text first:
- `build 36`, `GRCh36`, `hg18`
- `build 37`, `GRCh37`, `hg19`
- `build 38`, `GRCh38`, `hg38`

Known examples in the shared data:
- 23andMe v2/v3 exports declare build 36
- 23andMe v4/v5 exports declare build 37
- Dynamic DNA GSAv3 export declares build 38

### VCF

Signals to infer:
- compressed or plain
- sample count
- phased vs unphased
- contig naming style
- assembly, when the header contains enough information

Phasing rule:
- any `GT` using `|` => `Some(true)`
- all parsed `GT` use `/` => `Some(false)`
- no usable `GT` seen => `None`

### CRAM / BAM

Signals to infer:
- subtype
- adjacent index presence
- reference requirement
- reference compatibility status

Reference compatibility states should be explicit, for example:
- `NotChecked`
- `MissingReference`
- `Compatible`
- `Md5MismatchButReadable`
- `Incompatible`

For CRAM this should be based on an actual open/decode attempt against the supplied reference, not just on filename matching. The existing alignment path already has a real MD5-mismatch fallback behavior documented in [`docs/architecture.md`](/Users/madhavajay/dev/biovault-app/main/bioscript/docs/architecture.md:54).

### Reference FASTA

Signals to infer:
- subtype
- adjacent `.fai`
- contig naming style
- assembly hint from filename/path

## Confidence Model

The inspector should carry confidence so callers can decide whether to trust or surface a result.

Suggested levels:
- `Authoritative`
- `StrongHeuristic`
- `WeakHeuristic`
- `Unknown`

Examples:
- VCF identified from `##fileformat=VCF` => `Authoritative`
- 23andMe detected from header text => `StrongHeuristic`
- 23andMe `platform_version` guessed only from a filename token like `v5` => `WeakHeuristic`
- 23andMe `platform_version` supported by a version-specific rsid panel match => `StrongHeuristic`

## Fixture Matrix

The test suite should cover both synthetic and real-world fixtures.

### Fixtures already in this repo

`rust/bioscript-formats/tests/fixtures`
- `mini.cram`, `mini.cram.crai`
- `mini.fa`, `mini.fa.fai`
- `mini_bad_ref.fa`, `mini_bad_ref.fa.fai`

Use these for:
- CRAM subtype detection
- CRAI/FAI detection
- reference compatibility and MD5-mismatch behavior

`old/examples/apol1/genotype_files`
- Dynamic DNA GSAv3-style plain-text files
- filenames include `GSAv3` and `GRCh38`
- header includes build 38 wording

Use these for:
- extended genotype layout
- build 38 inference
- chip/version filename heuristics

### Shared test data

Repo-local shared dataset now lives under [`test-data/`](/Users/madhavajay/dev/biovault-app/main/bioscript/test-data:1), populated by [`tools/fetch_test_data.sh`](/Users/madhavajay/dev/biovault-app/main/bioscript/tools/fetch_test_data.sh:1) as symlinks into the shared cache at `~/.bioscript/cache/test-data/`:

- `23andme/v2/.../23data20100526.txt.zip`
- `23andme/v3/.../huE4DAE4_20120522224129.txt.zip`
- `23andme/v4/.../genome__v4_Full_2016.txt.zip`
- `23andme/v5/.../genome_hu50B3F5_v5_Full.zip`
- `dynamicdna/.../100001_X_X_GSAv3-DTC_GRCh38-07-12-2025.txt.zip`
- `1k-genomes/ref/...fa`, `.fai`
- `1k-genomes/aligned/...cram`, `.crai`

These are the right fixtures for real-world detection coverage:
- 23andMe vendor + version + build heuristics
- Dynamic DNA vendor + GSAv3 + GRCh38 heuristics
- large real CRAM/reference/index presence and compatibility checks

The fetch tool and manifest now live in this repo:
- [`tools/fetch_test_data.sh`](/Users/madhavajay/dev/biovault-app/main/bioscript/tools/fetch_test_data.sh:1)
- [`tools/sources.yaml`](/Users/madhavajay/dev/biovault-app/main/bioscript/tools/sources.yaml:1)

That makes the same cache-backed layout available in local development, CI, and downstream repos.

## Tests To Add

### Current coverage

- generated `.vcf` and `.zip` containing VCF are covered in `bioscript-formats`
- real-world zipped 23andMe fixtures are readable through `GenotypeStore`
- real-world zipped Dynamic DNA fixture is readable through `GenotypeStore`
- `inspect_file(...)` metadata assertions now cover:
  - AncestryDNA
  - FamilyTreeDNA
  - Genes for Good
  - Dynamic DNA
  - CRAM index detection
  - phased VCF detection
- the `bioscript inspect` CLI has a smoke test in `bioscript-cli`
- existing mini CRAM fixtures continue to validate reference/index behavior

### Next

- version-specific rsid-panel heuristics for 23andMe
- source detection tests for more vendors from curated real files
- build inference tests for 36 / 37 / 38 across more sources
- missing-index and mismatched-reference tests for alignments
- JSON output mode for `bs inspect`

### Later

- AncestryDNA
- FamilyTreeDNA
- Living DNA
- MyHeritage
- deCODEme
- ambiguous or malformed files
- BAM fixtures

## Shared Test-Data Layout

To avoid duplicate downloads across repos, the shared fetch tool should populate a cache under:

```text
~/.bioscript/cache/test-data/
```

Each consuming repo can then expose a local `test-data/` tree as symlinks into that cache by running its local `tools/fetch_test_data.sh`.

This gives us:
- one physical copy of large CRAM/FASTA fixtures
- stable local paths for tests
- the option to reuse the same data across Bioscript, ExVitae, and Biovault

## Non-Goals

The inspection layer should not:
- fully parse or normalize every file up front
- silently auto-fix mismatches
- claim certainty where only weak hints exist

The file inspection API should not mutate inspected input files, but the fetch tooling is explicitly allowed to materialize a local `test-data/` tree as symlinks into the shared cache.

Its job is to describe the file accurately enough that the caller can decide what to do next.
