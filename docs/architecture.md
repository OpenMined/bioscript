# BioScript Architecture

This document describes how bioscript is assembled and what each dependency does. It is intended for contributors who need to reason about runtime, genotype sourcing, and CRAM decoding.

## Top-level Layout

```
bioscript/
  rust/                Rust workspace (the executable + core crates)
    bioscript-cli/       binary entry point (`bs`)
    bioscript-runtime/   orchestrates assay execution + timeouts
    bioscript-formats/   genotype sourcing (txt, zip, vcf, cram)
    bioscript-core/      shared types: VariantSpec, GenomicLocus, etc.
    bioscript-schema/    validation of variant/catalogue YAML
    bioscript-ffi/       C-ABI bindings for embedding
  monty/               submodule — Pythonic language runtime
  noodles/             submodule — bioinformatics I/O (our fork)
  bs                   shell shim → `rust/target/…/bioscript`
  assays/              user-authored assay scripts + variant YAMLs
```

Two submodules carry third-party code. Everything else is ours.

## Dependencies

### `monty` — Pythonic runtime

- Upstream: `git@github.com:pydantic/monty.git`
- Vendored as a git submodule at `bioscript/monty`
- Provides the Python-like language that assay authors write (`classify_apol1()`, `count_char()`, etc.)
- Executed via `bioscript-runtime`, which owns the Monty interpreter and enforces per-assay time and memory budgets

### `noodles` — CRAM / BAM / FASTA / CSI / CRAI

- Upstream: `git@github.com:zaeleus/noodles.git`
- **We ship a fork** at `git@github.com:madhavajay/noodles.git`, branch `madhava/streaming-slice-records`, vendored as a submodule at `bioscript/noodles`
- Pulled into the build via `[patch.crates-io]` in `rust/Cargo.toml` — the registry versions of `noodles-cram`, `noodles-core`, `noodles-sam`, etc. are all redirected to the submodule so the crate graph stays unified

The fork carries two upstreamable patches to `noodles-cram`:

1. **`Slice::records_while(...)`** — a public streaming iterator over a CRAM slice's records with an early-termination callback. Upstream only exposes `Slice::records()` which returns `Vec<Record>`, forcing the entire slice (often ~10 000 records) to be decoded even for a single-base locus query. This was the root cause of the original apparent "hang" on chr22 APOL1 lookups.
2. **`validate_reference_md5` flag** — allows callers to opt out of the strict slice-level reference MD5 check when the supplied FASTA disagrees with the CRAM's encoding reference. Matches `samtools`'s default lenient behavior.

Both changes also expose `Slice::header()` and `ReferenceSequenceContext` publicly, which upstream keeps `pub(crate)`.

### `lexical-util` — vendored pin

`rust/vendor/lexical-util` contains a checked-in copy of the `lexical-util` crate (a transitive dep of the Rust number-parsing stack). It is wired in via `[patch.crates-io]` in `rust/Cargo.toml`. Unlike `noodles`, this is *not* a fork with new functionality — it's a pin to a known-good source, isolating us from upstream churn. Treat it as read-only; if it ever needs refreshing, re-copy from a clean crates.io download and re-pin.

## CRAM Read Path (`bioscript-formats/src/alignment.rs`)

Given a single SNP or indel locus, the read path is:

1. Build a `fasta::Repository` from the indexed reference FASTA.
2. Open the CRAM with its CRAI index.
3. **Filter CRAI entries** by reference sequence id + interval overlap. Only containers whose alignment spans touch the locus survive.
4. For each surviving container: `seek` to its byte offset, `read_container(...)`, then iterate its slices and skip any whose landmark is not in the selected set.
5. For each selected slice:
   - `decode_blocks()` — decompresses the CRAM blocks once per slice.
   - `records_while(..., validate_reference_md5 = true, on_record)` — streams records one at a time. For each record we construct an `AlignmentRecord` (start/end/sequence/CIGAR) and either skip it (outside the interval), forward it to the caller, or **stop** (once `record.start > locus.end`, since slices are coordinate-sorted).
6. On `reference sequence checksum mismatch`, the call is retried with `validate_reference_md5 = false`, a loud warning is written to `stderr`, and decoding proceeds. Results may be wrong at positions where the supplied FASTA actually differs from the encoding reference — the warning tells the user to investigate.

Compared to calling upstream `Slice::records()` directly, the streaming path turns decoding ~10 000 records into decoding ~40 — roughly three orders of magnitude less work per locus.

## Performance Expectations

For a single SNP lookup on an aligned whole-genome CRAM:

- `samtools view -T ref.fa file.cram region` — ~40 ms
- bioscript via this path — ~200 ms (three APOL1 loci total)
- previous `Slice::records()` path — >5 minutes, often killed by the assay timeout

The integration tests in `bioscript-formats/tests/file_formats.rs` enforce a 5 s ceiling per single-locus CRAM lookup to catch regressions (e.g. if the streaming path silently breaks and we fall back to full-slice decoding).

## Testing

- `./test-assays.sh` — end-to-end: runs every assay against every input in `test-data/`, reports pass/fail + timings.
- `cargo test -p bioscript-formats --test file_formats` — integration tests including:
  - `cram_apol1_snp_lookup_is_fast_and_correct` — correctness + wall-time budget
  - `cram_md5_mismatch_is_tolerated_and_returns_correct_result` — asserts the fallback path produces a sensible genotype and a read depth within tolerance of `samtools mpileup`

Tests skip gracefully when large test fixtures (the NYGC 1kG CRAM, GRCh38 FASTA) aren't present locally.

## Updating the `noodles` Submodule

When making changes to the noodles fork:

```sh
cd bioscript/noodles
# edit, commit on madhava/streaming-slice-records
git push origin madhava/streaming-slice-records
cd ..
git add noodles
git commit -m "bump noodles submodule"
```

When upstream eventually accepts the streaming + MD5-bypass patches, drop the `[patch.crates-io]` entries in `rust/Cargo.toml`, remove the submodule, and depend on the registry version directly.
