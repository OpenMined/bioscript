# Fix Kestrel VNtyper FASTQ Parity

## Goal

Fix the `kestrel-rs` VNtyper FASTQ parity blocker so the Rust implementation
matches the Java Kestrel outputs closely enough for BioScript VNtyper FASTQ
classification, TSV fingerprint, and report JSON parity.

## Current Evidence

The BioScript opt-in parity gate currently fails:

```sh
BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1 \
PYTHONPATH=python:ports/vntyper/bioscript \
python -m unittest ports.vntyper.tests.test_native_fastq_pipeline_gate.VntyperNativeFastqPipelineGateTests.test_native_fastq_pipeline_with_native_kestrel_and_bcftools_matches_expected_classification
```

Observed failure:

- Positive FASTQ case: Rust emits `2417` rows, Java expected output has `3737`.
- Negative FASTQ case: Rust emits `2322` rows, Java expected output has `4897`.
- Negative classification is wrong: Rust reports `High_Precision`, Java
  expected classification is `negative`.

The existing vendor-level gate is:

```sh
cd /home/linux/dev/bioscript/workspace1/vendor/rust/kestrel-rs
KESTREL_RUN_VNTYPER_FASTQ_PARITY=1 CC=cc AR=ar \
  cargo test -p kestrel --test vntyper_fastq_parity -- --nocapture
```

## Work Plan

1. Create a dedicated branch in
   `/home/linux/dev/bioscript/workspace1/vendor/rust/kestrel-rs`.
2. Run the normal Kestrel Rust test baseline before editing code.
3. Run the existing opt-in VNtyper FASTQ parity gate and save the failing
   evidence.
4. Add or tighten a focused test that reproduces the VNtyper false-positive /
   row-count mismatch at the smallest practical scope.
5. Compare Java Kestrel artifacts and Rust Kestrel artifacts for the same FASTQ
   inputs:
   - retained Rust VCF,
   - Java expected VCF,
   - shared/missing/extra record keys,
   - GDP and DP differences,
   - passing VNtyper-filter rows.
6. Fix the Rust Kestrel implementation in the vendor repo. Prefer matching Java
   Kestrel semantics over adding BioScript-side filters.
7. Verify:
   - normal Kestrel Rust tests pass,
   - new focused regression test passes,
   - vendor VNtyper FASTQ parity gate passes or has only explicitly accepted
     non-behavioral differences,
   - BioScript native FASTQ parity gate passes from the workspace root.

## Non-Goals

- Do not hide the parity gap in BioScript post-processing.
- Do not change VNtyper expected fixtures unless Java Kestrel evidence proves
  the fixture is wrong.
- Do not weaken parity assertions just to make the gate green.

## Status

- Branch created in `vendor/rust/kestrel-rs`: `fix/vntyper-fastq-parity`.
- Baseline before edits:
  - `CC=cc AR=ar cargo test --workspace` passed.
  - Opt-in VNtyper FASTQ parity failed:
    - Negative: Rust `2322` records vs Java expected `4897`.
    - Positive: Rust `2417` records vs Java expected `3737`.
- Added a focused Rust regression in `crates/kestrel/src/runner.rs`:
  `graph_haplotypes_assembles_overlapping_kmer_path_without_full_read`.
  This covers a k-mer graph path that is not backed by one full read sequence.
- Replaced the temporary read-backed haplotype path with a Kestrel-style
  k-mer branch traversal using `KmerAligner` saved states and
  `HaplotypeContainer`.
- Added a bounded repeat/sequence guard so repeated k-mer branches cannot
  restore forever.
- Added a focused active-region regression in
  `crates/kestrel/src/activeregion/mod.rs`:
  `active_region_detector_splits_repetitive_peaks_at_last_stable_valley`.
  This covers the VNtyper `N-S`-like repetitive profile that the original Rust
  port missed. Before the fix, Rust merged the two Java regions into one large
  active region; after the fix, it splits them at `(4, 43)` and `(60, 94)`.
- Ported Java's right-scan peak/valley fallback into Rust active-region
  detection.
- Matched Java's saved-state capacity tie behavior in `KmerAligner`: when
  equal minimum-depth saved states compete for removal, Java's linked stack
  removes the newest equal-depth state, not the oldest. Added
  `kmer_aligner_capacity_removes_newest_equal_min_depth_like_java_stack`.
- Added runner-side deduplication for cloned saved-state haplotypes and saved
  branch states. Java saved states share `MaxAlignmentScoreNode` objects and
  suppress already-built haplotypes through shared `haplotypeBuilt` flags; Rust
  deep-clones those nodes, so duplicate haplotypes/states need explicit
  suppression.
- Refactored `KmerAligner` trace nodes to shared `Rc<TraceNode>` references so
  saved alignment states keep Java-like shared traceback structure instead of
  deep-cloning large trace graphs on every state save.
- Added a reduced N-S insertion regression in `crates/kestrel/src/runner.rs`:
  `graph_haplotypes_recovers_reduced_vntyper_ns_insertion_branch`. This reduced
  static-count graph recovers Java's `sample-N-S-61-72` insertion sequence:
  `TGGGGGGGCGGTGGAGCCCGGGGCCGGGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACC`.
- Rechecked Java `KmerHashSet.HashElement`: it copies k-mer arrays when adding
  elements. The earlier mutable bucket-history hypothesis was wrong. Rust's
  exact `HashSet<KmerKey>` repeat detection is the correct model for this path,
  and keeping exact detection is what lets the reduced insertion regression pass.
- Vendor work is committed in `vendor/rust/kestrel-rs`:
  `63bbbe4 Fix Kestrel VNtyper graph traversal parity`.

## Verification So Far

These pass after the Rust runner change:

```sh
cd /home/linux/dev/bioscript/workspace1/vendor/rust/kestrel-rs
CC=cc AR=ar cargo test -p kestrel runner::tests:: -- --nocapture
CC=cc AR=ar cargo test -p kestrel align::tests:: -- --nocapture
CC=cc AR=ar cargo test -p kestrel active_region_detector_splits_repetitive_peaks_at_last_stable_valley -- --nocapture
CC=cc AR=ar cargo test -p kestrel kmer_aligner_capacity_removes_newest_equal_min_depth_like_java_stack -- --nocapture
CC=cc AR=ar cargo test -p kestrel add_unique_haplotype_skips_duplicate_sequence_and_alignment -- --nocapture
CC=cc AR=ar cargo test -p kestrel graph_haplotypes_recovers_reduced_vntyper_ns_insertion_branch -- --nocapture
CC=cc AR=ar cargo test -p kanalyze hash_is_deterministic -- --nocapture
CC=cc AR=ar cargo test -p kanalyze inserts_contains_removes_and_clones_independently -- --nocapture
CC=cc AR=ar cargo test --workspace
```

The focused tests above were re-run after the reduced insertion fix and pass.
`CC=cc AR=ar cargo test --workspace` passed earlier in this branch after the
active-region and saved-state changes; it was not re-run after the latest
reduced insertion regression.

The opt-in VNtyper FASTQ parity gate still fails after the reduced insertion
fix:

```sh
rm -rf /tmp/kestrel-vntyper-parity-current
KESTREL_RUN_VNTYPER_FASTQ_PARITY=1 \
KESTREL_VNTYPER_PARITY_OUT=/tmp/kestrel-vntyper-parity-current \
CC=cc AR=ar cargo test -p kestrel --test vntyper_fastq_parity -- --nocapture
```

Current failed counts:

- Positive: Rust `1804` records vs Java expected `3737`.
  Shared `1770`, missing `1967`, extra `34`.
- Negative: Rust `2217` records vs Java expected `4897`.
  Shared `2135`, missing `2762`, extra `82`.

The failed artifacts are retained under
`/tmp/kestrel-vntyper-parity-current`.

The retained positive FASTQ artifacts still show Rust missing Java's
`N-S:86 G>GGGTGGAGCCCGGGGCCGG` VCF record under the parity harness's bounded
`max_haplotypes=2` / `max_aligner_states=2` configuration, even though the
reduced static N-S regression emits the insertion under Java default-like
`10/15` traversal caps.

A positive FASTQ probe with `KESTREL_VNTYPER_MAX_ALIGNER_STATES=10` and
`KESTREL_VNTYPER_MAX_HAPLOTYPES=15` was started to test Java-effective caps but
was interrupted after running beyond a minute. The lingering cargo/test
processes were stopped before committing.

## Why This Was Missed

The original Rust unit tests did not include a repetitive VNtyper-like
active-region count profile. They covered simpler count drops/recoveries and
runner graph assembly, but not Java's repeated peak/valley fallback in
`ActiveRegionDetector.scanRight`. That allowed the Rust port to pass unit tests
while incorrectly merging Java's two `N-S` active regions into one large region.

The new active-region regression reproduces that missing Java behavior directly
from a reduced `N-S` profile.

## Current Blocker

The active-region split is fixed, and a reduced static N-S graph now emits the
Java insertion haplotype. Full VNtyper FASTQ parity remains blocked because the
parity harness still runs Rust with `max_haplotypes=2` and
`max_aligner_states=2`, while Java's CLI path appears to reset those caps to
builder defaults after `setMaxRepeatCount(0)` reconstructs
`KmerAlignmentBuilder`.

Current observed behavior:

- Reduced static regression at `10/15`: Rust emits the expected insertion
  branch and the test passes.
- Full FASTQ parity at `2/2`: Rust still misses low-depth Java records,
  including `N-S:86 G>GGGTGGAGCCCGGGGCCGG`.
- Full FASTQ parity at `10/15`: not confirmed. A positive-case probe ran longer
  than a minute and was interrupted.

The remaining work is therefore not in BioScript post-processing, BCFtools, or
Samtools. It is in Kestrel Rust's Java-cap parity/performance behavior:

- Decide whether Rust should intentionally emulate Java runner ordering, where
  `setMaxRepeatCount` rebuilds the alignment builder after aligner/haplotype
  caps are set.
- If Java-effective defaults are required, fix saved-state traversal
  performance enough for `10/15` FASTQ parity to complete.
- Keep the reduced N-S insertion test as the fast inner loop before repeating
  broad FASTQ probes.

## Current Thinking

The earlier bucket-history repeat hypothesis was disproved. Java
`KmerHashSet.HashElement` copies k-mer arrays on insertion, so Java repeat
detection is exact k-mer membership, not mutable bucket history. Rust should
keep exact `KmerHashSet::insert(kmer.clone())` cycle detection.

The most useful reduced target is now covered by a passing test:

```text
sample-N-S-61-72
CIGAR: 20=1X6=18I4=1X1=1X20
VCF:   N-S:86 G>GGGTGGAGCCCGGGGCCGG
```

Avoid these paths unless a smaller unit test justifies them:

- Do not widen `region_sequence_limit(...)` again; full-reference, `2*k`, and
  `1.5*k` guards were tried and reverted.
- Do not just disable repeat detection globally; it causes unacceptable
  traversal growth.
- Do not reintroduce bucket-based repeat detection; Java does not work that way.
- Do not rerun broad FASTQ parity loops as the primary debug loop unless the
  next change is specifically about cap parity or state traversal performance.

## Completion Audit

Objective from the original request:

1. Create a branch in `/home/linux/dev/bioscript/workspace1/vendor/rust/kestrel-rs`.
2. Write this markdown first in the workspace root.
3. Confirm normal tests pass before behavior changes.
4. Add a test for the VNtyper/Kestrel problem.
5. Fix the Rust code.
6. Verify the result against both original Java Kestrel and new Rust Kestrel.

Current evidence:

- Branch: done. Current Kestrel branch is `fix/vntyper-fastq-parity`.
- Markdown: done. This file is
  `/home/linux/dev/bioscript/workspace1/fix-kestrel.md`.
- Baseline tests: done before edits. `CC=cc AR=ar cargo test --workspace`
  passed before behavior changes.
- Current focused tests: done. Runner, aligner, active-region, and kanalyze
  focused tests pass after the reduced insertion fix.
- Reduced regressions: done for two confirmed misses. The test
  `active_region_detector_splits_repetitive_peaks_at_last_stable_valley`
  reproduces the Java right-scan peak/valley fallback that the original Rust
  port lacked. The test
  `graph_haplotypes_recovers_reduced_vntyper_ns_insertion_branch` now
  reproduces and recovers the Java N-S insertion branch.
- Rust fix: partially done. The active-region split for the reduced `N-S`
  profile is fixed, saved-state equal-depth pruning now matches Java's
  linked-stack tie behavior, duplicate saved-state haplotypes/branches are
  suppressed in the Rust runner, and saved alignment states share traceback
  nodes instead of deep-cloning them.
- Java/Rust verification: partially done. The reduced `N-S` active regions now
  match Java and the reduced static insertion branch is recovered, but full
  VNtyper FASTQ parity still fails. Java still emits haplotypes and VCF records
  that Rust does not under the Rust harness's bounded `2/2` caps.

Not complete:

- The opt-in vendor parity gate still fails:

  ```sh
  KESTREL_RUN_VNTYPER_FASTQ_PARITY=1 \
  KESTREL_VNTYPER_PARITY_OUT=/tmp/kestrel-vntyper-parity-peak \
  CC=cc AR=ar cargo test -p kestrel --test vntyper_fastq_parity -- --nocapture
  ```

- The BioScript native FASTQ gate remains blocked until Kestrel Rust matches
  Java's haplotype/VCF output:

  ```sh
  BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1 \
  PYTHONPATH=python:ports/vntyper/bioscript \
  python -m unittest ports.vntyper.tests.test_native_fastq_pipeline_gate.VntyperNativeFastqPipelineGateTests.test_native_fastq_pipeline_with_native_kestrel_and_bcftools_matches_expected_classification
  ```

The current blocker is therefore not in BioScript or BCFtools/Samtools. It is
inside `kestrel-rs` haplotype graph traversal, saved-state pruning, or aligner
continuation/performance behavior.
