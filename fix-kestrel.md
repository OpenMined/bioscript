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

## 2026-05-15 Update: Empirical confirmation of cap-reset and traversal divergence

### Java cap-reset bug confirmed empirically
Running the Java jar (`vendor/rust/kestrel-rs/kestrel/lib/kestrel.jar`) against
the negative FASTQ with three different cap settings produces **byte-identical**
output:

```sh
java -jar kestrel.jar -k 20 --maxalignstates 2  --maxhapstates 2  ...   # md5 cb0ed3...
java -jar kestrel.jar -k 20 --maxalignstates 10 --maxhapstates 15 ...   # md5 cb0ed3...
java -jar kestrel.jar -k 20 --maxalignstates 40 --maxhapstates 40 ...   # md5 cb0ed3...
```

All three produce the same 4897 records that the expected fixture contains
(after sorting). This proves Java's CLI silently runs at `DEFAULT_MAX_STATE=10`
/ `DEFAULT_MAX_HAPLOTYPES=15` regardless of the flags, because
`ActiveRegionDetector.setMaxRepeatCount(int)` calls `initAlignmentBuilder()`
which constructs a fresh `KmerAlignmentBuilder` with default caps, throwing
away the user-supplied `setMaxAlignerState` / `setMaxHaplotypes`.

So the parity test's `2/2` defaults are wrong — Java's expected output was
generated at effective `10/15`.

### Even at matching 10/15 caps, Rust still emits ~70 % more records
Running Rust at `KESTREL_VNTYPER_MAX_ALIGNER_STATES=10` /
`KESTREL_VNTYPER_MAX_HAPLOTYPES=15` (negative case, release mode, ~8.5 min):

- Rust: 8269 records vs Java 4897 (shared 4272, missing 625, extra 3997).
- Rust per-record type distribution skews heavily toward insertions:
  `del:112, ins:2589, snp:5568` vs Java `del:75, ins:390, snp:4432`.
- Rust GDP bucket distribution has tons of low-GDP records (`1:679, 2-5:620,
  6-20:747, 21-100:529, >100:5694`) while Java has almost everything in
  `>100:4878` and only `2:21-100, 8:6-20, 9:21-100, 2:2-5` outside.
- Active-region detection counts match almost exactly (Rust 980 vs Java 976).

So the parity gap is in **haplotype graph traversal**, not active-region
detection.

### Per-region haplotype-count distributions diverge sharply
After instrumenting `[KDBG-BUILD]` in `build_forward_haplotypes` /
`build_reverse_haplotypes`:

- Java max haplotypes per region = **8**. Distribution peaks at 4 (237 regions)
  and 7 (201 regions).
- Rust max haplotypes per region = **15** (the cap). **501 of ~993 regions hit
  the cap**, generating thousands of unique haplotype keys per region.

For the worst Rust region `J-R:4-119`:
- Rust: 219,920 outer iters, 4040 raw emits, 3771 unique emitted, 15 in
  container. Save attempts 1,689,188 / accepts 302,576 / rejects 1,386,612
  (18 % accept rate).
- Java (same region): 446 save attempts, 408 rejects, 38 successful saves, 28
  evictions, **0 haplotypes emitted** ("Built 0 haplotypes (fwd)").

So Java's traversal never produces any trace that reaches `refLength - 1` with
positive score for this region (the chain stays empty even after 38 restore
cycles). Rust's traversal reaches end-of-region thousands of times.

### Findings on what is NOT the cause
- Toggling the runner-side `saved_states: HashSet<SavedBranchKey>` dedup off
  (`KESTREL_DISABLE_STATE_DEDUP=1`) does not change the result — keys never
  collide, so the dedup is a no-op for this workload.
- Toggling `region_sequence_limit` off (`KESTREL_DISABLE_SEQ_LIMIT=1`) makes
  the divergence **worse** (higher iter counts).
- `Base::ALL` ordering matches Java's A,C,G,T order.
- `state_min_depth`, save-rejection logic, `remove_min_state` tie behaviour,
  `add_base` return semantics, and `record_max_node` all match Java
  line-for-line.
- `KmerHashSet::insert` (Rust) and Java `KmerHashSet.add(int[])` both copy
  k-mers on insertion (no mutable-bucket-history difference).
- `extend_kmer` / `kUtil.append` produce byte-identical encoded k-mers.

### Active region retry: the missing piece

A direct comparison of active-region traces in `J-R` finally exposed the
biggest divergence: **Java retries overlapping active regions from
`refCountIndex + 1` whenever haplotype assembly returns zero (or wildtype-only)
haplotypes**. Rust's pipeline does not. Java's `KestrelRunner.exec` walks
`refCountIndex` one base at a time when haps fail; Rust's
`detect_active_regions` returns a static list and the runner consumes each
region exactly once.

For the `J-R` reference Java tries five overlapping active regions —
`4-119`, `11-119`, `18-119`, `19-60`, `41-119` — and rejects the first four
because their wider spans hit cycles before reaching the right anchor. Only
`J-R:41-119` succeeds and produces the 8 haplotypes that yield the 9 expected
VCF records. Rust's detector emits only `J-R:4-119`, accepts it (since
Rust's traversal happens to reach the right end), produces 15 noisy
haplotypes whose minimum k-mer depths are low, and emits a different mix of
VCF records.

So the missing fix is at the detector–runner interface, not (only) inside the
haplotype graph:

1. Replicate Java's `KestrelRunner.exec` flow: each iteration of the main
   `REF_SEARCH` loop tries one candidate region. Build haplotypes for it
   immediately. If the result is empty or wildtype-only, advance
   `refCountIndex` by 1; otherwise skip past the region. This must be done
   for both right-anchor and left-anchor scans.
2. Implement Java's `setMaxRepeatCount`-driven cap reset (already added as
   `apply_java_cli_cap_reset` in `run_pipeline`).
3. Keep the haplotype trim, capacity, and dedup logic as-is.

The second-order question — why Rust's `J-R:4-119` produces 15 haplotypes
where Java's produces 0 — likely resolves on its own once Java-style
overlap-retry is in place, because Java's narrower retry region
`J-R:41-119` is exactly the region whose haplotypes match the expected VCF.
If Rust starts emitting from `J-R:41-119`, the wider `J-R:4-119` is no
longer the only candidate and the noisy haplotype set should match Java
without any change to graph traversal.

### Current Blocker

The active-region split is fixed, and a reduced static N-S graph now emits the
Java insertion haplotype. Full VNtyper FASTQ parity remains blocked because:

1. The parity harness used the wrong caps (Java's `2/2` is silently `10/15`).
   **Fixed** with `apply_java_cli_cap_reset` in `run_pipeline`.
2. Rust never applied Java's default `kmercount:5` post-count filter. Java's
   `KestrelRunnerBase.getCountModule()` adds the filter whenever
   `minKmerCount > 0`; Rust kept the field on the config but never applied
   it. **Fixed** with `MemoryCountMap::with_min_count` /
   `IkcCountMap::with_min_count` + `KmerCounter::retain`. Also updated the
   parity test to use `min_kmer_count=5` (Java's effective default) instead
   of `1`.
3. Active-region detector didn't retry overlapping regions when haplotype
   assembly produced 0 / wildtype-only haplotypes. **Fixed** with
   `ActiveRegionDetector::detect_from_counts_with`, a callback-driven
   variant that mirrors Java's `REF_SEARCH` loop.

After these three fixes the negative VNtyper FASTQ case now produces 7062
records vs Java 4897 (shared 4335, missing 562, extra 2727). That is a 33%
reduction in extras from the pre-fix state of 4040 extras. The test now
completes in ~93s instead of ~520s. K-mer counts and per-step choose_branch
decisions now match Java's trace line-for-line for the J-R:4-119 region.

### Remaining gap (in progress)

Even with the kmercount filter Rust still emits more haplotypes per region
than Java for wide repetitive regions. Example: 4-5:3-88 — Java assembles 0
haplotypes and retries with narrower 4-5:48-88; Rust assembles 6 haplotypes
from 4-5:3-88 and never reaches 4-5:48-88. Save attempts/accepts:

- Java 4-5:3-88: 503 attempts, 466 rejects, 37 accepts (93% reject), 0 haps.
- Rust 4-5:3-88: 12,745 attempts, 4,454 rejects, 8,291 accepts (35% reject),
  6 haps.

So Rust's saved-state acceptance rate is still much higher than Java's
despite matching k-mer counts and matching choose_branch decisions on the
first ~20 inner iterations. The candidates for the remaining work:

- Investigate whether Rust's saved alignment matrices accumulate scores in
  a way that lets a later restored state propagate higher scores than
  Java's, allowing more chain entries to record max alignments.
- Check whether Rust's haplotype container or `MaxAlignmentScoreNode` chain
  retains nodes that Java naturally drops via shared-mutable
  `haplotypeBuilt` flag semantics.
- Verify whether Java's CountModule has an additional filter (e.g. read
  length minimum, segment cutoff) that is being applied to FASTQ input
  before counting.

### Quantitative progress summary

| Step | Negative parity result | Test time |
|------|------------------------|-----------|
| Initial state (2/2 caps, no kmercount) | 2322 vs 4897, missing 2762, extra 82 | ~10 min |
| 10/15 caps (no kmercount, no overlap retry) | 8269 vs 4897, missing 625, extra 3997 | ~8 min |
| 10/15 caps + overlap retry | 8376 vs 4897, missing 561, extra 4040 | ~8 min |
| 10/15 caps + overlap retry + kmercount:5 | 7062 vs 4897, missing 562, extra 2727 | ~93 s |
| 10/5 caps + overlap retry + kmercount:5 (manual test) | 4563 vs 4897, missing 1371, extra 1037 | ~93 s |

The kmercount filter alone closed ~33 % of the gap and cut test time by ~5×.
Forcing `max_haplotypes=5` closes the gap further but undershoots Java's
record count — that knob is therefore not the right fix on its own. The
remaining work is in the haplotype graph traversal itself: Rust's accept
rate during state save (~35–75 %) needs to converge to Java's ~90 %, and
Rust's `MaxAlignmentScoreNode` chain emissions per region need to drop
from ~1750 to Java's ~5–8.

### What is verified clean

- `apply_java_cli_cap_reset` (replicates Java's CLI cap-reset bug). Empirical
  proof: Java jar at `--maxalignstates 2,10,40` produces byte-identical
  output md5 `cb0ed3...`, matching the expected fixture sorted.
- `KmerCounter::retain` + `MemoryCountMap::with_min_count` /
  `IkcCountMap::with_min_count` (replicates Java's kmercount:5 default).
  Verified: for k-mer `GGCGGTGGAGCCCGGGGCCA` in the negative FASTQ, manual
  occurrence count is 6 (1 fwd + 5 revComp); kanalyze CLI without
  `-rduplicate` returns 1 fwd + 5 revComp = 6; Java in-runtime sums to 5
  because the forward occurrence (count=1) is dropped by `kmercount:5`,
  giving 0 + 5 = 5; Rust now matches Java when `min_kmer_count=5`.
- `ActiveRegionDetector::detect_from_counts_with` callback API (replicates
  Java `REF_SEARCH` overlap retry). Verified by inspection of Java trace.
- Per-step choose_branch decision parity for the first 20+ inner iterations
  of the J-R:4-119 region. Verified via `KESTREL_TRACE_REGION` trace
  comparison to Java's `Saving state` log lines.

### Best lead for the next session

Inner-loop break-cause counters added under `KESTREL_DEBUG_BUILD`. For
J-R:4-119 in Rust (26,894 outer iters):

- `cycle_breaks=1256` (4.7 %)
- `choose_none_breaks=3601` (13.4 %)
- `add_base_false_breaks=17,871` (66.5 %)
- `seq_limit_breaks=4166` (15.5 %)

Java for the same region: 11 cycle breaks out of 38 outer iters = 29 %.
So Rust's cycle-break rate per inner iter is ~6× lower than Java's. The
dominant Rust exit path is `addBase returns false`, which fires when
`max_pot_score < max_alignment_score`. With Rust's chain growing to ~1753
unique entries vs Java's ~9, Rust's `max_alignment_score` likely rises
faster than Java's during a traversal, causing addBase to return false
earlier and the outer loop to restart more often. Each restart begins
from a saved state, generating more saves and continuing the explosion.

So the remaining work is to understand why Rust's
`MaxAlignmentScoreNode` chain accumulates more entries than Java's per
unit of traversal. Candidates:

- Rust's `record_max_node` fires for both align- and gap_con-matrix
  end-of-row positives. Verify Java emits at exactly the same conditions.
- Java's `MaxAlignmentScoreNode` linked list is mutated in place via
  shared `haplotypeBuilt` flags; Rust deep-clones on save_state. The
  runner-side `emitted` dedup catches duplicates at emission time but
  does not prune the chain itself, so a long chain may persist across
  many restore_state cycles and contribute to chain-driven `addBase`
  early-exits.
- A focused unit test that constructs a static count map for J-R:4-119
  and steps add_base / save_state / restore_state until the chain hits
  the expected refLength position would isolate this. The data inputs
  needed for that test are: J-R reference (already in /tmp/jr.fa) and
  the post-kmercount-filter count map for the J-R region.

### `region_sequence_limit` experiments

Added two diagnostic knobs that change the loose default
`region_len + peak_scan + k_size`:

- `KESTREL_MED_SEQ_LIMIT=1` (limit = `region_len + peak_scan`): 6818
  records (extra=2481, missing=560, ins=409 vs Java's 390 — closest yet
  to Java's insertion count). Insertions drop from 1300 to 409 with
  this knob.
- `KESTREL_TIGHT_SEQ_LIMIT=1` (limit = `region_len`): 6532 records
  (extra=2185, missing=550, ins=0). Insertions vanish entirely.

So Java's natural addBase-driven exit appears to cap consensus length
near `region_len + peak_scan` for this dataset. The default
`region_len + peak_scan + k_size` ceiling is too loose by ~20 bases and
that extra rope is exactly what fuels Rust's deletion-like haplotype
traversal through MUC1 repeats. These knobs are off by default so the
existing N-S regression test (`graph_haplotypes_recovers_reduced_vntyper_ns_insertion_branch`)
remains in scope; once the root cause of the over-extension is found
they should become unnecessary.

### Focused diagnostic test

`crates/kestrel/tests/jr_traversal.rs` runs the haplotype graph
assembly on the real post-`kmercount:5` k-mer count map for the
negative VNtyper FASTQ (committed as
`crates/kestrel/tests/fixtures/jr_counts.tsv`, ~603 KB, 25,299 unique
k-mers). It assembles `J-R:4-119` and asserts that the result should
match Java's 0 haplotypes. The test currently fails with 15 Rust
haplotypes. Gate it behind `KESTREL_RUN_JR_DIAGNOSTIC=1` so it does
not block normal test runs.

The new `[KDBG-ITER-END]` trace adds per-iter consensus length,
max-alignment-score, and saved-state stack-size at the end of each
outer iter (first 5 only). For J-R:4-119 in isolation:

```
iter=1 consensus_len=80  max_align_score=536  stack_size=10
iter=2 consensus_len=117 max_align_score=940  stack_size=10
iter=3 consensus_len=100 max_align_score=728  stack_size=10
iter=4 consensus_len=117 max_align_score=980  stack_size=10
iter=5 consensus_len=117 max_align_score=960  stack_size=10
```

First successful emit lands at iter 481 with `consensus_len=117`.
So early iters all build chain entries that `trim_haplotypes` removes
(consensus does not end-anchor on `ref[100..120]`), and the saved-
state stack never drains. The cycle break rate per iter is ~4.7% in
Rust vs ~29% in Java — Rust's saved-state stack stays churning while
Java's saturates. The investigation has gone as far as code reading,
empirical experiments, and per-iter diagnostics can take it without
side-by-side Java instrumentation. The next step is either to add a
custom JVM agent that prints Java's `maxAlignmentScoreNode` chain
contents per addBase, or to write a Rust-only emulator that mirrors
Java's exact stack and chain-handling and bisects against the
observed Rust trace.

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

## 2026-05-15: Shared `haplotype_built` flag (correct but not sufficient)

### Hypothesis

Java's `MaxAlignmentScoreNode.haplotypeBuilt` mutates via reference semantics:
when `getHaplotypes` walks the chain and sets `node.haplotypeBuilt = true`, the
flag is observed by *every* saved snapshot that references the same node. Rust
deep-cloned the chain via `Box<MaxAlignmentScoreNode>::clone`, so each
snapshot got an isolated `haplotype_built: bool`. The hypothesis: a node
emitted in iter N is re-emitted in iter N+M when its containing chain is
restored from a snapshot taken before iter N.

### Implementation

Changed `MaxAlignmentScoreNode.haplotype_built` from `bool` to
`Rc<Cell<bool>>` (`vendor/rust/kestrel-rs/crates/kestrel/src/align/mod.rs`).
Cloning a `MaxAlignmentScoreNode` now `Rc::clone`s the flag, so every snapshot
of a node observes mutations made by any other snapshot. This matches Java's
reference semantics exactly.

### Result

- Compiles. All node tests still pass.
- `KESTREL_RUN_JR_DIAGNOSTIC=1 cargo test -p kestrel --test jr_traversal`:
  unchanged. Still produces 15 haplotypes for J-R:4-119. Java produces 0.
- Negative VNtyper FASTQ parity: still 7062 actual vs 4897 expected, 2727
  extras, 562 missing — identical to before the fix.

### Why it didn't move parity

The J-R diagnostic counters confirm:

```
[KDBG-BUILD] fwd region J-R:4-119 iters=26894 raw_emits=1753
unique_emitted=1753 container=15
```

`raw_emits == unique_emitted`. The runner-level `emitted` HashSet (keyed by
sequence + cigar) sees zero duplicates. Every one of the 1753 haplotypes
emitted across 26,894 outer iters has a distinct (sequence, cigar). So they
come from 1753 *different* chain terminal nodes — not 1753 re-emissions of
the same node. `haplotype_built` sharing has no effect when every emit is
already a fresh node.

The remaining gap is therefore in **chain generation**, not chain emission.
Rust generates 1753 distinct chain terminal positions; Java generates 0
that survive `trim_haplotypes`. Both `trim_haplotypes` implementations are
byte-equivalent (verified). The divergence is upstream — Rust's outer
iterations explore far more chain configurations than Java's.

### Side-by-side save-event match for first 20 inner iters

Manual comparison of Java's `trace.log` `Saving state` events against Rust's
`[KDBG-CHOOSE]` traces for J-R:4-119 first chain build:

| iter | kmer (start)         | depths (A,C,G,T)        | java saves                      | rust saves                      | match |
| ---- | -------------------- | ----------------------- | ------------------------------- | ------------------------------- | ----- |
| 1.1  | GGGGCGGTGGAGCCCGGGGC | 6, 21382, 1600, 1572    | A(6), G(1600), T(1572)          | A(6), G(1600), T(1572)          | ✓     |
| 1.2  | GGGCGGTGGAGCCCGGGGCC | 5, 35, 21499, 0         | A(5), C(35)                     | A(5), C(35)                     | ✓     |
| 1.3  | GGCGGTGGAGCCCGGGGCCG | 29, 23, 26513, 0        | C(23), A(29)                    | C(23), A(29)                    | ✓     |
| 1.4  | GCGGTGGAGCCCGGGGCCGG | 18, 25154, 1021, 24     | A(18), G(1021), T(24)           | A(18), G(1021), T(24)           | ✓     |
| 1.5  | CGGTGGAGCCCGGGGCCGGC | 12, 26661, 59, 27       | A(12), G(59), T(27)             | A(12), G(59), T(27)             | ✓     |
| 1.6  | GGTGGAGCCCGGGGCCGGCC | 16, 216, 197, 26536     | A(16), G(197), C(216)           | A(16), G(197), C(216)           | ✓     |
| 1.7  | GTGGAGCCCGGGGCCGGCCT | 8, 0, 26633, 0          | A(8)                            | A(8)                            | ✓     |
| 1.8  | TGGAGCCCGGGGCCGGCCTG | 8, 5849, 21662, 21      | A(8), C(5849), T(21)            | A(8), C(5849), T(21)            | ✓     |
| 1.9  | GGAGCCCGGGGCCGGCCTGG | 56, 308, 544, 20471     | A(56), C(308), G(544)           | A(56), C(308), G(544)           | ✓     |

Every save event in the first 9 inner iters matches Java byte-for-byte
(kmer, depth, order). Stack eviction events also match — Java removes
`min=5` before save 12, Rust does the same. The divergence emerges *somewhere
past iter 1.9*, but the per-iter trace shows identical save attempt streams
for the early iters.

### Java's stack drains via rejection; Rust's stays full

Java for J-R:4-119:
- 446 total save attempts.
- 38 accepted (10 initial + 28 evictions).
- 408 rejected (stack at capacity, proposed `min_depth` ≤ stack min).
- 38 outer iters, drained to empty.

Rust for J-R:4-119:
- 164,140 total save attempts.
- 40,582 accepted.
- 123,558 rejected.
- 26,894 outer iters, stack remained at cap=10 throughout.

Reject ratio: Java 91.5%, Rust 75.3%. Rust accepts 3× more frequently per
attempt. With ~1.51 accepts/iter and 1 restore/iter, Rust's net stack growth
is +0.51 per iter — capped at 10 by eviction. Java's net is ~0/iter (1
accept ≈ 1 restore), eventually draining when later iters produce shorter
chains that don't refill saves at the same rate.

The 3× acceptance-rate divergence must come from differences in the
`min_depth` proposed at save time vs the stack's current minimum. But the
first-9-iter trace shows identical proposed `min_depth` values, so the
divergence must emerge later (deeper in the chain, or after a different
restore path is taken).

### Next steps

The chain-building algorithm itself is byte-equivalent for at least the first
20 inner iters. The divergence must emerge later in the same outer iter OR on
the first restore. The remaining instrumentation gap is to **dump Java's
saves for iters 2-10+ and compare against Rust's** — pinning down the exact
inner-iter where Java rejects but Rust accepts, or vice versa. With 26,894
Rust iters vs 38 Java iters, the divergence is somewhere in those first ~38
iters that Java terminates with. After that, Rust's extra iters are purely
exploring paths that Java has already excluded.

## 2026-05-15: Initial `min_depth` and runner-level state dedup

Two more checks ruled out, both no-op for the parity numbers:

### Initial `min_depth` reverse-complement fix

`build_forward_haplotypes` and `build_reverse_haplotypes` initialized
`min_depth` from `counter.get(&kmer)` only — forward strand only. Java does
`counter.get(kmer) + counter.get(revKmer)` when `countReverseKmers` is true.
Switched both call sites to use `kmer_depth(...)` so the initial value adds
the reverse-complement count.

Result: parity numbers unchanged (still 7062 vs 4897 expected, 2727 extras,
562 missing). The initial value is quickly overwritten by lower depths from
chain progression, so the off-by-one start was masked.

### Runner-level `SavedBranchKey` HashSet dedup

`save_alignment_state` keys every save attempt by `(kmer, next_base,
consensus)` and skips duplicates via a `HashSet<SavedBranchKey>` that
persists for the lifetime of the build (never cleared). Java has no such
filter.

Wrapped the dedup in a `KESTREL_DISABLE_STATE_DEDUP=1` opt-out and re-ran:

- J-R diagnostic: identical 26,894 iters, 1753 raw emits, 15 haps.
- Negative parity: identical 7062 vs 4897, 2727 extras, 562 missing.

So the runner-level dedup is *not* the source of the divergence — the
duplicate keys never actually fire in J-R.

### Matrix and weight inspection

Verified Rust vs Java match on:
- `AlignmentWeight` defaults: `match=10, mismatch=-10, gap_open=-40,
  gap_extend=-4, init=0, new_gap=gap_open+gap_extend=-44`.
- Align-table candidate score formula: `source.score + (match or mismatch)`.
- Ref-gap-table candidate scores: `align→ref_gap = +new_gap`, `ref_gap→
  ref_gap = +gap_extend`, `con_gap→ref_gap = +new_gap`.
- Con-gap-table candidate scores: `align_next→con_gap = +new_gap`,
  `ref_gap_next→con_gap = +new_gap`, `con_gap_next→con_gap = +gap_extend`.
- `trace_branch` order: Rust iterates `[align, ref_gap, con_gap]`
  candidates; Java does the same. Tie-broken branches prepend in the same
  order.
- `record_max_node` gating: both use `maxScore >= maxAlignmentScore &&
  maxScore > 0`. `next` is `null/None` if strictly greater, else the
  existing chain head.
- `allow_end_deletion` setting: `left_end || right_end`. For J-R
  diagnostic (start=4, end=100), both ends are bounded so allow_end is
  false in both ports.
- `KmerHashSet.insert` / `KmerHashSet.add`: both return `true` if inserted,
  `false` if already present. No semantic difference.

### Status after this session

Confirmed bug fixes in this session:

1. `MaxAlignmentScoreNode.haplotype_built` now shares its `Cell<bool>`
   across clones via `Rc`, matching Java's reference semantics.
2. Initial `min_depth` now includes the reverse-complement count when
   `count_reverse_kmers` is set.
3. `KESTREL_DISABLE_STATE_DEDUP` env var gates the runner-level
   `SavedBranchKey` HashSet so future investigations can bisect it cleanly.

None of the three closed the parity gap. The numbers are persistently
**7062 actual vs 4897 expected, 2727 extras, 562 missing** on the negative
VNtyper FASTQ test. The extras are biased toward low-GDP records
(`gdp_buckets`: 2-5: 232, 6-20: 523, 21-100: 528 in Rust vs Java's 2, 8, 9)
while the missing are concentrated at GDP=970 high-coverage insertions
(`G→GGGTGGAGCCCGGGGCCGG` repeated across E-N, N-R, O-N, R-M, F-N at
position 26).

### Remaining investigative angles

The algorithm appears textually byte-equivalent in:

- Matrix score formulas (3 tables, 3 source-table transitions each).
- `record_max_node` chain-extension/reset semantics.
- `trace_branch` tie-broken candidate ordering.
- `save_state` rejection and `removeMinState` eviction policies.
- Cycle detection via `KmerHashSet`.
- `kmer_depth` (forward + optional reverse).
- `trim_haplotypes` end-kmer-mismatch removal.
- `get_haplotypes` `haplotype_built` skip-on-rebuild (now shared via Rc).

Yet Rust's J-R outer-iter count is 707× Java's (26894 vs 38), and the
overall variant set differs in both directions (extras + missing). The
divergence must be in:

1. **Matrix data flow across iters.** Specifically, the `matrix_col_*`
   `Vec<Option<TraceNodeRef>>` snapshots at save time — these are deep
   `Vec::clone`d but the inner `Rc<TraceNode>` are shared. Need to verify
   that the matrix state at restore matches Java byte-for-byte (the swap
   of `next` → current happens at end of `add_base`; if the snapshot
   captures before the swap, the matrices look different).
2. **`addBase` return value.** Java's `addBase` returns
   `maxPotScore >= maxAlignmentScore && maxPotScore > 0`. Rust's equivalent
   is the same formula. But `maxPotScore` is accumulated DURING the add_base
   call. If Rust accumulates an extra contribution somewhere Java doesn't
   (e.g., an additional max-of-candidate within the loop), Rust's iter would
   keep returning true longer, leading to longer chains and more saves.
3. **The `record_max_node` call at the deletion bottom-row.** In Java this
   is gated by `allowEndDeletion`; in Rust the same gating exists. But the
   `record_max_node` for the ALIGN-table bottom (line 1124, no gate) fires
   unconditionally on `Some(node)` — if Rust's matrix update produces a
   non-None bottom-row node where Java's is `ZERO_NODE`, Rust would record
   max where Java would not.

Next session pursue (3) — instrument Rust to log
`matrix_col_align_next[ref_length - 1]` per iter, run with
`KESTREL_TRACE_REGION=J-R:4-119`, and check the FIRST iter where Rust's
bottom-row is `Some` while Java's would be ZERO_NODE.

## 2026-05-15: Cycle confirmation in Rust outer iters

Extended `KESTREL_TRACE_ITER_MAX=50` and dumped iter-end stats for J-R:4-119:

```
iter=1  consensus_len=80  max_align_score=536.0  stack_size=10  min_depth=17943
iter=2  consensus_len=117 max_align_score=940.0  stack_size=10  min_depth=1526
iter=3  consensus_len=100 max_align_score=728.0  stack_size=10  min_depth=21
iter=4  consensus_len=117 max_align_score=980.0  stack_size=10  min_depth=21
iter=5  consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=18
iter=6  consensus_len=115 max_align_score=886.0  stack_size=9   min_depth=7
iter=7  consensus_len=100 max_align_score=728.0  stack_size=10  min_depth=1600
iter=8  consensus_len=117 max_align_score=980.0  stack_size=10  min_depth=562
iter=9  consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=80
iter=10 consensus_len=117 max_align_score=940.0  stack_size=10  min_depth=22
iter=11 consensus_len=117 max_align_score=920.0  stack_size=9   min_depth=6
iter=12 consensus_len=108 max_align_score=814.0  stack_size=10  min_depth=222
iter=13 consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=127
iter=14 consensus_len=117 max_align_score=940.0  stack_size=10  min_depth=100
iter=15 consensus_len=118 max_align_score=940.0  stack_size=9   min_depth=17
...
iter=24 consensus_len=107 max_align_score=774.0  stack_size=8   min_depth=6
iter=25 consensus_len=80  max_align_score=536.0  stack_size=10  min_depth=988
iter=26 consensus_len=117 max_align_score=940.0  stack_size=10  min_depth=988
iter=27 consensus_len=100 max_align_score=728.0  stack_size=10  min_depth=21
iter=28 consensus_len=117 max_align_score=980.0  stack_size=10  min_depth=21
iter=29 consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=18
iter=30 consensus_len=115 max_align_score=886.0  stack_size=9   min_depth=7
iter=31 consensus_len=100 max_align_score=728.0  stack_size=10  min_depth=988
iter=32 consensus_len=117 max_align_score=980.0  stack_size=10  min_depth=562
iter=33 consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=80
iter=34 consensus_len=117 max_align_score=940.0  stack_size=10  min_depth=22
iter=35 consensus_len=117 max_align_score=920.0  stack_size=9   min_depth=6
iter=36 consensus_len=108 max_align_score=814.0  stack_size=10  min_depth=222
iter=37 consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=127
iter=38 consensus_len=117 max_align_score=940.0  stack_size=10  min_depth=100
iter=39 consensus_len=118 max_align_score=940.0  stack_size=9   min_depth=17
iter=40 consensus_len=117 max_align_score=960.0  stack_size=10  min_depth=43
```

**Iters 25-40 are a near-perfect structural mirror of iters 1-16**:

| Iter A | Iter B | consensus_len | max_align_score | stack_size |
| ------ | ------ | ------------- | --------------- | ---------- |
| 1      | 25     | 80            | 536.0           | 10         |
| 2      | 26     | 117           | 940.0           | 10         |
| 3      | 27     | 100           | 728.0           | 10         |
| 4      | 28     | 117           | 980.0           | 10         |
| 5      | 29     | 117           | 960.0           | 10         |
| 6      | 30     | 115           | 886.0           | 9          |
| 7      | 31     | 100           | 728.0           | 10         |
| 8      | 32     | 117           | 980.0           | 10         |
| 9      | 33     | 117           | 960.0           | 10         |
| 10     | 34     | 117           | 940.0           | 10         |
| 11     | 35     | 117           | 920.0           | 9          |
| 12     | 36     | 108           | 814.0           | 10         |
| 13     | 37     | 117           | 960.0           | 10         |
| 14     | 38     | 117           | 940.0           | 10         |
| 15     | 39     | 118           | 940.0           | 9          |

Only `min_depth` differs between corresponding rows. Everything else matches.

This is conclusive evidence that **Rust's saved-state stack is recycling
the same kmer/consensus configurations**, leading to the same chain
shapes being rebuilt across cycles. Java's stack drains after 38 iters
because its saves saturate (no more new save opportunities); Rust's saves
keep refilling because each "cycle" of 15-16 iters produces enough new
saves to keep the stack at ~10.

The MUC1 reference is highly repetitive, so different alt-branch kmers
do converge to the same chain shapes. Java's algorithm avoids this somehow
— either by saving fewer alt branches or by skipping branches that lead
to already-explored configurations.

### Hypothesis for the root cause

Rust's `KmerHashSet` cycle detection works per-outer-iter (each restore
gets a fresh clone of kmer_hash from save time). So within a single outer
iter, repeated kmers are caught. But across outer iters, the cycle
detection doesn't apply — iter 25's path can re-traverse kmers visited
by iter 1.

Java's algorithm somehow doesn't have this property — maybe Java's
saveState stores fewer alt branches, or Java's chain-extension semantics
differ on score ties.

Concrete next-session experiment: **Force Rust's `saved_states` HashSet
dedup to be per-outer-iter** (clear at start of each outer iter or hash
by chain-shape rather than kmer/consensus). Currently the HashSet
persists for the entire build but only keys by (kmer, next_base,
consensus) — it doesn't catch alt branches that lead to the same chain
shape via different intermediate kmers.

Alternatively, **add a chain-shape dedup at the haplotype emission
level**: when emit produces a haplotype whose final (chain_length,
chain_score, end_kmer) matches a previous emission's, skip it. This
would catch the cycle without requiring deeper algorithm changes.

### Aggressive dedup experiment (KESTREL_AGGRESSIVE_STATE_DEDUP=1)

Tested by hashing save keys by `(kmer, next_base)` only — dropping the
consensus suffix that currently distinguishes alt branches converging at
the same kmer via different intermediate paths.

J-R diagnostic results (with aggressive dedup):
- Outer iters: 26,894 → 283 (99% reduction — cycle confirmed).
- Raw emits: 1753 → 11.
- Haplotypes: 15 → 11.

Full negative VNtyper FASTQ parity (with aggressive dedup):
- Actual records: 7062 → 9359 (WORSE, 32% more).
- Extras: 2727 → 5020.
- Missing: 562 → 558 (mostly unchanged).

**Conclusion**: The cycle hypothesis is confirmed for J-R, but
aggressive dedup is the wrong fix. It prunes save attempts that are
legitimately distinct in OTHER regions, causing different chains to win
the eviction race and producing different (often worse) variant calls.
The right fix must distinguish "cycle-driving alt branches" from
"legitimate distinct alt branches with different downstream consensus",
which requires algorithm-level instead of save-key-level discrimination.

### Session summary (2026-05-15)

**Confirmed and committed fixes** (none alone closes parity gap):

1. `MaxAlignmentScoreNode.haplotype_built` now `Rc<Cell<bool>>` —
   shared across snapshot clones to match Java's reference semantics.
2. Initial `min_depth` in `build_forward_haplotypes` and
   `build_reverse_haplotypes` adds the reverse-complement count to
   match Java's `countReverseKmers` behavior.
3. Three opt-in escape-hatch env vars for future investigations:
   `KESTREL_DISABLE_STATE_DEDUP`, `KESTREL_AGGRESSIVE_STATE_DEDUP`,
   `KESTREL_TRACE_ITER_MAX`.

**Bug status**: The parity gap remains at 7062 actual vs 4897 expected,
2727 extras, 562 missing. The root cause is identified as Rust's
outer-iter cycle on repetitive regions like MUC1 J-R: saved alt branches
converging at the same kmer via different consensus paths cause the
saved-state stack to refill faster than it drains, producing 700×
more outer iters than Java.

**Confirmed byte-equivalent vs Java** (extensive verification this session):
- All `AlignmentWeight` defaults and derived values.
- All matrix transition scores (align, gap_ref, gap_con tables).
- `trace_branch` candidate ordering.
- `record_max_node` chain-extension/reset semantics.
- `save_state` rejection and `removeMinState` eviction policies
  (verified for J-R iters 1-9).
- Cycle detection via `KmerHashSet`.
- `kmer_depth` (forward + reverse).
- `trim_haplotypes` end-kmer-mismatch removal.
- `get_haplotypes` `haplotype_built` skip-on-rebuild (now shared via Rc).
- `addBase` true/false return formula.

**Next session priority**: Either (a) JVM-side instrument Java to dump its
exact saved-state stack contents per inner-iter and bisect against Rust's,
OR (b) attempt a targeted fix in `record_max_node` that detects "alt branches
producing the same trace-shape tail" and skips chain extension when the new
node's trace_node tail matches an already-emitted node's tail (a
chain-shape-aware variant of the `haplotype_built` flag).
