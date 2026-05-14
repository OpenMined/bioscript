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
