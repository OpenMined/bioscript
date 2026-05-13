# Upstream Test Plan

BioScript vendors upstream libraries as reference material, but should port only
focused tests for the compatibility subset it claims to support.

## Vendored Sources

| Project | Path | Use |
| --- | --- | --- |
| `pysam` | `vendor/python/pysam` | Alignment, CRAM/BAM, VCF API reference and targeted test ports. |
| `pyfaidx` | `vendor/python/pyfaidx` | FASTA lookup and slicing API reference and targeted test ports. |

## CLI Reference Sources

Do not vendor `htslib`, `samtools`, or `bcftools` yet. The first compatibility
slice is API-shaped (`from bioscript import pysam` and `pyfaidx`), so upstream
Python tests give the most direct coverage. Add CLI repositories later if one
of these becomes true:

- a failing parity case requires htslib/samtools fixture-generation behavior
- BioScript starts emulating a CLI command surface
- pysam upstream tests require source-level htslib/samtools context that cannot
  be captured in a small BioScript-owned fixture

## `pysam` First Test Candidates

Use `vendor/python/pysam/tests/AlignmentFile_test.py` as the initial
source for parity cases.

Smallest useful targets:

- `BasicTestBAMFromFetch.setUp`: open `AlignmentFile(..., "rb")` and call
  `list(self.samfile.fetch())`.
- `BasicTestBAMFromFetch.testARqname`: read `query_name`.
- `BasicTestBAMFromFetch.testARpos`: read `reference_start`.
- `BasicTestBAMFromFetch.testARmapq`: read `mapping_quality`.
- `BasicTestBAMFromFetch.testARcigarstring`: read `cigarstring`.
- `BasicTestBAMFromFetch.testARseq`: read `query_sequence`.
- Region fetch comparisons around `fetch('chr1', start=1000, end=2000)`.

These tests should be ported to tiny BioScript-owned fixtures rather than
depending on the full upstream test harness.

## `pyfaidx` First Test Candidates

Use `vendor/python/pyfaidx/tests/test_feature_bounds_check.py` as the
initial source for parity cases.

Smallest useful targets:

- `test_blank_string`: `seq[0:0]` returns an empty string.
- `test_slice_from_beginning`: first bases through `[:4]`.
- `test_fetch_reversed_coordinates`: reversed coordinates fail.
- `test_fetch_keyerror`: missing contig fails.

The current Rust scaffold already covers blank slices, beginning slices, normal
middle slices, reversed coordinates, and missing contigs against a tiny local
FASTA fixture.
