# `bioscript.pysam` Support Matrix

Import form:

```python
from bioscript import pysam
```

This module is a BioScript-supported subset of `pysam`, backed by Rust native
code. Unsupported APIs should fail with explicit compatibility errors.

## First Slice

| API | Status | Notes |
| --- | --- | --- |
| `pysam.AlignmentFile(path, "rc", reference_filename=...)` | Scaffolded | Rust object and mode validation exist; CRAM fetch backend is pending. |
| `pysam.AlignmentFile(path, "rb")` | Scaffolded | Rust object and mode validation exist; BAM backend is pending. |
| `AlignmentFile.fetch(contig, start, stop)` | Initial CRAM support | Rust and BioScript runtime can stream local CRAM fixtures when `reference_filename` is supplied. |
| `AlignedSegment.query_name` | Scaffolded | Rust field exists. Backend population is pending. |
| `AlignedSegment.reference_name` | Initial CRAM support | Populated from the fetch contig. |
| `AlignedSegment.reference_start` | Initial CRAM support | Converted to pysam-style 0-based start from BioScript alignment records. |
| `AlignedSegment.reference_end` | Initial CRAM support | Populated from BioScript alignment records. |
| `AlignedSegment.query_sequence` | Scaffolded | Rust field exists. Backend population is pending. |
| `AlignedSegment.mapping_quality` | Scaffolded | Rust field exists. Backend population is pending. |
| `AlignedSegment.cigarstring` | Initial CRAM support | Populated from the BioScript alignment CIGAR operations. |
| `AlignedSegment.is_unmapped` | Initial CRAM support | Populated from BioScript alignment records. |
| `AlignedSegment.is_reverse` | Scaffolded | Rust field exists. Backend population is pending. |

## Explicitly Unsupported Initially

| API | Behavior |
| --- | --- |
| Write modes such as `"w"`, `"wb"`, `"wc"` | Return unsupported mode error. |
| Mutating reads or headers | Return unsupported feature error. |
| Remote URLs | Return unsupported feature error unless a future sandbox policy allows them. |
| Tags and auxiliary fields | Return unsupported feature error until needed by assays. |
| Full pileup API | Deferred until read iteration and APOL1 parity are stable. |
| Full htslib compatibility | Not a goal for the first slice. |

## Test Sources

Use upstream `pysam` source and tests as reference material under
`vendor/python/pysam` once vendored. Port focused tests for:

- `AlignmentFile.fetch` region behavior
- coordinate conventions
- read attribute names and values
- unsupported mode behavior
