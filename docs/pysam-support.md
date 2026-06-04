# `bioscript.pysam` Support Matrix

Import form:

```python
from bioscript import pysam
```

This module is a BioScript-supported subset of `pysam`, backed by Rust native
code. Local BAM/CRAM fetches route through `htslib-rs` alignment helpers, and
unsupported APIs should fail with explicit compatibility errors.

## First Slice

| API | Status | Notes |
| --- | --- | --- |
| `pysam.AlignmentFile(path, "rc", reference_filename=...)` | Initial support | Local indexed CRAM fetches use `htslib-rs`; `reference_filename` is required. |
| `pysam.AlignmentFile(path, "rb")` | Initial support | Local indexed BAM fetches use `htslib-rs` associated BAI/CSI lookup. |
| `AlignmentFile.fetch(contig, start, stop)` | Initial BAM/CRAM support | Requires explicit 0-based `start` and half-open `stop`; converts to HTSlib 1-based inclusive regions internally. |
| `AlignedSegment.query_name` | Initial BAM/CRAM support | Populated from the read name when present. |
| `AlignedSegment.reference_name` | Initial BAM/CRAM support | Populated from the fetch contig for mapped reads. |
| `AlignedSegment.reference_start` | Initial BAM/CRAM support | Converted back to pysam-style 0-based start. |
| `AlignedSegment.reference_end` | Initial BAM/CRAM support | Derived from reference-consuming CIGAR operations. |
| `AlignedSegment.query_sequence` | Initial BAM/CRAM support | Populated from the read sequence when present. |
| `AlignedSegment.mapping_quality` | Initial BAM/CRAM support | Populated from the read mapping quality when present. |
| `AlignedSegment.cigarstring` | Initial BAM/CRAM support | Populated from CIGAR operations. |
| `AlignedSegment.is_unmapped` | Initial BAM/CRAM support | Populated from SAM flags. |
| `AlignedSegment.is_reverse` | Initial BAM/CRAM support | Populated from SAM flags. |

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
