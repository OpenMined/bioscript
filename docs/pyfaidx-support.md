# `bioscript.pyfaidx` Support Matrix

Import form:

```python
from bioscript import pyfaidx
```

This module is a BioScript-supported subset of `pyfaidx`, backed by Rust native
FASTA/FAI access.

## First Slice

| API | Status | Notes |
| --- | --- | --- |
| `pyfaidx.Fasta(path)` | Initial Rust support | `bioscript-libs` can load local FASTA contents with `Fasta::from_path`; runtime/Python constructor binding is still pending. |
| `fasta["22"]` | Initial Rust support | `bioscript-libs` can look up loaded contigs by name. Runtime/Python `[]` binding is pending. |
| `fasta["22"][start:stop]` | Initial Rust support | `FastaRecord::slice` implements 0-based exclusive slicing. Runtime/Python `[]` binding is pending. |
| `str(fasta["22"][start:stop])` | Planned | Python wrapper/runtime conversion still pending. |

## Explicitly Unsupported Initially

| API | Behavior |
| --- | --- |
| FASTA mutation/write APIs | Return unsupported feature error. |
| Remote FASTA URLs | Return unsupported feature error unless a future sandbox policy allows them. |
| Indexed large FASTA access | Deferred; current Rust scaffold loads local FASTA contents directly. |
| Full `pyfaidx.Sequence` behavior | Deferred until needed by assays. |

## Test Sources

Use upstream `pyfaidx` source and tests as reference material under
`vendor/python/pyfaidx` once vendored. Port focused tests for:

- contig lookup
- slicing coordinate behavior
- string conversion
- out-of-bounds and invalid slice errors
