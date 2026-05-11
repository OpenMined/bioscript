# Python Parity Testing

BioScript library shims should be testable from normal Python and from the
BioScript runtime. The goal is to let authors prototype with the same import
shape that BioScript supports:

```python
from bioscript import pysam
```

## Backends

The future Python package should support three backend modes:

| Backend | Purpose |
| --- | --- |
| `rust` | Use the Rust native shim exposed through PyO3 or an equivalent extension. |
| `python` | Delegate to the real Python library, such as installed `pysam`, when available. |
| `auto` | Prefer Rust native shim, fall back only where explicitly allowed by tests. |

Backend selection can be controlled by an environment variable such as:

```text
BIOSCRIPT_BACKEND=rust
BIOSCRIPT_BACKEND=python
BIOSCRIPT_BACKEND=auto
```

## Test Strategy

Each compatibility test should run the same high-level case against every
available backend:

1. Real Python library, when installed.
2. Python package using Rust native shim.
3. BioScript/Monty runtime using `from bioscript import ...`.

Tests should compare observable behavior, not internal implementation details.
For example, a `pysam.AlignmentFile.fetch` parity test should compare read
coordinates and selected read attributes for a tiny fixture region.

## Upstream Tests

Do not gate BioScript on entire upstream suites at first. Instead:

- vendor upstream source for reference
- identify the smallest upstream tests that cover supported APIs
- port focused tests into BioScript-owned test files
- link comments back to upstream files or test names

This keeps compatibility deliberate and avoids accidentally promising the whole
surface of large libraries such as `pysam`.

