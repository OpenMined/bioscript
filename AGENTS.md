# BioScript Agent Notes

## Source Size Heuristic

Keep first-party production Rust source files at or below 500 lines. This applies
to files under `rust/bioscript-*/src/**/*.rs`.

The 500-line rule does not apply to:

- integration tests and unit-test modules
- vendored code and patched path dependencies
- generated code, if any is added later

Put substantial test coverage in separate test modules under `tests/` so the
production limit measures production code, not test scaffolding. Test files
should still be split when they mix unrelated behavior or become hard to scan.

When a production file grows past 500 lines, split it before adding more
behavior. Temporary exceptions must be listed in this file under
`Current Refactor Backlog`; the source-size guard reads that list and fails when
it drifts from the code.
