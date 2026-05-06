# BioScript Agent Notes

## Source Size Heuristic

Keep first-party production Rust source files at or below 500 lines. This applies
to files under `rust/bioscript-*/src/**/*.rs`.

When editing BioScript Rust, prefer adding behavior to a small, named module
whose filename describes the responsibility. If a file is approaching 500 lines,
split it along a real domain boundary before adding more code. Do not satisfy
the guard by creating arbitrary numbered chunks or `*_part_*` files.

The 500-line rule does not apply to:

- integration tests and unit-test modules
- vendored code and patched path dependencies
- generated code, if any is added later

Put substantial test coverage in separate test modules under `tests/` so the
production limit measures production code, not test scaffolding. Test files
should still be split when they mix unrelated behavior or become hard to scan.

When a production file grows past 500 lines, split it before adding more
behavior. Keep the include list in the parent file short and logical, and leave
file names meaningful enough that future agents can find the right place to edit.
