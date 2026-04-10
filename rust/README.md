# Rust BioScript

This workspace is the starting point for a Rust-native BioScript runtime.

The first slice includes:
- domain primitives for DNA/genotype handling
- APOL1 classification logic
- a Monty-backed host runtime that exposes safe genomics functions into sandboxed Python-like code

This is intentionally separate from the existing Python implementation so the
secure execution path can evolve independently.

