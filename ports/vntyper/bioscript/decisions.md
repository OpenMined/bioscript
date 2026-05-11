# VNtyper BioScript Port Decisions

## Public API Shape

Use a step-oriented API for the port internals:

- `vntyper_regions.region_string(...)`
- `vntyper_commands.plan_bam_pipeline(...)`
- `vntyper_port.process_kestrel_vcf(...)`
- `vntyper_port.build_report_json(...)`

A later `vntyper.run(config)` convenience wrapper can call these steps once the
minimal BAM path has parity. The step-oriented shape keeps tests focused and
lets BioScript expose only the native/library surface needed by each stage.

## Kestrel Resolution

Use the vendored Kestrel source under `ports/vntyper/kestrel` as the reference,
but do not assume a built JAR exists there. The first runnable adapter accepts a
configured JAR path and defaults command plans to:

```text
ports/vntyper/kestrel/kestrel.jar
```

The native Rust Kestrel spike comes after external-tool parity.

## Table Operations

Keep pandas-like operations VNtyper-local for now. The first BioScript port uses
plain lists of dictionaries and small helper functions. Add a shared
`bioscript.table` module only if another port needs the same operations or the
VNtyper implementation starts duplicating generic table logic.

## References

Read the VNtyper MUC1 motif reference from the upstream submodule for the first
milestone:

```text
ports/vntyper/vntyper/reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa
```

Copy references into BioScript-owned fixtures only for tiny deterministic tests
or if upstream reference layout becomes unstable.
