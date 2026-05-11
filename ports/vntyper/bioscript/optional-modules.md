# VNtyper Optional Module Triage

The minimal BioScript VNtyper path remains:

```text
BAM -> MUC1 read extraction -> Kestrel VCF -> classification -> TSV/JSON/HTML
```

Optional upstream modules are intentionally not part of the first runnable path.

## adVNTR

Status: defer execution, keep report surface.

Reasoning:
- Upstream treats adVNTR as an independent confirmation caller.
- The BioScript report JSON can already carry `advntr_variants`, compute an
  adVNTR algorithm result, and emit a cross-match summary.
- Running adVNTR needs its own external tool/reference setup and expected test
  outputs.

Next work:
- Add an external `bioscript.advntr` command planner only after the Kestrel BAM
  path has parity.
- Add tiny adVNTR row fixtures for report-only tests.
- Add integration tests only when adVNTR references and outputs are available.

## SHARK

Status: defer.

Reasoning:
- SHARK is not required for the core MUC1 frameshift call.
- It adds another external dependency and output contract before the primary
  Kestrel path is proven.

Next work:
- Read upstream `vntyper/modules/shark`.
- Document the exact command/API surface.
- Decide whether it belongs in BioScript libs or remains an external wrapper.

## Cohort Summaries

Status: defer until single-sample parity.

Reasoning:
- Cohort output depends on stable per-sample JSON/TSV contracts.
- Building it before single-sample parity would lock in unstable report fields.

Next work:
- Define a stable single-sample report schema.
- Add a pure Python/BioScript aggregation helper over report JSON files.

## Mutation Counter

Status: defer.

Reasoning:
- It is not needed for the minimal pathogenic frameshift classification path.
- It should be evaluated after Kestrel/adVNTR output parity is clear.

Next work:
- Inventory upstream mutation-counter inputs and outputs.
- Add fixture-level tests before adding runtime wrappers.
