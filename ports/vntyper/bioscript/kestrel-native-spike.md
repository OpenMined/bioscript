# Kestrel Native Rust Feasibility Spike

Outcome: keep Kestrel behind the external JVM adapter for the first runnable
VNtyper BioScript milestone.

## Evidence

- Vendored Kestrel source is present at `ports/vntyper/kestrel`.
- Source size is non-trivial: 91 Java files and about 24,955 lines under
  `ports/vntyper/kestrel/src`.
- Main package areas include:
  - `counter`
  - `activeregion`
  - `align`
  - `refreader`
  - `runner`
  - `variant`
  - `varfilter`
  - `writer/vcf`
  - `hapwriter/sam`
- The repository includes an Ant `build.xml` and dependency JARs under `lib`,
  but there is no built `ports/vntyper/kestrel/kestrel.jar` in the submodule.
- No Kestrel Java test source files were found in the vendored tree.
- The BioScript side does not yet have large expected VNtyper Kestrel VCF/TSV
  outputs for regression comparison.

## Decision

Do not start a native Rust Kestrel port yet.

The external adapter is the practical first target because it lets BioScript
validate the VNtyper pipeline contract before reimplementing a large local
assembly and variant-calling engine. A native port should happen only after the
external-tool-backed path has parity fixtures that can detect behavioral drift.

## Native-Port Entry Points Later

If/when parity fixtures exist, port in this order:

1. `counter`: k-mer count representation and lookup.
2. `refreader`: reference window parsing for the VNTR motif dictionary.
3. `activeregion`: active-region detection and haplotype candidates.
4. `align`: bounded alignment with VNtyper's Kestrel settings.
5. `variant`: insertion/deletion/SNV call representation.
6. `writer/vcf`: reproduce the exact VCF fields consumed by VNtyper.
7. `hapwriter/sam`: reproduce optional SAM output only if report/IGV parity
   requires it.

## Required Before Reopening

- Build or configure a JVM Kestrel JAR for local integration tests.
- Generate expected `output.vcf`, `output_indel.vcf`, `kestrel_pre_result.tsv`,
  and `kestrel_result.tsv` for at least one positive and one negative fixture.
- Add an integration test that runs the external Kestrel adapter and verifies
  those outputs.
