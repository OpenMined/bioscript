# Observation Schema

This schema defines the normalized row format for a single BioScript variant observation.

Use it when BioScript needs to emit results that are:

- easy to store in a database table
- easy to write as TSV
- easy to aggregate across cohorts, assays, and input file types
- detailed enough to generate human reports later

The observation row answers: "What did we observe for this participant at this assay variant?"

It does not answer the higher-level interpretation question. Assay-specific results such as APOL1 `G1/G0` should be emitted as a separate assay interpretation that derives from one or more observations.

## Schema Identity

```yaml
schema: "bioscript:observation:1.0"
version: "1.0"
```

For row-oriented outputs, the `schema` and `version` can be file-level metadata rather than repeated in every TSV row. For JSON Lines, include both fields per object.

## Minimal Example

```yaml
participant_id: "NA12878"
assay_id: "APOL1"
assay_version: "1.0"
variant_key: "G2_SITE"

rsid: "rs71785313"

assembly: "GRCh38"
chrom: "22"
pos_start: 36266000
pos_end: 36266005

ref: "TTATAA"
alt: "<DEL:6>"
kind: "del"

match_status: "found"
coverage_status: "covered"
call_status: "called"

genotype: "0/1"
genotype_display: "ID"
zygosity: "het"

ref_count: 1
alt_count: 1
depth: 38
genotype_quality: null
allele_balance: null

outcome: "variant"

evidence_type: "mpileup"
evidence_raw: "22:36266000-36266005 depth=38 deletion_reads=18 ref_reads=20"

facets: "self_reported_ethnicity=African Caribbean;cohort=CariGenetics"
```

## Columns

Required identity fields:

- `participant_id`: participant/sample identifier.
- `assay_id`: stable assay or panel identifier.
- `assay_version`: assay or panel version.
- `variant_key`: assay-local variant key, stable within the assay.

Variant identity:

- `rsid`: dbSNP rsID, or `null` when unavailable.
- `assembly`: reference assembly label, for example `GRCh37`, `GRCh38`, or `unknown`.
- `chrom`: chromosome or contig name.
- `pos_start`: 1-based inclusive start.
- `pos_end`: 1-based inclusive end.
- `ref`: reference allele or reference sequence.
- `alt`: alternate allele or symbolic alternate such as `<DEL:6>`.
- `kind`: variant class.

Observation status:

- `match_status`: whether a data source record matched the target.
- `coverage_status`: whether the target locus was covered or checkable.
- `call_status`: whether a usable genotype call was produced.

Call fields:

- `genotype`: normalized VCF-style genotype using allele indexes.
- `genotype_display`: display genotype in assay-friendly alleles or symbols.
- `zygosity`: simplified zygosity class.

Quantitative evidence:

- `ref_count`: count of reference-supporting calls or reads.
- `alt_count`: count of alternate-supporting calls or reads.
- `depth`: total usable depth.
- `genotype_quality`: genotype quality, when available.
- `allele_balance`: alternate allele fraction, when available.

Report outcome:

- `outcome`: normalized downstream result bucket for this observation.

Evidence:

- `evidence_type`: source of evidence used for this row.
- `evidence_raw`: compact raw evidence text for audit/debugging.

Facets:

- `facets`: optional semicolon-separated key-value metadata for filtering and aggregation.

## Controlled Vocabularies

### `kind`

Allowed values:

- `snp`
- `ins`
- `del`
- `mnv`
- `sv`
- `cnv`
- `fusion`
- `hybrid`

Use `snp` for single nucleotide substitutions. Use `mnv` for multi-nucleotide substitutions. Use symbolic `alt` values such as `<DEL:6>` when that is clearer than a biological sequence for reporting.

### `match_status`

Allowed values:

- `found`
- `not_found`
- `ambiguous`

`found` means a relevant input record or evidence source matched the target variant. `not_found` means no matching record/evidence was found. `ambiguous` means multiple or conflicting matches were found and the row should not be treated as a confident call.

### `coverage_status`

Allowed values:

- `covered`
- `not_covered`
- `unknown`

For direct-to-consumer text exports, absence of an rsID is usually `unknown` unless the platform's marker design is known well enough to say `not_covered`. For CRAM/BAM/mpileup evidence, zero usable depth at the locus is `not_covered`.

### `call_status`

Allowed values:

- `called`
- `no_call`
- `low_quality`
- `filtered`

`called` means the row has a usable genotype. `no_call` means the locus was observed or considered but no genotype could be assigned. `low_quality` and `filtered` preserve cases where a call-like signal exists but should not be promoted to a confident result.

### `genotype`

Allowed values for first-pass diploid reporting:

- `0/0`
- `0/1`
- `1/1`
- `./.`

Future versions may add phased genotypes and multi-allelic genotypes. Keep `genotype_display` for assay-specific display forms such as `ID`, `DD`, or `AG`.

### `zygosity`

Allowed values:

- `hom_ref`
- `het`
- `hom_alt`
- `unknown`

### `outcome`

Allowed values:

- `variant`
- `reference`
- `no_call`
- `not_covered`
- `unknown`
- `partial`

`outcome` is the reporting bucket for the observation:

- `variant`: alternate allele observed.
- `reference`: reference genotype observed.
- `no_call`: source was considered, but no genotype call was available.
- `not_covered`: locus was not covered or not present.
- `unknown`: insufficient information to classify.
- `partial`: incomplete evidence, such as one allele found where two are expected.

### `evidence_type`

Allowed values:

- `vcf_record`
- `gvcf_block`
- `mpileup`
- `cram_slice`
- `split_reads`
- `read_pairs`
- `imputed`
- `null`

Use `null` when no raw evidence was available or retained.

## TSV Header

Recommended TSV column order:

```text
participant_id	assay_id	assay_version	variant_key	rsid	assembly	chrom	pos_start	pos_end	ref	alt	kind	match_status	coverage_status	call_status	genotype	genotype_display	zygosity	ref_count	alt_count	depth	genotype_quality	allele_balance	outcome	evidence_type	evidence_raw	facets
```

TSV rules:

- represent nulls as empty fields
- keep `facets` as semicolon-separated `key=value` pairs
- escape tabs and newlines in `evidence_raw` before writing TSV
- do not put human prose in `facets`; use `evidence_raw` or a separate report field

## Facets

`facets` is intentionally a single nullable string for database and TSV compatibility.

Format:

```text
key=value;key=value
```

Examples:

```text
self_reported_ethnicity=African Caribbean;cohort=CariGenetics
cohort=TestBatch1
```

Rules:

- keys should use lowercase snake_case
- values should not contain tabs or newlines
- semicolons in values should be percent-encoded as `%3B`
- equals signs in values should be percent-encoded as `%3D`

## Relationship To Assay Interpretations

An observation is not the final assay result when the assay combines multiple variants.

For APOL1:

- `G1_SITE_1` observation records the `rs73885319` call.
- `G1_SITE_2` observation records the `rs60910145` call.
- `G2_SITE` observation records the `rs71785313` deletion call.
- A separate assay interpretation derives the final APOL1 status, such as `G1/G0`, `G2/G1`, or `G2/G2`.

This keeps auditability and reporting clean:

- observations are low-level evidence rows
- interpretations are assay-level conclusions
- rendered reports can show both
