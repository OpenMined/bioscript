# Variant Catalogue Schema

For the assay/runtime design, including direct catalogue members and Python
analysis assets, see [variant-catalogue.md](variant-catalogue.md).

This schema describes a compact, auditable source catalogue for projects with
many variants. It is intended for curation and packaging. Current BioScript
runtime lookup still consumes normal `bioscript:variant:1.0` variant manifests;
a packaging step can expand this catalogue into per-variant YAML files.

Use `bioscript:variant-catalogue:1.0` when a project has more than about 10
variants, or when repeated provenance and per-file boilerplate make review
harder than a tabular source.

## Shape

```yaml
schema: "bioscript:variant-catalogue:1.0"
version: "1.0"
name: "thalassemia-variants"

variants:
  source: "variants.tsv"
  format: "tsv"
  key: "variant_id"
  columns:
    variant_id:
      role: "id"
      required: true
    rsid:
      role: "identifier.rsid"
    alts:
      role: "alleles.alts"
      list_separator: "|"
    grch38_pos:
      role: "coordinates.grch38.pos"
      type: "integer"

findings:
  source: "findings.tsv"
  format: "tsv"
  key: "variant_id"
  columns:
    variant_id:
      role: "variant.id"
      required: true
    finding_id:
      role: "finding.id"
    alt:
      role: "finding.alt"
    notes:
      role: "finding.notes"

provenance:
  sources:
    - id: "ithagenes"
      kind: "database"
      label: "IthaGenes"
      url: "https://www.ithanet.eu/db/ithagenes?action=list"
    - id: "dbsnp"
      kind: "database"
      label: "dbSNP"
      url_template: "https://www.ncbi.nlm.nih.gov/snp/{rsid}"
```

## Required Fields

- `schema`: must be `bioscript:variant-catalogue:1.0`
- `version`: must be `1.0`
- `name`
- `variants.source`

`variants.format` and `findings.format`, when present, must be `tsv`.
`variants.columns` and `findings.columns` are recommended for real catalogues so
tools can validate and interpret TSV columns without relying on hardcoded column
names.

## Recommended Files

Keep files separate during curation:

- `variants.yaml`: catalogue manifest and shared provenance
- `variants.tsv`: biological identity and matching fields
- `findings.tsv`: one interpretation row per finding

Generate per-variant `bioscript:variant:1.0` YAML as a packaging/build artifact.

## `variants.tsv`

Recommended columns:

```text
variant_id	name	gene	rsid	aliases	kind	ref	alts	observed_alts	grch37_chrom	grch37_pos	grch37_start	grch37_end	grch38_chrom	grch38_pos	grch38_start	grch38_end
```

Use `|` inside list cells such as `aliases`, `alts`, and `observed_alts`.
Use `pos` columns for SNVs and `start`/`end` columns for spans.

Declare the semantic mapping in `variants.columns`. Recommended roles include:

```text
id
name
gene
identifier.rsid
identifier.aliases
alleles.kind
alleles.ref
alleles.alts
alleles.observed_alts
coordinates.grch37.chrom
coordinates.grch37.pos
coordinates.grch37.start
coordinates.grch37.end
coordinates.grch38.chrom
coordinates.grch38.pos
coordinates.grch38.start
coordinates.grch38.end
```

## `findings.tsv`

Recommended columns:

```text
variant_id	finding_id	schema	alt	label	summary	notes
```

Additional source-specific columns such as `itha_id`, `functionality`,
`phenotype`, `transcript_hgvs`, `genomic_hgvs`, and `ncbi_spdi` are allowed as
curation data. The packaging tool can combine those atoms into standard variant
finding entries.

Each `findings.tsv.variant_id` should match a row in `variants.tsv`. Each
allele-specific `alt` should be present in the corresponding variant row's
`alts`, unless `alt` is `*`.

Declare the semantic mapping in `findings.columns`. Recommended roles include:

```text
variant.id
finding.id
finding.schema
finding.alt
finding.label
finding.summary
finding.notes
source.itha_id
source.functionality
source.phenotype
source.transcript_hgvs
source.genomic_hgvs
source.ncbi_spdi
```
