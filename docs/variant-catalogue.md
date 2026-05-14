# Variant Catalogue Design Note

Variant catalogues are for assays and panels that need to test many variants
without making reviewers inspect hundreds of near-identical YAML files. The
curated source stays as plain text TSV plus a small YAML manifest. Build tools
can still expand the catalogue into normal `bioscript:variant:1.0` YAML files
when a package or older runtime path needs per-variant manifests.

This is a design note for the intended BioScript flow. The catalogue manifest
schema already exists as `bioscript:variant-catalogue:1.0`; direct assay member
support for catalogues is the next BioScript schema/runtime change.

## Why This Exists

Small assays are easiest to review as hand-written variant YAML:

```yaml
members:
  - kind: "variant"
    path: "variants/rs73885319.yaml"
    version: "1.0"
```

That breaks down for catalogues such as thalassemia, where IthaGenes/dbSNP style
sources can produce hundreds of rows. Repeating the same provenance, source
URLs, and schema boilerplate in every YAML file makes review noisier. A TSV is
more compact, can be sorted, can be checked with command-line tools, and keeps
the source curation auditable.

The runtime still needs structured observations and controlled file access. The
intended flow is:

1. BioScript reads the assay.
2. BioScript expands catalogue members into variant lookup work.
3. BioScript gathers one observation row per variant for the participant.
4. BioScript mounts a small virtual filesystem for the analysis script.
5. BioScript writes observations to a conventional work path.
6. BioScript runs the analysis Python with a structured `bioscript.context`.
7. Python reads TSV inputs with `bioscript.read_tsv(path)`.
8. Python joins observations to attached TSV assets and emits report rows.

## Required BioScript Schema Change

Assay and panel `members` should accept a catalogue member:

```yaml
members:
  - kind: "variant-catalogue"
    path: "catalogue/variants.yaml"
    version: "1.0"
```

Today, `bioscript:assay:1.0` assay members accept `kind: "variant"` only, so the
direct member above is the required next schema change. The member should point
to the catalogue YAML manifest, not directly to `variants.tsv`, because the YAML
manifest holds schema identity, version, table locations, keys, and shared
provenance.

The catalogue manifest shape is:

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
    name:
      role: "name"
    gene:
      role: "gene"
    rsid:
      role: "identifier.rsid"
    aliases:
      role: "identifier.aliases"
      list_separator: "|"
    kind:
      role: "alleles.kind"
    ref:
      role: "alleles.ref"
    alts:
      role: "alleles.alts"
      list_separator: "|"
    observed_alts:
      role: "alleles.observed_alts"
      list_separator: "|"
    grch37_chrom:
      role: "coordinates.grch37.chrom"
    grch37_pos:
      role: "coordinates.grch37.pos"
      type: "integer"
    grch38_chrom:
      role: "coordinates.grch38.chrom"
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
    schema:
      role: "finding.schema"
    alt:
      role: "finding.alt"
    label:
      role: "finding.label"
    summary:
      role: "finding.summary"
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

## Column Mapping

The manifest should declare how TSV columns map to BioScript concepts. This
keeps the TSV compact while avoiding hidden assumptions about column names.
Tools may provide default mappings for conventional columns, but real catalogues
should encode the mapping explicitly.

Each entry under `variants.columns` or `findings.columns` is keyed by the TSV
column name:

```yaml
variants:
  source: "variants.tsv"
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
```

Column mapping fields:

- `role`: semantic target for this TSV column
- `required`: whether validation should fail if the column is missing or empty
- `type`: optional primitive validation, such as `string` or `integer`
- `list_separator`: separator for list-valued cells, normally `|`

Recommended `variants.columns` roles:

- `id`
- `name`
- `gene`
- `identifier.rsid`
- `identifier.aliases`
- `alleles.kind`
- `alleles.ref`
- `alleles.alts`
- `alleles.observed_alts`
- `coordinates.grch37.chrom`
- `coordinates.grch37.pos`
- `coordinates.grch37.start`
- `coordinates.grch37.end`
- `coordinates.grch38.chrom`
- `coordinates.grch38.pos`
- `coordinates.grch38.start`
- `coordinates.grch38.end`

Recommended `findings.columns` roles:

- `variant.id`
- `finding.id`
- `finding.schema`
- `finding.alt`
- `finding.label`
- `finding.summary`
- `finding.notes`
- `source.itha_id`
- `source.functionality`
- `source.phenotype`
- `source.transcript_hgvs`
- `source.genomic_hgvs`
- `source.ncbi_spdi`

Unknown columns may still exist as source curation data, but executable tools
should only depend on columns declared in the manifest. Validation should check
that required mapped columns are present in the TSV header, integer columns parse
as integers when populated, and list-valued columns are split consistently.

The assay should also allow analysis assets, which are files made available to
the Python analysis by id:

```yaml
analyses:
  - id: "thalassemia_status"
    kind: "bioscript"
    path: "thalassemia.py"
    output_format: "tsv"
    derived_from:
      - "catalogue/variants.yaml"
    assets:
      - id: "variants"
        path: "catalogue/variants.tsv"
      - id: "findings"
        path: "catalogue/findings.tsv"
    emits:
      - key: "thalassemia_status"
        label: "Thalassemia status"
        value_type: "string"
        format: "badge"
```

`assets` are not variant members. They are analysis inputs. The variant
catalogue member tells BioScript what to look up; the pipeline files and assets
give Python the extra curation tables it needs to classify the observed
variants.

## Analysis File Sandbox

Analysis scripts should not receive arbitrary host filesystem paths. BioScript
should mount a small virtual filesystem with conventional paths:

```text
/input
/input/pipeline
/work
/output
```

Root meanings:

- `/input`: read-only run inputs, starting with the genotype/sample input file
- `/input/pipeline`: read-only assay package files, including YAML, TSV assets,
  and analysis code
- `/work`: read/write runtime workspace for BioScript-generated intermediate
  files such as observations
- `/output`: write/preserved output files and auxiliary artifacts

Initial conventional paths:

```text
/input/genotypes
/input/participant.tsv
/input/pipeline/assay.yaml
/input/pipeline/catalogue.yaml
/input/pipeline/variants.tsv
/input/pipeline/findings.tsv
/work/observations.tsv
/work/analysis.jsonl
/output/results.tsv
```

`/input/participant.tsv` and `/work/analysis.jsonl` are optional future
conventions. They can be absent; scripts should check with `bioscript.exists()`.

Access rules:

- `/input` and `/input/pipeline` are read-only
- `/work` is read/write and intended for intermediate files
- `/output` is write/preserved; files here are report artifacts
- host paths, path traversal, and network paths are not available to analysis
  scripts

`/output/results.tsv` is the canonical report output when the script writes a
file. Other files under `/output` are auxiliary artifacts. If a future runtime
captures returned rows/dicts directly, it should reject ambiguous scripts that
both return a primary result and write `/output/results.tsv`.

## `bioscript.context`

BioScript should expose a structured context object so scripts can interrogate
the run without hardcoding every path:

```python
bioscript.context["participant_id"]
bioscript.context["input_files"]
bioscript.context["pipeline_files"]
bioscript.context["assets"]
bioscript.context["observations_file"]
bioscript.context["output_file"]
```

Example context shape:

```python
{
    "participant_id": "genome_hu50B3F5_v5_Full",
    "input_files": {
        "genotypes": "/input/genotypes",
    },
    "pipeline_files": {
        "assay": "/input/pipeline/assay.yaml",
        "catalogue": "/input/pipeline/catalogue.yaml",
    },
    "assets": {
        "variants": "/input/pipeline/variants.tsv",
        "findings": "/input/pipeline/findings.tsv",
    },
    "observations_file": "/work/observations.tsv",
    "output_file": "/output/results.tsv",
}
```

Compatibility variables may still be injected during migration:

```python
participant_id = bioscript.context["participant_id"]
observations_file = bioscript.context["observations_file"]
output_file = bioscript.context["output_file"]
asset_paths = bioscript.context["assets"]
```

The context object should be read-only from the script's point of view. It can
also expose helper methods later, but a dictionary shape is enough for the first
implementation.

## Thalassemia Use Case

The thalassemia project is the motivating case:

- `catalogue/variants.yaml` stores shared catalogue metadata and provenance.
- `catalogue/variants.tsv` stores variant identity, rsIDs, aliases, alleles,
  coordinates, and source curation columns.
- `catalogue/findings.tsv` stores what a variant contributes to: for example
  alpha-globin, beta-globin, or other globin/modifier categories, with source
  notes such as IthaGenes functionality and phenotype.
- `thalassemia.py` reads `observations_file` plus the TSV assets and emits a
  compact classification row.

This lets the assay report a derived classification such as:

- no catalogue thalassemia variant observed
- possible beta carrier
- possible biallelic HBB risk
- possible alpha silent carrier or trait
- possible HbH/Hb Bart spectrum

Those labels require careful clinical curation. The important first step is the
data path: variants become observations, then Python joins observations to
attached TSV assets and emits rows. The exact classification rules need separate
review and should not be inferred from row count alone.

## Minimal Python Analysis

This is intentionally a hello-world style analysis. It proves how the injected
files work, but the real thalassemia logic needs proper allele, phase, gene, and
evidence handling.

```python
import bioscript


observations = bioscript.read_tsv(bioscript.context["observations_file"])
variants = bioscript.read_tsv(bioscript.context["assets"]["variants"])
findings = bioscript.read_tsv(bioscript.context["assets"]["findings"])

# Minimal proof of shape only. Real work should:
# - match observations to variant_id via rsid, aliases, and coordinates
# - verify the observed allele matches the finding allele
# - handle deletions, insertions, complex alleles, no-calls, and phase
# - classify alpha-globin, beta-globin, and modifier findings separately
# - preserve provenance and source notes in emitted rows

if len(observations) > 0 and len(variants) > 0 and len(findings) > 0:
    status = "catalogue_assets_loaded"
else:
    status = "catalogue_assets_missing_or_empty"

bioscript.write_tsv(bioscript.context["output_file"], [
    {
        "participant_id": bioscript.context["participant_id"],
        "thalassemia_status": status,
    }
])
```

The analysis script receives virtual file paths as strings through
`bioscript.context`. `observations_file` is a TSV created by BioScript for the
current participant. `assets` is a dictionary whose keys are the ids from
`analyses[].assets`.

BioScript should provide `bioscript.read_tsv(path)` as a built-in helper. TSV is
part of the catalogue asset contract, so every analysis should not need to carry
its own parser. The helper should return a list of dictionaries keyed by the TSV
header row, use tab as the delimiter, preserve empty cells as empty strings, and
skip empty lines. Later versions may add typed parsing based on catalogue column
metadata, but the default helper should stay predictable and string-based.

## Packaging And Review

Keep these files separate during curation:

```text
catalogue/variants.yaml
catalogue/variants.tsv
catalogue/findings.tsv
```

Generated per-variant YAML files are useful as a packaging artifact, especially
while older runtime paths still expect `kind: "variant"` members. They should not
be the primary review surface for large catalogues.

The long-term runtime path should be better than generated YAML for large
catalogues: BioScript can stream or batch the TSV-backed catalogue directly into
lookups, then pass observations and assets to the analysis script.
