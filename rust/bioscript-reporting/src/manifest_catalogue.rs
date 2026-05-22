use std::{collections::BTreeMap, path::PathBuf};

use bioscript_core::{GenomicLocus, VariantKind, VariantSpec};
use bioscript_schema::VariantManifest;

use super::{
    ManifestWorkspace, VariantManifestTask, matches_variant_manifest_filters, yaml_string,
    yaml_string_sequence,
};

pub(super) fn load_variant_catalogue_tasks(
    workspace: &impl ManifestWorkspace,
    path: &str,
    filters: &[String],
) -> Result<Vec<VariantManifestTask>, String> {
    let value = workspace.load_yaml(path)?;
    let variants = value
        .get("variants")
        .and_then(serde_yaml::Value::as_mapping)
        .ok_or_else(|| format!("{path} is missing variants table declaration"))?;
    let source = variants
        .get(serde_yaml::Value::String("source".to_owned()))
        .and_then(serde_yaml::Value::as_str)
        .ok_or_else(|| format!("{path} variants.source missing"))?;
    let source_path = workspace.resolve(path, source)?;
    let table_text = workspace.load_text(&source_path)?;
    let rows = parse_tsv_rows(&table_text)?;
    let columns = CatalogueColumns::from_manifest(variants);
    let catalogue_name =
        yaml_string(&value, "name").unwrap_or_else(|| "variant-catalogue".to_owned());
    let catalogue_tags = yaml_string_sequence(&value, "tags")
        .into_iter()
        .filter_map(|value| value.as_str().map(ToOwned::to_owned))
        .collect::<Vec<_>>();

    rows.into_iter()
        .enumerate()
        .map(|(idx, row)| {
            catalogue_row_task(path, &catalogue_name, &catalogue_tags, &columns, idx, &row)
        })
        .filter_map(|result| match result {
            Ok(task)
                if matches_variant_manifest_filters(
                    &task.manifest,
                    &task.manifest_path,
                    filters,
                ) =>
            {
                Some(Ok(task))
            }
            Ok(_) => None,
            Err(err) => Some(Err(err)),
        })
        .collect()
}

fn catalogue_row_task(
    catalogue_path: &str,
    catalogue_name: &str,
    catalogue_tags: &[String],
    columns: &CatalogueColumns,
    idx: usize,
    row: &BTreeMap<String, String>,
) -> Result<VariantManifestTask, String> {
    let variant_id = columns
        .value(row, "id")
        .or_else(|| columns.value(row, "identifier.rsid"))
        .filter(|value| !value.is_empty())
        .ok_or_else(|| {
            format!(
                "{catalogue_path} variants.tsv row {} missing variant id",
                idx + 2
            )
        })?;
    let name = columns
        .value(row, "name")
        .filter(|value| !value.is_empty())
        .map_or_else(
            || format!("{catalogue_name}-{variant_id}"),
            ToOwned::to_owned,
        );
    let mut rsids = split_list(
        columns.value(row, "identifier.rsid"),
        columns.separator("identifier.rsid"),
    );
    rsids.extend(split_list(
        columns.value(row, "identifier.aliases"),
        columns.separator("identifier.aliases"),
    ));
    let mut deduped_rsids = Vec::with_capacity(rsids.len());
    for rsid in rsids {
        if !deduped_rsids.iter().any(|seen| seen == &rsid) {
            deduped_rsids.push(rsid);
        }
    }
    let rsids = deduped_rsids;
    let alternates = split_list(
        columns.value(row, "alleles.alts"),
        columns.separator("alleles.alts"),
    );
    let observed_alternates = {
        let observed = split_list(
            columns.value(row, "alleles.observed_alts"),
            columns.separator("alleles.observed_alts"),
        );
        if observed.is_empty() {
            alternates.clone()
        } else {
            observed
        }
    };
    let spec = VariantSpec {
        rsids,
        grch37: catalogue_locus(columns, row, "grch37")?,
        grch38: catalogue_locus(columns, row, "grch38")?,
        reference: columns.value(row, "alleles.ref").map(ToOwned::to_owned),
        alternate: alternates.first().cloned(),
        observed_alternates,
        kind: columns
            .value(row, "alleles.kind")
            .map(catalogue_variant_kind),
        ..VariantSpec::default()
    };
    if !spec.has_rsids() && !spec.has_coordinates() {
        return Err(format!(
            "{catalogue_path} variants.tsv row {} has no rsid or coordinates",
            idx + 2
        ));
    }
    let mut tags = catalogue_tags.to_vec();
    if let Some(gene) = columns.value(row, "gene").filter(|value| !value.is_empty()) {
        let gene_tag = format!("gene:{gene}");
        if !tags.iter().any(|tag| tag == &gene_tag) {
            tags.push(gene_tag);
        }
    }
    let manifest = VariantManifest {
        path: PathBuf::from(format!("{catalogue_path}#{variant_id}")),
        name,
        tags,
        spec,
    };
    Ok(VariantManifestTask {
        manifest_path: manifest.path.display().to_string(),
        manifest,
    })
}

fn catalogue_locus(
    columns: &CatalogueColumns,
    row: &BTreeMap<String, String>,
    assembly: &str,
) -> Result<Option<GenomicLocus>, String> {
    let chrom_role = format!("coordinates.{assembly}.chrom");
    let Some(chrom) = columns
        .value(row, &chrom_role)
        .filter(|value| !value.is_empty())
    else {
        return Ok(None);
    };
    let pos_role = format!("coordinates.{assembly}.pos");
    if let Some(pos) = columns
        .value(row, &pos_role)
        .filter(|value| !value.is_empty())
    {
        let pos = parse_catalogue_i64(pos, &pos_role)?;
        return Ok(Some(GenomicLocus {
            chrom: chrom.to_owned(),
            start: pos,
            end: pos,
        }));
    }
    let start_role = format!("coordinates.{assembly}.start");
    let end_role = format!("coordinates.{assembly}.end");
    match (
        columns
            .value(row, &start_role)
            .filter(|value| !value.is_empty()),
        columns
            .value(row, &end_role)
            .filter(|value| !value.is_empty()),
    ) {
        (Some(start), Some(end)) => Ok(Some(GenomicLocus {
            chrom: chrom.to_owned(),
            start: parse_catalogue_i64(start, &start_role)?,
            end: parse_catalogue_i64(end, &end_role)?,
        })),
        (None, None) => Ok(None),
        _ => Err(format!("expected both {start_role} and {end_role}")),
    }
}

fn parse_catalogue_i64(value: &str, role: &str) -> Result<i64, String> {
    value
        .parse::<i64>()
        .map_err(|err| format!("invalid integer for {role}: {value}: {err}"))
}

fn catalogue_variant_kind(value: &str) -> VariantKind {
    match value {
        "snv" | "snp" => VariantKind::Snp,
        "deletion" => VariantKind::Deletion,
        "insertion" => VariantKind::Insertion,
        "indel" => VariantKind::Indel,
        _ => VariantKind::Other,
    }
}

fn split_list(value: Option<&str>, separator: Option<&str>) -> Vec<String> {
    let Some(value) = value else {
        return Vec::new();
    };
    let separator = separator.unwrap_or("|");
    value
        .split(separator)
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(ToOwned::to_owned)
        .collect()
}

fn parse_tsv_rows(text: &str) -> Result<Vec<BTreeMap<String, String>>, String> {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let Some(header_line) = lines.next() else {
        return Ok(Vec::new());
    };
    let headers = header_line
        .split('\t')
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    if headers.is_empty() || headers.iter().any(|header| header.trim().is_empty()) {
        return Err("TSV header contains an empty column name".to_owned());
    }
    Ok(lines
        .map(|line| {
            let cells = line.split('\t').collect::<Vec<_>>();
            headers
                .iter()
                .enumerate()
                .map(|(idx, header)| {
                    (
                        header.clone(),
                        cells.get(idx).copied().unwrap_or_default().to_owned(),
                    )
                })
                .collect()
        })
        .collect())
}

#[derive(Debug, Default)]
struct CatalogueColumns {
    role_to_column: BTreeMap<String, String>,
    role_to_separator: BTreeMap<String, String>,
}

impl CatalogueColumns {
    fn from_manifest(mapping: &serde_yaml::Mapping) -> Self {
        let mut columns = Self::default();
        if let Some(column_mappings) = mapping
            .get(serde_yaml::Value::String("columns".to_owned()))
            .and_then(serde_yaml::Value::as_mapping)
        {
            for (column, config) in column_mappings {
                let Some(column_name) = column.as_str() else {
                    continue;
                };
                let Some(config) = config.as_mapping() else {
                    continue;
                };
                let Some(role) = config
                    .get(serde_yaml::Value::String("role".to_owned()))
                    .and_then(serde_yaml::Value::as_str)
                else {
                    continue;
                };
                columns
                    .role_to_column
                    .insert(role.to_owned(), column_name.to_owned());
                if let Some(separator) = config
                    .get(serde_yaml::Value::String("list_separator".to_owned()))
                    .and_then(serde_yaml::Value::as_str)
                {
                    columns
                        .role_to_separator
                        .insert(role.to_owned(), separator.to_owned());
                }
            }
        }
        columns.insert_default("id", "variant_id");
        columns.insert_default("name", "name");
        columns.insert_default("gene", "gene");
        columns.insert_default("identifier.rsid", "rsid");
        columns.insert_default("identifier.aliases", "aliases");
        columns.insert_default("alleles.kind", "kind");
        columns.insert_default("alleles.ref", "ref");
        columns.insert_default("alleles.alts", "alts");
        columns.insert_default("alleles.observed_alts", "observed_alts");
        columns.insert_default("coordinates.grch37.chrom", "grch37_chrom");
        columns.insert_default("coordinates.grch37.pos", "grch37_pos");
        columns.insert_default("coordinates.grch37.start", "grch37_start");
        columns.insert_default("coordinates.grch37.end", "grch37_end");
        columns.insert_default("coordinates.grch38.chrom", "grch38_chrom");
        columns.insert_default("coordinates.grch38.pos", "grch38_pos");
        columns.insert_default("coordinates.grch38.start", "grch38_start");
        columns.insert_default("coordinates.grch38.end", "grch38_end");
        columns
    }

    fn insert_default(&mut self, role: &str, column: &str) {
        self.role_to_column
            .entry(role.to_owned())
            .or_insert_with(|| column.to_owned());
    }

    fn value<'a>(&self, row: &'a BTreeMap<String, String>, role: &str) -> Option<&'a str> {
        self.role_to_column
            .get(role)
            .and_then(|column| row.get(column))
            .map(String::as_str)
    }

    fn separator(&self, role: &str) -> Option<&str> {
        self.role_to_separator.get(role).map(String::as_str)
    }
}
