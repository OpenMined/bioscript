use std::collections::BTreeMap;

use bioscript_core::{Assembly, VariantObservation};

pub const MANIFEST_ROW_TSV_HEADERS: [&str; 13] = [
    "kind",
    "name",
    "path",
    "tags",
    "participant_id",
    "backend",
    "matched_rsid",
    "assembly",
    "genotype",
    "ref_count",
    "alt_count",
    "depth",
    "evidence",
];

pub fn variant_row(
    path: &str,
    name: &str,
    tags: &[String],
    observation: &VariantObservation,
    participant_id: &str,
) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("kind".to_owned(), "variant".to_owned());
    row.insert("name".to_owned(), name.to_owned());
    row.insert("path".to_owned(), path.to_owned());
    row.insert("tags".to_owned(), tags.join(","));
    row.insert("backend".to_owned(), observation.backend.clone());
    row.insert("participant_id".to_owned(), participant_id.to_owned());
    row.insert(
        "matched_rsid".to_owned(),
        observation.matched_rsid.clone().unwrap_or_default(),
    );
    row.insert(
        "assembly".to_owned(),
        observation
            .assembly
            .map(assembly_row_value)
            .unwrap_or_default(),
    );
    row.insert(
        "genotype".to_owned(),
        observation.genotype.clone().unwrap_or_default(),
    );
    row.insert(
        "ref_count".to_owned(),
        observation
            .ref_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "alt_count".to_owned(),
        observation
            .alt_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "depth".to_owned(),
        observation
            .depth
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "raw_counts".to_owned(),
        serde_json::to_string(&observation.raw_counts).unwrap_or_default(),
    );
    row.insert("evidence".to_owned(), observation.evidence.join(" | "));
    row
}

pub fn render_manifest_rows_tsv(rows: &[BTreeMap<String, String>]) -> String {
    let mut out = MANIFEST_ROW_TSV_HEADERS.join("\t");
    out.push('\n');
    for row in rows {
        let line = MANIFEST_ROW_TSV_HEADERS
            .iter()
            .map(|header| {
                row.get(*header)
                    .cloned()
                    .unwrap_or_default()
                    .replace('\t', " ")
            })
            .collect::<Vec<_>>()
            .join("\t");
        out.push_str(&line);
        out.push('\n');
    }
    out
}

pub fn render_manifest_trace_tsv(rows: &[BTreeMap<String, String>]) -> String {
    let mut trace = String::from("step\tline\tcode\n");
    for (idx, row) in rows.iter().enumerate() {
        trace.push_str(&(idx + 1).to_string());
        trace.push('\t');
        trace.push_str(&(idx + 1).to_string());
        trace.push('\t');
        trace.push_str(&row.get("path").cloned().unwrap_or_default());
        trace.push('\n');
    }
    trace
}

fn assembly_row_value(assembly: Assembly) -> String {
    match assembly {
        Assembly::Grch37 => "grch37".to_owned(),
        Assembly::Grch38 => "grch38".to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::{MANIFEST_ROW_TSV_HEADERS, render_manifest_rows_tsv, render_manifest_trace_tsv};

    #[test]
    fn renders_manifest_rows_with_cli_header_order() {
        let mut row = std::collections::BTreeMap::new();
        row.insert("kind".to_owned(), "variant".to_owned());
        row.insert("name".to_owned(), "APOE\trs429358".to_owned());
        row.insert("path".to_owned(), "assets/APOE/rs429358.yaml".to_owned());
        row.insert("genotype".to_owned(), "TT".to_owned());

        let text = render_manifest_rows_tsv(&[row]);
        let mut lines = text.lines();
        assert_eq!(lines.next().unwrap(), MANIFEST_ROW_TSV_HEADERS.join("\t"));
        assert_eq!(
            lines.next().unwrap(),
            "variant\tAPOE rs429358\tassets/APOE/rs429358.yaml\t\t\t\t\t\tTT\t\t\t\t"
        );
    }

    #[test]
    fn renders_manifest_trace_with_cli_format() {
        let mut first = std::collections::BTreeMap::new();
        first.insert("path".to_owned(), "assets/APOE/rs429358.yaml".to_owned());
        let mut second = std::collections::BTreeMap::new();
        second.insert("path".to_owned(), "assets/APOE/rs7412.yaml".to_owned());

        assert_eq!(
            render_manifest_trace_tsv(&[first, second]),
            "step\tline\tcode\n1\t1\tassets/APOE/rs429358.yaml\n2\t2\tassets/APOE/rs7412.yaml\n"
        );
    }
}
