use super::helpers::{
    class_cell, html_escape, json_field_as_tsv, render_table_start, table_column_class, value_str,
};
use std::fmt::Write as _;
pub(super) fn render_observation_table(
    out: &mut String,
    observations: &[serde_json::Value],
    show_participant_id: bool,
) {
    render_observation_filters(out);
    let all_headers = [
        "participant_id",
        "outcome",
        "rsid",
        "gene",
        "ref_alt",
        "genotype_display",
        "genotype",
        "zygosity",
        "assembly",
        "chrom",
        "pos_start",
        "pos_end",
        "kind",
        "ref_count",
        "alt_count",
        "depth",
        "genotype_quality",
        "allele_balance",
        "evidence_type",
        "evidence_raw",
        "match_quality",
        "match_notes",
        "facets",
        "assay_id",
        "assay_version",
        "variant_key",
        "match_status",
        "coverage_status",
        "call_status",
        "source",
    ];
    let show_counts = observations.iter().any(observation_has_quantitative_depth);
    let show_genotype_quality = observations
        .iter()
        .any(|observation| !json_field_as_tsv(observation.get("genotype_quality")).is_empty());
    let show_facets = observations
        .iter()
        .any(|observation| !json_field_as_tsv(observation.get("facets")).is_empty());
    let show_match_quality = observations.iter().any(|observation| {
        !json_field_as_tsv(observation.get("match_quality")).is_empty()
            || !json_field_as_tsv(observation.get("match_notes")).is_empty()
    });
    let show_imputed_reference_note = observations
        .iter()
        .any(observation_is_imputed_vcf_reference);
    let show_weak_indel_note = observations.iter().any(observation_is_weak_indel_match);
    let headers = all_headers
        .iter()
        .copied()
        .filter(|header| show_participant_id || *header != "participant_id")
        .filter(|header| {
            show_counts
                || !matches!(
                    *header,
                    "ref_count" | "alt_count" | "depth" | "allele_balance"
                )
        })
        .filter(|header| show_genotype_quality || *header != "genotype_quality")
        .filter(|header| show_match_quality || !matches!(*header, "match_quality" | "match_notes"))
        .filter(|header| show_facets || *header != "facets")
        .collect::<Vec<_>>();
    render_table_start(out, "observations-table", &headers);
    for observation in observations {
        let _ = write!(
            out,
            "<tr class=\"{}\" data-observation=\"{}\" data-participant=\"{}\">",
            observation_row_class(observation),
            observation_filter_group(observation),
            html_escape(value_str(observation, "participant_id"))
        );
        for header in &headers {
            render_observation_cell(out, observation, header);
        }
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></div>");
    if show_imputed_reference_note {
        out.push_str("<p class=\"muted observation-note\">* In variant-only VCF inputs, absent queried variant rows are shown as imputed reference genotypes. This is usually appropriate for variant-only VCFs, but it may be wrong if the VCF omits loci for another reason.</p>");
    }
    if show_weak_indel_note {
        out.push_str("<p class=\"muted observation-note\">* Indel calls from consumer genotype files are weak matches: the file reports an insertion/deletion token at the marker, but does not provide sequence-resolved evidence for the exact deletion allele.</p>");
    }
}

pub(super) fn observation_filter_group(observation: &serde_json::Value) -> &'static str {
    match observation_row_class(observation) {
        "row-reference" => "reference",
        "row-missing" => "missing",
        _ => "variant",
    }
}

pub(super) fn render_observation_filters(out: &mut String) {
    out.push_str("<div class=\"level-filter\"><span class=\"muted\">Observations:</span>");
    for (outcome, label) in [
        ("variant", "Show variants"),
        ("reference", "Show reference"),
        ("missing", "Show missing"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-observation-filter=\"{outcome}\" checked onchange=\"applyTableFilters('observations-table')\"> {label}</label>"
        );
    }
    out.push_str("<button class=\"filter-action\" type=\"button\" onclick=\"setObservationFilterGroup(true)\">Show all</button><button class=\"filter-action\" type=\"button\" onclick=\"setObservationFilterGroup(false)\">Hide all</button>");
    out.push_str("<label><input type=\"checkbox\" onchange=\"toggleDebug(this.checked)\"> Show debug</label>");
    out.push_str("</div>");
}

pub(super) fn observation_has_quantitative_depth(observation: &serde_json::Value) -> bool {
    ["ref_count", "alt_count", "depth", "allele_balance"]
        .iter()
        .any(|key| !json_field_as_tsv(observation.get(*key)).is_empty())
}

pub(super) fn observation_row_class(observation: &serde_json::Value) -> &'static str {
    let outcome = observation
        .get("outcome")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if outcome == "variant" {
        "row-variant"
    } else if outcome == "reference" {
        "row-reference"
    } else if observation
        .get("call_status")
        .and_then(serde_json::Value::as_str)
        != Some("called")
        || observation
            .get("match_status")
            .and_then(serde_json::Value::as_str)
            == Some("not_found")
    {
        "row-missing"
    } else {
        ""
    }
}

pub(super) fn render_observation_cell(
    out: &mut String,
    observation: &serde_json::Value,
    header: &str,
) {
    let cell_class = table_column_class(header);
    if header == "outcome" {
        let mut value = json_field_as_tsv(observation.get(header));
        if (value == "reference" && observation_is_imputed_vcf_reference(observation))
            || (value == "variant" && observation_is_weak_indel_match(observation))
        {
            value.push('*');
        }
        let _ = write!(
            out,
            "<td class=\"{}\">{}</td>",
            cell_class,
            html_escape(&value)
        );
        return;
    }
    if header == "ref_alt" {
        ref_alt_cell(out, observation);
        return;
    }
    if header == "allele_balance" {
        let value = observation
            .get(header)
            .and_then(serde_json::Value::as_f64)
            .map(|value| format!("{value:.2}"))
            .unwrap_or_default();
        class_cell(out, &value, cell_class);
        return;
    }
    if header == "source" {
        let source = observation
            .get("source")
            .unwrap_or(&serde_json::Value::Null);
        let url = source
            .get("url")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        if url.is_empty() {
            let _ = write!(out, "<td class=\"{cell_class}\"></td>");
        } else {
            let _ = write!(
                out,
                "<td class=\"{}\"><a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">Source</a></td>",
                cell_class,
                html_escape(url)
            );
        }
        return;
    }
    if header == "genotype_display" {
        let outcome = observation
            .get("outcome")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let value = json_field_as_tsv(observation.get(header));
        if matches!(outcome, "variant" | "observed_alt" | "unknown_alt") {
            let alt = observation
                .get("alt")
                .and_then(serde_json::Value::as_str)
                .unwrap_or_default();
            let _ = write!(
                out,
                "<td class=\"genotype-hit {}\">{}</td>",
                cell_class,
                highlight_allele(&value, alt)
            );
            return;
        }
    }
    let _ = write!(
        out,
        "<td class=\"{}\">{}</td>",
        cell_class,
        html_escape(&json_field_as_tsv(observation.get(header)))
    );
}

fn ref_alt_cell(out: &mut String, observation: &serde_json::Value) {
    let value = observation_ref_alt(observation);
    let escaped_value = html_escape(&value);
    let _ = write!(
        out,
        "<td class=\"mono ref-alt-cell\" title=\"{escaped_value}\"><span class=\"truncate-cell\">{escaped_value}</span></td>"
    );
}

pub(super) fn observation_is_imputed_vcf_reference(observation: &serde_json::Value) -> bool {
    observation
        .get("evidence_raw")
        .and_then(serde_json::Value::as_str)
        .is_some_and(|evidence| {
            evidence.contains("imputed reference genotype from absent variant-only VCF record")
        })
}

pub(super) fn observation_is_weak_indel_match(observation: &serde_json::Value) -> bool {
    observation
        .get("match_quality")
        .and_then(serde_json::Value::as_str)
        == Some("weak")
}

pub(super) fn observation_ref_alt(observation: &serde_json::Value) -> String {
    let ref_allele = observation
        .get("ref")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let alt_allele = observation
        .get("alt")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if ref_allele.is_empty() && alt_allele.is_empty() {
        String::new()
    } else {
        format!("{ref_allele}->{}", alt_allele.replace(',', "/"))
    }
}

pub(super) fn highlight_allele(value: &str, allele: &str) -> String {
    if value.is_empty() || allele.is_empty() {
        return html_escape(value);
    }
    if allele.chars().count() == 1 {
        let target = allele
            .chars()
            .next()
            .unwrap_or_default()
            .to_ascii_uppercase();
        let mut out = String::new();
        for ch in value.chars() {
            let escaped = html_escape(&ch.to_string());
            if ch.to_ascii_uppercase() == target {
                let _ = write!(out, "<span class=\"allele-hit\">{escaped}</span>");
            } else {
                out.push_str(&escaped);
            }
        }
        return out;
    }
    let escaped_value = html_escape(value);
    let escaped_allele = html_escape(allele);
    escaped_value.replace(
        &escaped_allele,
        &format!("<span class=\"allele-hit\">{escaped_allele}</span>"),
    )
}
