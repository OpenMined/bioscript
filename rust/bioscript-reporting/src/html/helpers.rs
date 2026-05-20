use std::fmt::Write as _;
pub(super) fn render_table_start(out: &mut String, table_id: &str, headers: &[&str]) {
    let escaped_id = html_escape(table_id);
    let _ = write!(
        out,
        "<div class=\"table-tools\"><input type=\"search\" placeholder=\"Filter table\" data-filter-for=\"{escaped_id}\" oninput=\"scheduleTableFilter('{escaped_id}')\"></div><div class=\"table-wrap\"><table id=\"{escaped_id}\"><thead><tr>"
    );
    for (index, header) in headers.iter().enumerate() {
        let _ = write!(
            out,
            "<th class=\"{}\" onclick=\"sortTable('{}',{})\">{}<span class=\"sort-mark\"></span></th>",
            table_column_class(header),
            escaped_id,
            index,
            html_escape(&table_header_label(header))
        );
    }
    out.push_str("</tr></thead><tbody>");
}

pub(super) fn table_column_class(header: &str) -> &'static str {
    if is_debug_column(header) {
        "debug-col"
    } else {
        ""
    }
}

pub(super) fn is_debug_column(header: &str) -> bool {
    matches!(
        header,
        "allele_balance"
            | "genotype_quality"
            | "evidence_type"
            | "evidence_raw"
            | "match_quality"
            | "match_notes"
            | "assay_id"
            | "assay_version"
            | "variant_key"
            | "match_status"
            | "coverage_status"
            | "call_status"
    )
}

pub(super) fn table_header_label(header: &str) -> String {
    match header {
        "participant_id" => "Participant ID".to_owned(),
        "rsid" => "RSID".to_owned(),
        "gene" => "Gene".to_owned(),
        "ref" => "Ref".to_owned(),
        "alt" => "Alt".to_owned(),
        "ref_alt" | "Ref/Alt" => "Ref / Alt".to_owned(),
        "genotype_display" => "Genotype".to_owned(),
        "genotype" => "GT".to_owned(),
        "zygosity" => "Zygosity".to_owned(),
        "outcome" => "Outcome".to_owned(),
        "match_status" => "Match Status".to_owned(),
        "coverage_status" => "Coverage Status".to_owned(),
        "call_status" => "Call Status".to_owned(),
        "assembly" => "Assembly".to_owned(),
        "chrom" => "Chrom".to_owned(),
        "pos_start" => "Start".to_owned(),
        "pos_end" => "End".to_owned(),
        "kind" => "Kind".to_owned(),
        "ref_count" => "Ref Count".to_owned(),
        "alt_count" => "Alt Count".to_owned(),
        "depth" => "Depth".to_owned(),
        "genotype_quality" => "Genotype Quality".to_owned(),
        "allele_balance" => "Allele Balance".to_owned(),
        "evidence_type" => "Evidence Type".to_owned(),
        "source" => "Source".to_owned(),
        "evidence_raw" => "Evidence Raw".to_owned(),
        "match_quality" => "Match Quality".to_owned(),
        "match_notes" => "Match Notes".to_owned(),
        "facets" => "Facets".to_owned(),
        "assay_id" => "Assay ID".to_owned(),
        "assay_version" => "Assay Version".to_owned(),
        "variant_key" => "Variant Key".to_owned(),
        other => other.to_owned(),
    }
}

pub(super) fn render_table_end(out: &mut String) {
    out.push_str("</tbody></table></div>");
}

pub(super) fn table_cell(out: &mut String, value: &str) {
    class_cell(out, value, "");
}

pub(super) fn class_cell(out: &mut String, value: &str, class_name: &str) {
    if class_name.is_empty() {
        let _ = write!(out, "<td>{}</td>", html_escape(value));
    } else {
        let _ = write!(
            out,
            "<td class=\"{}\">{}</td>",
            class_name,
            html_escape(value)
        );
    }
}

pub(super) fn link_cell(out: &mut String, url: &str) {
    if url.is_empty() {
        out.push_str("<td></td>");
    } else {
        let escaped = html_escape(url);
        let _ = write!(
            out,
            "<td><a href=\"{escaped}\" target=\"_blank\" rel=\"noopener noreferrer\">source</a></td>"
        );
    }
}

pub(super) fn value_str<'a>(value: &'a serde_json::Value, key: &str) -> &'a str {
    value
        .get(key)
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
}

pub(super) fn join_string_array(value: Option<&serde_json::Value>) -> String {
    value
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_json::Value::as_str)
                .collect::<Vec<_>>()
                .join(", ")
        })
        .unwrap_or_default()
}

pub(super) fn join_drugs(finding: &serde_json::Value) -> String {
    finding
        .get("drugs")
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(|drug| drug.get("name").and_then(serde_json::Value::as_str))
                .collect::<Vec<_>>()
                .join(", ")
        })
        .unwrap_or_default()
}

pub(super) fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}

pub(super) fn repeat_notation(value: &str) -> String {
    if let Some((unit, repeat_count)) = repeat_unit_count(value) {
        return format!("{unit}[{repeat_count}]");
    }
    value.to_owned()
}

pub(super) fn padded_homopolymer_repeat_notation(value: &str) -> String {
    if let Some((unit, repeat_count)) = padded_homopolymer_repeat_unit_count(value) {
        return format!("{unit}[{repeat_count}]");
    }
    repeat_notation(value)
}

#[cfg(test)]
pub(super) fn allele_list_repeat_notation(value: &str) -> String {
    value
        .split(',')
        .map(repeat_notation)
        .collect::<Vec<_>>()
        .join(", ")
}

pub(super) fn repeat_notation_html(value: &str, reference: Option<&str>) -> String {
    let notation = repeat_notation(value);
    let class_name = repeat_delta_class(value, reference);
    if class_name.is_empty() || class_name == "repeat-same" {
        return html_escape(&notation);
    }
    format!(
        "<span class=\"{class_name}\">{}</span>",
        html_escape(&notation)
    )
}

pub(super) fn paired_repeat_notation(reference: &str, alternate: &str) -> String {
    if let Some(display) = paired_repeat_display(reference, alternate) {
        return display.alternate;
    }
    repeat_notation(alternate)
}

pub(super) fn paired_repeat_notation_html(reference: &str, alternate: &str) -> String {
    if let Some(display) = paired_repeat_display(reference, alternate) {
        if display.alt_class.is_empty() || display.alt_class == "repeat-same" {
            return html_escape(&display.alternate);
        }
        return format!(
            "<span class=\"{}\">{}</span>",
            display.alt_class,
            html_escape(&display.alternate)
        );
    }
    repeat_notation_html(alternate, Some(reference))
}

pub(super) fn paired_reference_repeat_notation(reference: &str, alternate: &str) -> String {
    paired_repeat_display(reference, alternate)
        .map(|display| display.reference)
        .unwrap_or_else(|| repeat_notation(reference))
}

pub(super) fn paired_reference_repeat_notation_html(reference: &str, alternate: &str) -> String {
    paired_repeat_display(reference, alternate)
        .map(|display| html_escape(&display.reference))
        .unwrap_or_else(|| repeat_notation_html(reference, None))
}

pub(super) fn padded_homopolymer_repeat_notation_html(value: &str, reference: &str) -> String {
    let notation = padded_homopolymer_repeat_notation(value);
    let class_name = padded_homopolymer_repeat_delta_class(value, reference);
    if class_name.is_empty() || class_name == "repeat-same" {
        return html_escape(&notation);
    }
    format!(
        "<span class=\"{class_name}\">{}</span>",
        html_escape(&notation)
    )
}

struct PairedRepeatDisplay {
    reference: String,
    alternate: String,
    alt_class: &'static str,
}

fn paired_repeat_display(reference: &str, alternate: &str) -> Option<PairedRepeatDisplay> {
    if reference.is_empty() || !reference.chars().all(is_sequence_char) {
        return None;
    }
    if let Some((unit, ref_count)) = exact_repeat_unit_count(reference, 2)
        && let Some(alt_count) = repeat_count_for_unit(alternate, &unit)
    {
        return Some(PairedRepeatDisplay {
            reference: repeat_token(&unit, ref_count),
            alternate: repeat_token(&unit, alt_count),
            alt_class: repeat_count_delta_class(alt_count, ref_count),
        });
    }
    anchored_suffix_repeat_display(reference, alternate)
}

fn anchored_suffix_repeat_display(reference: &str, alternate: &str) -> Option<PairedRepeatDisplay> {
    if alternate.is_empty() || !reference.ends_with(alternate) || reference.len() <= alternate.len()
    {
        return None;
    }
    let repeat_part = &reference[..reference.len() - alternate.len()];
    let (unit, ref_count) = exact_repeat_unit_count(repeat_part, 2)?;
    Some(PairedRepeatDisplay {
        reference: format!("{}{}", repeat_token(&unit, ref_count), alternate),
        alternate: alternate.to_owned(),
        alt_class: "repeat-less",
    })
}

pub(super) fn genotype_repeat_notation(value: &str, ref_allele: &str, alt_allele: &str) -> String {
    if value.is_empty() {
        return String::new();
    }
    if value.contains('/') || value.contains('|') {
        return value
            .split_inclusive(['/', '|'])
            .map(|token| {
                let separator = token
                    .chars()
                    .last()
                    .filter(|ch| matches!(ch, '/' | '|'))
                    .map(|ch| ch.to_string())
                    .unwrap_or_default();
                let allele = token.trim_end_matches(['/', '|']);
                format!(
                    "{}{}",
                    genotype_allele_repeat_notation(allele, ref_allele, alt_allele),
                    separator
                )
            })
            .collect::<String>();
    }
    if ref_allele.is_empty()
        || alt_allele.is_empty()
        || (ref_allele.len() == 1 && alt_allele.len() == 1)
    {
        return repeat_notation(value);
    }
    for (left, right) in [
        (ref_allele, ref_allele),
        (alt_allele, alt_allele),
        (ref_allele, alt_allele),
        (alt_allele, ref_allele),
    ] {
        if value.eq_ignore_ascii_case(&format!("{left}{right}")) {
            return format!("{}/{}", repeat_notation(left), repeat_notation(right));
        }
    }
    repeat_notation(value)
}

pub(super) fn genotype_repeat_notation_html(
    value: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> String {
    if value.is_empty() {
        return String::new();
    }
    if value.contains('/') || value.contains('|') {
        return value
            .split_inclusive(['/', '|'])
            .map(|token| {
                let separator = token
                    .chars()
                    .last()
                    .filter(|ch| matches!(ch, '/' | '|'))
                    .map(|ch| ch.to_string())
                    .unwrap_or_default();
                let allele = token.trim_end_matches(['/', '|']);
                format!(
                    "{}{}",
                    genotype_allele_repeat_notation_html(allele, ref_allele, alt_allele),
                    html_escape(&separator)
                )
            })
            .collect::<String>();
    }
    if ref_allele.is_empty()
        || alt_allele.is_empty()
        || (ref_allele.len() == 1 && alt_allele.len() == 1)
    {
        return repeat_notation_html(value, Some(ref_allele));
    }
    for (left, right) in [
        (ref_allele, ref_allele),
        (alt_allele, alt_allele),
        (ref_allele, alt_allele),
        (alt_allele, ref_allele),
    ] {
        if value.eq_ignore_ascii_case(&format!("{left}{right}")) {
            return format!(
                "{}/{}",
                genotype_allele_repeat_notation_html(left, ref_allele, alt_allele),
                genotype_allele_repeat_notation_html(right, ref_allele, alt_allele)
            );
        }
    }
    repeat_notation_html(value, Some(ref_allele))
}

fn genotype_allele_repeat_notation(allele: &str, ref_allele: &str, alt_allele: &str) -> String {
    if same_homopolymer_base(ref_allele, alt_allele).is_some()
        && padded_homopolymer_repeat_unit_count(allele).is_some()
    {
        return padded_homopolymer_repeat_notation(allele);
    }
    if allele.eq_ignore_ascii_case(ref_allele)
        && let Some(display) = paired_repeat_display(ref_allele, alt_allele)
    {
        return display.reference;
    }
    if let Some(display) = paired_repeat_display(ref_allele, allele) {
        return display.alternate;
    }
    repeat_notation(allele)
}

fn genotype_allele_repeat_notation_html(
    allele: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> String {
    if same_homopolymer_base(ref_allele, alt_allele).is_some()
        && padded_homopolymer_repeat_unit_count(allele).is_some()
    {
        return padded_homopolymer_repeat_notation_html(allele, ref_allele);
    }
    if allele.eq_ignore_ascii_case(ref_allele)
        && let Some(display) = paired_repeat_display(ref_allele, alt_allele)
    {
        return html_escape(&display.reference);
    }
    if paired_repeat_display(ref_allele, allele).is_some() {
        return paired_repeat_notation_html(ref_allele, allele);
    }
    repeat_notation_html(allele, Some(ref_allele))
}

fn is_sequence_char(ch: char) -> bool {
    matches!(ch.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N')
}

fn repeat_unit_count(value: &str) -> Option<(String, usize)> {
    if !value.chars().all(is_sequence_char) {
        return None;
    }
    exact_repeat_unit_count(value, 3)
}

fn exact_repeat_unit_count(value: &str, min_count: usize) -> Option<(String, usize)> {
    if value.is_empty() || !value.chars().all(is_sequence_char) {
        return None;
    }
    let chars = value.chars().collect::<Vec<_>>();
    for unit_len in 1..=chars.len().min(32) {
        if chars.len() % unit_len != 0 {
            continue;
        }
        let repeat_count = chars.len() / unit_len;
        if repeat_count < min_count {
            continue;
        }
        if chars
            .iter()
            .enumerate()
            .all(|(idx, ch)| ch.eq_ignore_ascii_case(&chars[idx % unit_len]))
        {
            let unit = chars[..unit_len]
                .iter()
                .map(|ch| ch.to_ascii_uppercase())
                .collect::<String>();
            return Some((unit, repeat_count));
        }
    }
    None
}

fn repeat_count_for_unit(value: &str, unit: &str) -> Option<usize> {
    if value.is_empty() {
        return Some(0);
    }
    if unit.is_empty() || !value.chars().all(is_sequence_char) {
        return None;
    }
    let value_chars = value.chars().collect::<Vec<_>>();
    let unit_chars = unit.chars().collect::<Vec<_>>();
    if value_chars.len() % unit_chars.len() != 0 {
        return None;
    }
    let repeat_count = value_chars.len() / unit_chars.len();
    if value_chars
        .iter()
        .enumerate()
        .all(|(idx, ch)| ch.eq_ignore_ascii_case(&unit_chars[idx % unit_chars.len()]))
    {
        Some(repeat_count)
    } else {
        None
    }
}

fn repeat_token(unit: &str, count: usize) -> String {
    let unit = unit.to_ascii_uppercase();
    if unit.chars().count() == 1 {
        format!("{unit}[{count}]")
    } else {
        format!("({unit})[{count}]")
    }
}

fn repeat_count_delta_class(count: usize, reference_count: usize) -> &'static str {
    if count < reference_count {
        "repeat-less"
    } else if count > reference_count {
        "repeat-more"
    } else {
        "repeat-same"
    }
}

fn padded_homopolymer_repeat_unit_count(value: &str) -> Option<(String, usize)> {
    let chars = value.chars().collect::<Vec<_>>();
    if !(3..=6).contains(&chars.len()) || !chars.iter().all(|ch| is_sequence_char(*ch)) {
        return None;
    }
    let first = chars.first()?.to_ascii_uppercase();
    if chars.iter().all(|ch| ch.to_ascii_uppercase() == first) {
        return Some((first.to_string(), chars.len() - 1));
    }
    None
}

pub(super) fn same_homopolymer_base(reference: &str, alternate: &str) -> Option<char> {
    let reference_base = homopolymer_base(reference)?;
    let alternate_base = homopolymer_base(alternate)?;
    if reference_base == alternate_base {
        Some(reference_base)
    } else {
        None
    }
}

pub(super) fn homopolymer_base(value: &str) -> Option<char> {
    let mut chars = value.chars();
    let first = chars.next()?.to_ascii_uppercase();
    if !is_sequence_char(first) {
        return None;
    }
    if chars.all(|ch| ch.to_ascii_uppercase() == first) {
        Some(first)
    } else {
        None
    }
}

fn repeat_delta_class(value: &str, reference: Option<&str>) -> &'static str {
    let Some(reference) = reference else {
        return "";
    };
    let Some((unit, count)) = repeat_unit_count(value) else {
        return "";
    };
    let Some((reference_unit, reference_count)) = repeat_unit_count(reference) else {
        return "";
    };
    if unit != reference_unit {
        return "";
    }
    if count < reference_count {
        "repeat-less"
    } else if count > reference_count {
        "repeat-more"
    } else {
        "repeat-same"
    }
}

fn padded_homopolymer_repeat_delta_class(value: &str, reference: &str) -> &'static str {
    let Some((unit, count)) = padded_homopolymer_repeat_unit_count(value) else {
        return "";
    };
    let Some((reference_unit, reference_count)) = padded_homopolymer_repeat_unit_count(reference)
    else {
        return "";
    };
    if unit != reference_unit {
        return "";
    }
    if count < reference_count {
        "repeat-less"
    } else if count > reference_count {
        "repeat-more"
    } else {
        "repeat-same"
    }
}

pub(super) fn html_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use super::{
        allele_list_repeat_notation, genotype_repeat_notation, genotype_repeat_notation_html,
        padded_homopolymer_repeat_notation, paired_reference_repeat_notation,
        paired_repeat_notation, paired_repeat_notation_html, repeat_notation, repeat_notation_html,
    };

    #[test]
    fn formats_homopolymer_and_tandem_repeats() {
        assert_eq!(repeat_notation("TTTTTTTTTTTTT"), "T[13]");
        assert_eq!(repeat_notation("CACACACACACACACACACACACA"), "CA[12]");
        assert_eq!(repeat_notation(&"CAG".repeat(36)), "CAG[36]");
    }

    #[test]
    fn leaves_short_or_non_repetitive_sequences_unchanged() {
        assert_eq!(repeat_notation("ATGCATGA"), "ATGCATGA");
        assert_eq!(repeat_notation("TT"), "TT");
    }

    #[test]
    fn formats_multiallelic_and_concatenated_genotype_display() {
        assert_eq!(
            allele_list_repeat_notation("TTTTTTTTTTTT,TTTTTTTTTTTTTT"),
            "T[12], T[14]"
        );
        assert_eq!(
            genotype_repeat_notation("TTTTTTTTTTTTTTTTTTTTTTTTT", "TTTTTTTTTTTTT", "TTTTTTTTTTTT"),
            "T[13]/T[12]"
        );
    }

    #[test]
    fn repeat_notation_html_marks_repeat_delta_from_reference() {
        assert_eq!(
            repeat_notation_html("TTTTTTTTTTTT", Some("TTTTTTTTTTTTT")),
            r#"<span class="repeat-less">T[12]</span>"#
        );
        assert_eq!(
            repeat_notation_html("TTTTTTTTTTTTTT", Some("TTTTTTTTTTTTT")),
            r#"<span class="repeat-more">T[14]</span>"#
        );
    }

    #[test]
    fn padded_homopolymer_notation_omits_shared_anchor_base() {
        assert_eq!(padded_homopolymer_repeat_notation("AAAA"), "A[3]");
        assert_eq!(padded_homopolymer_repeat_notation("AAA"), "A[2]");
        assert_eq!(
            genotype_repeat_notation("AAA/AAA", "AAAA", "AAA"),
            "A[2]/A[2]"
        );
    }

    #[test]
    fn genotype_repeat_notation_html_marks_repeat_delta_from_reference() {
        assert_eq!(
            genotype_repeat_notation_html(
                "TTTTTTTTTTTTT/TTTTTTTTTTTT",
                "TTTTTTTTTTTT",
                "TTTTTTTTTTTTT"
            ),
            r#"<span class="repeat-more">T[13]</span>/T[12]"#
        );
        assert_eq!(
            genotype_repeat_notation_html("AAA/AAAA", "AAAA", "AAA"),
            r#"<span class="repeat-less">A[2]</span>/A[3]"#
        );
    }

    #[test]
    fn paired_repeat_notation_detects_shared_motifs_and_anchored_deletions() {
        assert_eq!(
            paired_reference_repeat_notation("ATTCTGTCATTCTGTC", "ATTCTGTC"),
            "(ATTCTGTC)[2]"
        );
        assert_eq!(
            paired_repeat_notation("ATTCTGTCATTCTGTC", "ATTCTGTC"),
            "(ATTCTGTC)[1]"
        );
        assert_eq!(
            paired_repeat_notation_html("ATTCTGTCATTCTGTC", "ATTCTGTC"),
            r#"<span class="repeat-less">(ATTCTGTC)[1]</span>"#
        );
        assert_eq!(
            paired_reference_repeat_notation("CACACACACACACACAC", "C"),
            "(CA)[8]C"
        );
        assert_eq!(paired_repeat_notation("CACACACACACACACAC", "C"), "C");
        assert_eq!(
            genotype_repeat_notation("ATTCTGTC/ATTCTGTC", "ATTCTGTCATTCTGTC", "ATTCTGTC"),
            "(ATTCTGTC)[1]/(ATTCTGTC)[1]"
        );
        assert_eq!(
            genotype_repeat_notation(
                "CACACACACACACACAC/CACACACACACACACAC",
                "CACACACACACACACAC",
                "C"
            ),
            "(CA)[8]C/(CA)[8]C"
        );
        assert_eq!(
            genotype_repeat_notation_html(
                "CACACACACACACACAC/CACACACACACACACAC",
                "CACACACACACACACAC",
                "C"
            ),
            "(CA)[8]C/(CA)[8]C"
        );
    }
}
