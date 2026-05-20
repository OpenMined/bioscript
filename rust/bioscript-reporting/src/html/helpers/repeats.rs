use std::{cmp::Ordering, fmt::Write as _};

use super::html_escape;

pub(crate) fn repeat_notation(value: &str) -> String {
    if let Some((unit, repeat_count)) = repeat_unit_count(value) {
        return format!("{unit}[{repeat_count}]");
    }
    value.to_owned()
}

pub(crate) fn padded_homopolymer_repeat_notation(value: &str) -> String {
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

pub(crate) fn repeat_notation_html(value: &str, reference: Option<&str>) -> String {
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

pub(crate) fn paired_repeat_notation(reference: &str, alternate: &str) -> String {
    if let Some(display) = paired_repeat_display(reference, alternate) {
        return display.alternate;
    }
    repeat_notation(alternate)
}

pub(crate) fn paired_repeat_notation_html(reference: &str, alternate: &str) -> String {
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

pub(crate) fn paired_reference_repeat_notation(reference: &str, alternate: &str) -> String {
    paired_repeat_display(reference, alternate)
        .map_or_else(|| repeat_notation(reference), |display| display.reference)
}

pub(crate) fn paired_reference_repeat_notation_html(reference: &str, alternate: &str) -> String {
    paired_repeat_display(reference, alternate).map_or_else(
        || repeat_notation_html(reference, None),
        |display| html_escape(&display.reference),
    )
}

pub(crate) fn padded_homopolymer_repeat_notation_html(value: &str, reference: &str) -> String {
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

pub(crate) fn genotype_repeat_notation(value: &str, ref_allele: &str, alt_allele: &str) -> String {
    if value.is_empty() {
        return String::new();
    }
    if value.contains('/') || value.contains('|') {
        return value
            .split_inclusive(['/', '|'])
            .fold(String::new(), |mut formatted, token| {
                let separator = token
                    .chars()
                    .last()
                    .filter(|ch| matches!(ch, '/' | '|'))
                    .map(|ch| ch.to_string())
                    .unwrap_or_default();
                let allele = token.trim_end_matches(['/', '|']);
                let _ = write!(
                    formatted,
                    "{}{}",
                    genotype_allele_repeat_notation(allele, ref_allele, alt_allele),
                    separator
                );
                formatted
            });
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

pub(crate) fn genotype_repeat_notation_html(
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
            .fold(String::new(), |mut formatted, token| {
                let separator = token
                    .chars()
                    .last()
                    .filter(|ch| matches!(ch, '/' | '|'))
                    .map(|ch| ch.to_string())
                    .unwrap_or_default();
                let allele = token.trim_end_matches(['/', '|']);
                let _ = write!(
                    formatted,
                    "{}{}",
                    genotype_allele_repeat_notation_html(allele, ref_allele, alt_allele),
                    html_escape(&separator)
                );
                formatted
            });
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
                .map(char::to_ascii_uppercase)
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
    match count.cmp(&reference_count) {
        Ordering::Less => "repeat-less",
        Ordering::Equal => "repeat-same",
        Ordering::Greater => "repeat-more",
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

pub(crate) fn same_homopolymer_base(reference: &str, alternate: &str) -> Option<char> {
    let reference_base = homopolymer_base(reference)?;
    let alternate_base = homopolymer_base(alternate)?;
    if reference_base == alternate_base {
        Some(reference_base)
    } else {
        None
    }
}

fn homopolymer_base(value: &str) -> Option<char> {
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
    repeat_count_delta_class(count, reference_count)
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
    repeat_count_delta_class(count, reference_count)
}
