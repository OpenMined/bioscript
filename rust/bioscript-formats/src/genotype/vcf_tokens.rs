use super::normalize_genotype;

pub(crate) fn genotype_from_vcf_gt(
    gt: &str,
    reference: &str,
    alternates: &[&str],
) -> Option<String> {
    if matches!(gt.trim(), "" | "." | "./." | ".|.") {
        return Some("--".to_owned());
    }

    let cleaned = gt.trim().replace('|', "/");
    let parts: Vec<&str> = cleaned.split('/').collect();
    if parts.len() != 2 || parts.contains(&".") {
        return Some("--".to_owned());
    }

    let ref_token = vcf_reference_token(reference, alternates);
    let mut out = String::new();
    for part in parts {
        let idx = part.parse::<usize>().ok()?;
        if idx == 0 {
            out.push_str(&ref_token);
        } else {
            let alt = alternates.get(idx - 1)?;
            out.push_str(&vcf_alt_token(reference, alt));
        }
    }

    Some(normalize_genotype(&out))
}

pub(crate) fn vcf_reference_token(reference: &str, alternates: &[&str]) -> String {
    let mut saw_shorter = false;
    let mut saw_longer = false;

    for alt in alternates {
        if is_symbolic_vcf_alt(alt) {
            continue;
        }
        match alt.len().cmp(&reference.len()) {
            std::cmp::Ordering::Less => saw_shorter = true,
            std::cmp::Ordering::Greater => saw_longer = true,
            std::cmp::Ordering::Equal => {}
        }
    }

    match (saw_shorter, saw_longer) {
        (true, false) => "I".to_owned(),
        (false, true) => "D".to_owned(),
        _ => normalize_sequence_token(reference),
    }
}

pub(crate) fn vcf_alt_token(reference: &str, alternate: &str) -> String {
    if is_symbolic_vcf_alt(alternate) {
        return "--".to_owned();
    }
    match alternate.len().cmp(&reference.len()) {
        std::cmp::Ordering::Less => "D".to_owned(),
        std::cmp::Ordering::Greater => "I".to_owned(),
        std::cmp::Ordering::Equal => normalize_sequence_token(alternate),
    }
}

pub(crate) fn is_symbolic_vcf_alt(alternate: &str) -> bool {
    let trimmed = alternate.trim();
    trimmed.starts_with('<') && trimmed.ends_with('>')
}

pub(crate) fn normalize_sequence_token(value: &str) -> String {
    value.trim().to_ascii_uppercase()
}
