fn validate_url_string(
    value: &str,
    path: &str,
    require_origin_only: bool,
    issues: &mut Vec<Issue>,
) {
    let normalized = if require_origin_only {
        normalize_origin(value)
    } else {
        normalize_download_url(value)
    };
    if let Err(message) = normalized {
        issues.push(Issue {
            severity: Severity::Error,
            path: path.to_owned(),
            message,
        });
    }
}

fn normalize_origin(value: &str) -> Result<String, String> {
    let url = Url::parse(value).map_err(|err| format!("invalid URL: {err}"))?;
    if !matches!(url.scheme(), "http" | "https") {
        return Err("expected http or https origin".to_owned());
    }
    if url.host_str().is_none() {
        return Err("origin is missing host".to_owned());
    }
    if url.path() != "/" || url.query().is_some() || url.fragment().is_some() {
        return Err("expected origin only, without path, query, or fragment".to_owned());
    }
    let mut origin = format!("{}://{}", url.scheme(), url.host_str().unwrap_or_default());
    if let Some(port) = url.port() {
        let _ = write!(origin, ":{port}");
    }
    Ok(origin)
}

fn normalize_download_url(value: &str) -> Result<String, String> {
    let url = Url::parse(value).map_err(|err| format!("invalid URL: {err}"))?;
    if !matches!(url.scheme(), "http" | "https") {
        return Err("expected http or https URL".to_owned());
    }
    if url.host_str().is_none() {
        return Err("URL is missing host".to_owned());
    }
    let mut origin = format!("{}://{}", url.scheme(), url.host_str().unwrap_or_default());
    if let Some(port) = url.port() {
        let _ = write!(origin, ":{port}");
    }
    Ok(origin)
}

fn is_allowed_chromosome(value: &str) -> bool {
    matches!(value, "X" | "Y" | "MT")
        || value
            .parse::<u8>()
            .is_ok_and(|chrom| (1..=22).contains(&chrom))
}

fn is_base_allele(value: &str) -> bool {
    matches!(value, "A" | "C" | "G" | "T")
}

fn is_rsid(value: &str) -> bool {
    value.starts_with("rs") && value[2..].chars().all(|ch| ch.is_ascii_digit())
}

fn is_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .chars()
            .all(|ch| ch.is_ascii_hexdigit() && !ch.is_ascii_uppercase())
}

fn i64_at_mapping(mapping: &Mapping, key: &str) -> Option<i64> {
    mapping
        .get(Value::String(key.to_owned()))
        .and_then(Value::as_i64)
}

fn required_non_empty_string(root: &Value, path: &[&str]) -> Result<String, String> {
    scalar_at(root, path)
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| format!("{} missing or empty", path.join(".")))
}

fn render_single_manifest_errors(path: &Path, issues: &[Issue]) -> String {
    let mut out = format!("invalid manifest {}:\n", path.display());
    for issue in issues {
        let _ = writeln!(
            out,
            "  - [{}] {}: {}",
            issue.severity, issue.path, issue.message
        );
    }
    out
}

fn load_yaml(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))
}

fn require_const(root: &Value, path: &[&str], expected: &str, issues: &mut Vec<Issue>) {
    match scalar_at(root, path) {
        Some(actual) if actual == expected => {}
        Some(actual) => issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: format!("expected '{expected}', found '{actual}'"),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: "missing required field".to_owned(),
        }),
    }
}

fn require_path(root: &Value, path: &[&str], issues: &mut Vec<Issue>) {
    if value_at(root, path).is_none() {
        issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: "missing required field".to_owned(),
        });
    }
}

fn value_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut current = root;
    for key in path {
        let mapping = current.as_mapping()?;
        current = mapping.get(Value::String((*key).to_owned()))?;
    }
    Some(current)
}

fn mapping_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Mapping> {
    value_at(root, path)?.as_mapping()
}

fn scalar_at(root: &Value, path: &[&str]) -> Option<String> {
    value_at(root, path).and_then(|value| match value {
        Value::String(text) => Some(text.clone()),
        Value::Number(number) => Some(number.to_string()),
        _ => None,
    })
}

fn seq_of_strings(root: &Value, path: &[&str]) -> Option<Vec<String>> {
    value_at(root, path)?.as_sequence().map(|items| {
        items
            .iter()
            .filter_map(|item| item.as_str().map(ToOwned::to_owned))
            .collect()
    })
}
