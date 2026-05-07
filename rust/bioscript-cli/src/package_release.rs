fn package_zip_from_release_manifest(
    runtime_root: &Path,
    path: &Path,
    source_url: Option<&str>,
) -> Result<Option<PathBuf>, String> {
    if !is_package_release_path(path) || !path.exists() {
        return Ok(None);
    }
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read package release {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse package release {}: {err}", path.display()))?;
    let schema = yaml_string(&value, "schema");
    if schema.as_deref() != Some("bioscript:package-release:1.0") {
        return Ok(None);
    }
    let artifact = value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("artifact".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .ok_or_else(|| format!("package release {} is missing artifact", path.display()))?;
    let artifact_path = artifact
        .get(serde_yaml::Value::String("path".to_owned()))
        .and_then(serde_yaml::Value::as_str);
    let artifact_url = artifact
        .get(serde_yaml::Value::String("url".to_owned()))
        .and_then(serde_yaml::Value::as_str);
    let zip_path = if let Some(url) = artifact_url {
        download_package_url(runtime_root, url)?
    } else if let Some(relative) = artifact_path {
        if let Some(base_url) = source_url {
            download_package_url(runtime_root, &join_url(base_url, relative))?
        } else {
            path.parent()
                .ok_or_else(|| format!("package release has no parent: {}", path.display()))?
                .join(checked_relative_package_path(relative)?)
        }
    } else {
        return Err(format!(
            "package release {} artifact needs path or url",
            path.display()
        ));
    };
    if let Some(expected) = artifact
        .get(serde_yaml::Value::String("sha256".to_owned()))
        .and_then(serde_yaml::Value::as_str)
    {
        let actual = sha256_file(&zip_path)?;
        if actual != expected {
            return Err(format!(
                "package artifact sha256 mismatch for {}: expected {expected}, got {actual}",
                zip_path.display()
            ));
        }
    }
    Ok(Some(zip_path))
}

fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String(key.to_owned())))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

fn sha256_file(path: &Path) -> Result<String, String> {
    use sha2::{Digest, Sha256};

    let mut file = fs::File::open(path)
        .map_err(|err| format!("failed to open artifact {}: {err}", path.display()))?;
    let mut digest = Sha256::new();
    let mut buffer = vec![0_u8; 1024 * 64];
    loop {
        let n = std::io::Read::read(&mut file, &mut buffer)
            .map_err(|err| format!("failed to read artifact {}: {err}", path.display()))?;
        if n == 0 {
            break;
        }
        digest.update(&buffer[..n]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn join_url(base_url: &str, relative: &str) -> String {
    if relative.starts_with("https://") || relative.starts_with("http://") {
        return relative.to_owned();
    }
    let base = base_url.split('?').next().unwrap_or(base_url);
    match base.rsplit_once('/') {
        Some((prefix, _)) => format!("{prefix}/{relative}"),
        None => relative.to_owned(),
    }
}
