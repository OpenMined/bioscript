const PACKAGE_DESCRIPTOR: &str = "manifest.yaml";
const LEGACY_PACKAGE_DESCRIPTOR: &str = "bioscript-package.yaml";
const PACKAGE_CACHE_DIR: &str = ".bioscript-cache/packages";
const PACKAGE_DOWNLOAD_DIR: &str = ".bioscript-cache/downloads";
const MAX_PACKAGE_FILES: usize = 1000;
const MAX_PACKAGE_FILE_BYTES: u64 = 16 * 1024 * 1024;
const MAX_PACKAGE_TOTAL_BYTES: u64 = 64 * 1024 * 1024;

include!("package_release.rs");

fn prepare_package_entrypoint_from_arg(
    runtime_root: &Path,
    source: &Path,
) -> Result<PathBuf, String> {
    let source_text = source.to_string_lossy();
    let source_url = if is_package_url(&source_text) {
        Some(source_text.to_string())
    } else {
        None
    };
    let package_path = if let Some(url) = &source_url {
        download_package_url(runtime_root, url)?
    } else {
        source.to_path_buf()
    };
    if is_package_zip_path(&package_path) {
        let imported = import_package_zip(runtime_root, &package_path, None)?;
        Ok(imported.entrypoint)
    } else if is_package_release_path(&package_path) {
        match package_zip_from_release_manifest(runtime_root, &package_path, source_url.as_deref())?
        {
            Some(zip_path) => {
                let imported = import_package_zip(runtime_root, &zip_path, None)?;
                Ok(imported.entrypoint)
            }
            None => Ok(package_path),
        }
    } else {
        Ok(package_path)
    }
}

fn run_import_package(args: Vec<String>) -> Result<(), String> {
    let mut source: Option<PathBuf> = None;
    let mut root: Option<PathBuf> = None;
    let mut output_dir: Option<PathBuf> = None;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--root" => root = Some(PathBuf::from(next_arg(&mut iter, "--root")?)),
            "--output-dir" => {
                output_dir = Some(PathBuf::from(next_arg(&mut iter, "--output-dir")?));
            }
            value if value.starts_with('-') => return Err(format!("unexpected argument: {value}")),
            value if source.is_none() => source = Some(PathBuf::from(value)),
            value => return Err(format!("unexpected argument: {value}")),
        }
    }
    let source = source.ok_or(
        "usage: bioscript import-package <package.zip|https://.../package.zip> [--root <dir>] [--output-dir <dir>]",
    )?;
    let runtime_root = root
        .map_or_else(env::current_dir, Ok)
        .map_err(|err| format!("failed to get current directory: {err}"))?;
    let source_text = source.to_string_lossy();
    let source_url = if is_package_url(&source_text) {
        Some(source_text.to_string())
    } else {
        None
    };
    let package_path = if let Some(url) = &source_url {
        download_package_url(&runtime_root, url)?
    } else {
        absolutize(&runtime_root, &source)
    };
    let package_path =
        package_zip_from_release_manifest(&runtime_root, &package_path, source_url.as_deref())?
            .unwrap_or(package_path);
    let imported = import_package_zip(&runtime_root, &package_path, output_dir.as_deref())?;
    println!("root\t{}", imported.root.display());
    println!("entrypoint\t{}", imported.entrypoint.display());
    if let Some(name) = imported.name {
        println!("name\t{name}");
    }
    Ok(())
}

struct ImportedPackage {
    root: PathBuf,
    entrypoint: PathBuf,
    name: Option<String>,
}

fn import_package_zip(
    runtime_root: &Path,
    zip_path: &Path,
    output_dir: Option<&Path>,
) -> Result<ImportedPackage, String> {
    let replace_existing = output_dir.is_none();
    let target_root = output_dir.map_or_else(
        || package_cache_target(runtime_root, zip_path),
        |path| absolutize(runtime_root, path),
    );
    if target_root.exists() {
        if replace_existing {
            fs::remove_dir_all(&target_root).map_err(|err| {
                format!(
                    "failed to remove previous package import {}: {err}",
                    target_root.display()
                )
            })?;
        } else if target_root
            .read_dir()
            .map_err(|err| format!("failed to read output dir {}: {err}", target_root.display()))?
            .next()
            .is_some()
        {
            return Err(format!(
                "package output dir already exists and is not empty: {}",
                target_root.display()
            ));
        }
    }
    fs::create_dir_all(&target_root).map_err(|err| {
        format!(
            "failed to create package import dir {}: {err}",
            target_root.display()
        )
    })?;
    extract_package_zip(zip_path, &target_root)?;
    let descriptor = load_package_descriptor(&target_root)?;
    let entrypoint = target_root.join(&descriptor.entrypoint);
    let canonical_root = target_root.canonicalize().map_err(|err| {
        format!(
            "failed to resolve package root {}: {err}",
            target_root.display()
        )
    })?;
    let canonical_entrypoint = entrypoint.canonicalize().map_err(|err| {
        format!(
            "failed to resolve package entrypoint {}: {err}",
            entrypoint.display()
        )
    })?;
    if !canonical_entrypoint.starts_with(&canonical_root) {
        return Err(format!(
            "package entrypoint escapes package root: {}",
            descriptor.entrypoint.display()
        ));
    }
    match manifest_schema(&canonical_entrypoint)?.as_str() {
        "bioscript:panel:1.0"
        | "bioscript:assay:1.0"
        | "bioscript:variant:1.0"
        | "bioscript:variant" => {}
        other => {
            return Err(format!(
                "package entrypoint has unsupported schema '{other}'"
            ))
        }
    }
    Ok(ImportedPackage {
        root: canonical_root,
        entrypoint: canonical_entrypoint,
        name: descriptor.name,
    })
}

struct PackageDescriptor {
    entrypoint: PathBuf,
    name: Option<String>,
}

fn load_package_descriptor(root: &Path) -> Result<PackageDescriptor, String> {
    for name in [PACKAGE_DESCRIPTOR, LEGACY_PACKAGE_DESCRIPTOR] {
        let path = root.join(name);
        if path.exists() {
            let text = fs::read_to_string(&path).map_err(|err| {
                format!(
                    "failed to read package descriptor {}: {err}",
                    path.display()
                )
            })?;
            let value: serde_yaml::Value = serde_yaml::from_str(&text).map_err(|err| {
                format!(
                    "failed to parse package descriptor {}: {err}",
                    path.display()
                )
            })?;
            let schema = value
                .as_mapping()
                .and_then(|mapping| mapping.get(serde_yaml::Value::String("schema".to_owned())))
                .and_then(serde_yaml::Value::as_str)
                .ok_or_else(|| {
                    format!("package descriptor {} is missing schema", path.display())
                })?;
            if matches!(
                schema,
                "bioscript:panel:1.0"
                    | "bioscript:assay:1.0"
                    | "bioscript:variant:1.0"
                    | "bioscript:variant"
            ) {
                let package_name = value
                    .as_mapping()
                    .and_then(|mapping| mapping.get(serde_yaml::Value::String("name".to_owned())))
                    .and_then(serde_yaml::Value::as_str)
                    .map(ToOwned::to_owned);
                return Ok(PackageDescriptor {
                    entrypoint: PathBuf::from(PACKAGE_DESCRIPTOR),
                    name: package_name,
                });
            }
            if schema != "bioscript:package:1.0" {
                return Err(format!(
                    "package descriptor {} has unsupported schema '{schema}'",
                    path.display()
                ));
            }
            let entrypoint = value
                .as_mapping()
                .and_then(|mapping| mapping.get(serde_yaml::Value::String("entrypoint".to_owned())))
                .and_then(serde_yaml::Value::as_str)
                .ok_or_else(|| {
                    format!(
                        "package descriptor {} is missing entrypoint",
                        path.display()
                    )
                })?;
            let entrypoint = checked_relative_package_path(entrypoint)?;
            let name = value
                .as_mapping()
                .and_then(|mapping| mapping.get(serde_yaml::Value::String("name".to_owned())))
                .and_then(serde_yaml::Value::as_str)
                .map(ToOwned::to_owned);
            return Ok(PackageDescriptor { entrypoint, name });
        }
    }
    for candidate in ["panel.yaml", "assay.yaml", "variant.yaml"] {
        if root.join(candidate).exists() {
            return Ok(PackageDescriptor {
                entrypoint: PathBuf::from(candidate),
                name: None,
            });
        }
    }
    Err(format!(
        "package does not contain {PACKAGE_DESCRIPTOR}, {LEGACY_PACKAGE_DESCRIPTOR}, panel.yaml, assay.yaml, or variant.yaml"
    ))
}

fn extract_package_zip(zip_path: &Path, target_root: &Path) -> Result<(), String> {
    let file = fs::File::open(zip_path)
        .map_err(|err| format!("failed to open package zip {}: {err}", zip_path.display()))?;
    let mut archive = zip::ZipArchive::new(file)
        .map_err(|err| format!("failed to read package zip {}: {err}", zip_path.display()))?;
    if archive.len() > MAX_PACKAGE_FILES {
        return Err(format!(
            "package has too many entries: {} > {MAX_PACKAGE_FILES}",
            archive.len()
        ));
    }
    let mut seen = std::collections::BTreeSet::new();
    let mut total_size = 0_u64;
    for idx in 0..archive.len() {
        let mut entry = archive
            .by_index(idx)
            .map_err(|err| format!("failed to read package zip entry {idx}: {err}"))?;
        let Some(enclosed) = entry.enclosed_name() else {
            return Err(format!(
                "package zip entry has unsafe path: {}",
                entry.name()
            ));
        };
        let relative = checked_relative_package_path(&enclosed.to_string_lossy())?;
        if entry.is_dir() {
            fs::create_dir_all(target_root.join(relative)).map_err(|err| {
                format!(
                    "failed to create package directory from {}: {err}",
                    entry.name()
                )
            })?;
            continue;
        }
        if entry
            .unix_mode()
            .is_some_and(|mode| mode & 0o170_000 == 0o120_000)
        {
            return Err(format!("package zip entry is a symlink: {}", entry.name()));
        }
        if !is_allowed_package_file(&relative) {
            return Err(format!(
                "package zip entry has unsupported extension: {}",
                relative.display()
            ));
        }
        if !seen.insert(relative.clone()) {
            return Err(format!(
                "package zip contains duplicate path: {}",
                relative.display()
            ));
        }
        let size = entry.size();
        if size > MAX_PACKAGE_FILE_BYTES {
            return Err(format!(
                "package file {} exceeds {} bytes",
                relative.display(),
                MAX_PACKAGE_FILE_BYTES
            ));
        }
        total_size = total_size.saturating_add(size);
        if total_size > MAX_PACKAGE_TOTAL_BYTES {
            return Err(format!(
                "package contents exceed {MAX_PACKAGE_TOTAL_BYTES} bytes"
            ));
        }
        let target = target_root.join(&relative);
        if let Some(parent) = target.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create package dir {}: {err}", parent.display())
            })?;
        }
        let mut out = fs::File::create(&target)
            .map_err(|err| format!("failed to create package file {}: {err}", target.display()))?;
        std::io::copy(&mut entry, &mut out)
            .map_err(|err| format!("failed to extract package file {}: {err}", target.display()))?;
    }
    Ok(())
}

fn checked_relative_package_path(raw: &str) -> Result<PathBuf, String> {
    let path = Path::new(raw);
    if path.is_absolute() {
        return Err(format!("package path must be relative: {raw}"));
    }
    let mut out = PathBuf::new();
    for component in path.components() {
        match component {
            std::path::Component::Normal(part) => out.push(part),
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir
            | std::path::Component::RootDir
            | std::path::Component::Prefix(_) => {
                return Err(format!("package path escapes package root: {raw}"));
            }
        }
    }
    if out.as_os_str().is_empty() {
        return Err("package path is empty".to_owned());
    }
    Ok(out)
}

fn is_allowed_package_file(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name == PACKAGE_DESCRIPTOR || name == LEGACY_PACKAGE_DESCRIPTOR)
        || path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| {
                matches!(
                    ext.to_ascii_lowercase().as_str(),
                    "yaml" | "yml" | "py" | "md" | "txt" | "tsv" | "json" | "jsonl"
                        | "fa" | "fasta" | "fna" | "fai" | "vcf" | "tbi" | "bai"
                        | "crai"
                )
            })
}

fn package_cache_target(runtime_root: &Path, zip_path: &Path) -> PathBuf {
    let stem = zip_path
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or("package")
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' {
                ch
            } else {
                '-'
            }
        })
        .collect::<String>();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |duration| duration.as_nanos());
    runtime_root
        .join(PACKAGE_CACHE_DIR)
        .join(format!("{stem}-{nanos}"))
}

fn download_package_url(runtime_root: &Path, url: &str) -> Result<PathBuf, String> {
    if !url.starts_with("https://") {
        return Err("package URLs must use https://".to_owned());
    }
    let url_path = url.split('?').next().unwrap_or(url);
    let extension = Path::new(url_path)
        .extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if !matches!(extension.as_str(), "zip" | "yaml" | "yml") {
        return Err("package URL must point to a .zip, .yaml, or .yml file".to_owned());
    }
    let downloads = runtime_root.join(PACKAGE_DOWNLOAD_DIR);
    fs::create_dir_all(&downloads).map_err(|err| {
        format!(
            "failed to create package download dir {}: {err}",
            downloads.display()
        )
    })?;
    let file_name = url
        .split('?')
        .next()
        .unwrap_or(url)
        .rsplit('/')
        .next()
        .filter(|value| !value.is_empty())
        .unwrap_or("package.zip");
    let safe_name = file_name
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '.' | '-' | '_') {
                ch
            } else {
                '-'
            }
        })
        .collect::<String>();
    let target = downloads.join(safe_name);
    let status = std::process::Command::new("curl")
        .arg("-fL")
        .arg("--max-time")
        .arg("60")
        .arg("-o")
        .arg(&target)
        .arg(url)
        .status()
        .map_err(|err| format!("failed to run curl for package download: {err}"))?;
    if !status.success() {
        return Err(format!("package download failed for {url}"));
    }
    Ok(target)
}

fn is_package_url(value: &str) -> bool {
    value.starts_with("https://") || value.starts_with("http://")
}

fn is_package_zip_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("zip"))
}

fn is_package_release_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| matches!(ext.to_ascii_lowercase().as_str(), "yaml" | "yml"))
}

#[cfg(test)]
mod package_tests {
    use super::*;
    use std::{
        io::Write as _,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-package-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_zip(path: &Path, entries: &[(&str, &str)]) {
        let file = fs::File::create(path).unwrap();
        let mut writer = zip::ZipWriter::new(file);
        let options = zip::write::SimpleFileOptions::default();
        for (name, body) in entries {
            if name.ends_with('/') {
                writer.add_directory(*name, options).unwrap();
            } else {
                writer.start_file(*name, options).unwrap();
                writer.write_all(body.as_bytes()).unwrap();
            }
        }
        writer.finish().unwrap();
    }

    fn variant_manifest(name: &str) -> String {
        format!(
            r"
schema: bioscript:variant:1.0
name: {name}
gene: ABC
identifiers:
  rsids: [rs1]
"
        )
    }

    #[test]
    fn package_path_and_type_helpers_validate_inputs() {
        assert_eq!(
            checked_relative_package_path("./nested/manifest.yaml").unwrap(),
            PathBuf::from("nested/manifest.yaml")
        );
        assert!(checked_relative_package_path("../escape.yaml")
            .unwrap_err()
            .contains("escapes"));
        assert!(checked_relative_package_path("/absolute.yaml")
            .unwrap_err()
            .contains("must be relative"));
        assert!(checked_relative_package_path(".").unwrap_err().contains("empty"));

        assert!(is_allowed_package_file(Path::new("manifest.yaml")));
        assert!(is_allowed_package_file(Path::new("docs/readme.md")));
        assert!(is_allowed_package_file(Path::new("assets/muc1_motifs.fa")));
        assert!(is_allowed_package_file(Path::new("assets/reference.fa.fai")));
        assert!(!is_allowed_package_file(Path::new("bin/tool.sh")));
        assert!(is_package_url("https://example.test/pkg.zip"));
        assert!(is_package_url("http://example.test/pkg.zip"));
        assert!(is_package_zip_path(Path::new("PKG.ZIP")));
        assert!(is_package_release_path(Path::new("package.yml")));
        assert!(!is_package_release_path(Path::new("package.zip")));
    }

    #[test]
    fn package_descriptor_loads_supported_shapes_and_reports_errors() {
        let dir = temp_dir("descriptors");

        fs::write(dir.join("variant.yaml"), variant_manifest("Fallback")).unwrap();
        let descriptor = load_package_descriptor(&dir).unwrap();
        assert_eq!(descriptor.entrypoint, PathBuf::from("variant.yaml"));
        assert_eq!(descriptor.name, None);
        fs::remove_file(dir.join("variant.yaml")).unwrap();

        fs::write(dir.join("manifest.yaml"), variant_manifest("Direct")).unwrap();
        let descriptor = load_package_descriptor(&dir).unwrap();
        assert_eq!(descriptor.entrypoint, PathBuf::from("manifest.yaml"));
        assert_eq!(descriptor.name, Some("Direct".to_owned()));
        fs::remove_file(dir.join("manifest.yaml")).unwrap();

        fs::write(
            dir.join("manifest.yaml"),
            r"
schema: bioscript:package:1.0
name: Wrapped
entrypoint: nested/panel.yaml
",
        )
        .unwrap();
        let descriptor = load_package_descriptor(&dir).unwrap();
        assert_eq!(descriptor.entrypoint, PathBuf::from("nested/panel.yaml"));
        assert_eq!(descriptor.name, Some("Wrapped".to_owned()));

        fs::write(
            dir.join("manifest.yaml"),
            r"
schema: bioscript:package:1.0
name: Missing Entrypoint
",
        )
        .unwrap();
        assert!(descriptor_err(&dir).contains("missing entrypoint"));

        fs::write(
            dir.join("manifest.yaml"),
            r"
schema: bioscript:package:1.0
entrypoint: ../escape.yaml
",
        )
        .unwrap();
        assert!(descriptor_err(&dir).contains("escapes"));

        fs::write(dir.join("manifest.yaml"), "schema: unsupported\n").unwrap();
        assert!(descriptor_err(&dir).contains("unsupported schema"));

        fs::write(dir.join("manifest.yaml"), "name: missing schema\n").unwrap();
        assert!(descriptor_err(&dir).contains("missing schema"));

        fs::remove_file(dir.join("manifest.yaml")).unwrap();
        fs::write(
            dir.join(LEGACY_PACKAGE_DESCRIPTOR),
            r"
schema: bioscript:package:1.0
name: Legacy
entrypoint: assay.yaml
",
        )
        .unwrap();
        let descriptor = load_package_descriptor(&dir).unwrap();
        assert_eq!(descriptor.entrypoint, PathBuf::from("assay.yaml"));
        assert_eq!(descriptor.name, Some("Legacy".to_owned()));
        fs::remove_file(dir.join(LEGACY_PACKAGE_DESCRIPTOR)).unwrap();

        assert!(descriptor_err(&dir).contains("does not contain"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn import_package_zip_extracts_and_validates_entrypoint() {
        let dir = temp_dir("import");
        let zip_path = dir.join("pkg.zip");
        write_zip(
            &zip_path,
            &[
                (
                    "manifest.yaml",
                    r"
schema: bioscript:package:1.0
name: Test Package
entrypoint: nested/variant.yaml
",
                ),
                ("nested/", ""),
                ("nested/variant.yaml", &variant_manifest("Nested")),
                ("docs/readme.md", "hello"),
            ],
        );

        let imported = import_package_zip(&dir, &zip_path, Some(Path::new("out"))).unwrap();
        assert_eq!(imported.name, Some("Test Package".to_owned()));
        assert!(imported.entrypoint.ends_with("nested/variant.yaml"));
        assert!(imported.root.join("docs/readme.md").exists());

        let occupied = dir.join("occupied");
        fs::create_dir_all(&occupied).unwrap();
        fs::write(occupied.join("existing.txt"), "x").unwrap();
        assert!(import_err(&dir, &zip_path, Some(Path::new("occupied"))).contains("not empty"));

        let empty_output = dir.join("empty-output");
        fs::create_dir_all(&empty_output).unwrap();
        let imported_empty =
            import_package_zip(&dir, &zip_path, Some(Path::new("empty-output"))).unwrap();
        assert!(imported_empty.entrypoint.ends_with("nested/variant.yaml"));

        let cached = import_package_zip(&dir, &zip_path, None).unwrap();
        assert!(cached
            .root
            .starts_with(dir.join(PACKAGE_CACHE_DIR).canonicalize().unwrap()));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn extract_package_zip_rejects_bad_entries_and_bad_archives() {
        let dir = temp_dir("extract-errors");
        let bad_zip = dir.join("not.zip");
        fs::write(&bad_zip, "not a zip").unwrap();
        assert!(extract_package_zip(&bad_zip, &dir.join("out"))
            .unwrap_err()
            .contains("failed to read package zip"));

        let unsupported = dir.join("unsupported.zip");
        write_zip(&unsupported, &[("script.sh", "echo no")]);
        assert!(extract_package_zip(&unsupported, &dir.join("unsupported-out"))
            .unwrap_err()
            .contains("unsupported extension"));

        let unsafe_path = dir.join("unsafe.zip");
        write_zip(&unsafe_path, &[("../escape.yaml", "no")]);
        assert!(extract_package_zip(&unsafe_path, &dir.join("unsafe-out"))
            .unwrap_err()
            .contains("unsafe path"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn import_package_zip_reports_entrypoint_resolution_and_schema_errors() {
        let dir = temp_dir("import-entrypoint-errors");

        let missing_entrypoint = dir.join("missing-entrypoint.zip");
        write_zip(
            &missing_entrypoint,
            &[(
                "manifest.yaml",
                r"
schema: bioscript:package:1.0
entrypoint: missing.yaml
",
            )],
        );
        assert!(import_err(&dir, &missing_entrypoint, Some(Path::new("missing-out")))
            .contains("failed to resolve package entrypoint"));

        let bad_schema = dir.join("bad-schema.zip");
        write_zip(
            &bad_schema,
            &[
                (
                    "manifest.yaml",
                    r"
schema: bioscript:package:1.0
entrypoint: nested/custom.yaml
",
                ),
                ("nested/custom.yaml", "schema: custom:schema\nname: bad\n"),
            ],
        );
        assert!(import_err(&dir, &bad_schema, Some(Path::new("bad-schema-out")))
            .contains("unsupported schema"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn package_release_manifest_resolves_local_artifacts_and_hashes() {
        let dir = temp_dir("release");
        let zip_path = dir.join("pkg.zip");
        write_zip(&zip_path, &[("manifest.yaml", &variant_manifest("Release"))]);
        let digest = sha256_file(&zip_path).unwrap();
        let release_path = dir.join("package.yaml");
        fs::write(
            &release_path,
            format!(
                r"
schema: bioscript:package-release:1.0
artifact:
  path: pkg.zip
  sha256: {digest}
"
            ),
        )
        .unwrap();

        let resolved = package_zip_from_release_manifest(&dir, &release_path, None)
            .unwrap()
            .unwrap();
        assert_eq!(resolved, zip_path);

        fs::write(
            &release_path,
            r"
schema: bioscript:package-release:1.0
artifact:
  path: pkg.zip
  sha256: bad
",
        )
        .unwrap();
        assert!(package_zip_from_release_manifest(&dir, &release_path, None)
            .unwrap_err()
            .contains("sha256 mismatch"));

        fs::write(&release_path, "schema: bioscript:variant:1.0\n").unwrap();
        assert!(package_zip_from_release_manifest(&dir, &release_path, None)
            .unwrap()
            .is_none());
        assert!(package_zip_from_release_manifest(&dir, &dir.join("missing.yaml"), None)
            .unwrap()
            .is_none());

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn prepare_package_entrypoint_imports_zip_and_release_manifests() {
        let dir = temp_dir("prepare-entrypoint");
        let zip_path = dir.join("pkg.zip");
        write_zip(&zip_path, &[("manifest.yaml", &variant_manifest("Prepared"))]);

        let entrypoint = prepare_package_entrypoint_from_arg(&dir, &zip_path).unwrap();
        assert!(entrypoint.ends_with("manifest.yaml"));
        assert!(entrypoint.starts_with(dir.join(PACKAGE_CACHE_DIR).canonicalize().unwrap()));

        let digest = sha256_file(&zip_path).unwrap();
        let release_path = dir.join("package.yaml");
        fs::write(
            &release_path,
            format!(
                r"
schema: bioscript:package-release:1.0
artifact:
  path: pkg.zip
  sha256: {digest}
"
            ),
        )
        .unwrap();
        let release_entrypoint = prepare_package_entrypoint_from_arg(&dir, &release_path).unwrap();
        assert!(release_entrypoint.ends_with("manifest.yaml"));

        let plain = dir.join("panel.yaml");
        fs::write(&plain, "schema: bioscript:panel:1.0\n").unwrap();
        assert_eq!(prepare_package_entrypoint_from_arg(&dir, &plain).unwrap(), plain);

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn package_release_helpers_cover_yaml_strings_urls_and_download_validation() {
        let value: serde_yaml::Value = serde_yaml::from_str("name: Example\n").unwrap();
        assert_eq!(yaml_string(&value, "name"), Some("Example".to_owned()));
        assert_eq!(yaml_string(&value, "missing"), None);
        assert_eq!(
            join_url("https://example.test/releases/package.yaml?download=1", "pkg.zip"),
            "https://example.test/releases/pkg.zip"
        );
        assert_eq!(
            join_url("https://cdn.test/pkg.zip", "https://other.test/pkg.zip"),
            "https://other.test/pkg.zip"
        );

        let dir = temp_dir("download-validation");
        assert!(download_package_url(&dir, "http://example.test/pkg.zip")
            .unwrap_err()
            .contains("must use https"));
        assert!(download_package_url(&dir, "https://example.test/pkg.exe")
            .unwrap_err()
            .contains("must point"));
        let target = package_cache_target(&dir, Path::new("weird package!.zip"));
        assert!(target.display().to_string().contains("weird-package-"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn import_package_command_validates_arguments() {
        assert!(run_import_package(Vec::new()).unwrap_err().contains("usage"));
        assert!(run_import_package(vec!["--bad".to_owned()])
            .unwrap_err()
            .contains("unexpected argument"));
        assert!(run_import_package(vec!["one.zip".to_owned(), "two.zip".to_owned()])
            .unwrap_err()
            .contains("unexpected argument"));
        assert!(run_import_package(vec!["one.zip".to_owned(), "--root".to_owned()])
            .unwrap_err()
            .contains("--root requires"));
    }

    #[test]
    fn import_package_command_imports_local_zip_to_output_dir() {
        let dir = temp_dir("import-command");
        let zip_path = dir.join("pkg.zip");
        write_zip(
            &zip_path,
            &[
                (
                    "manifest.yaml",
                    r"
schema: bioscript:package:1.0
name: Import Command Package
entrypoint: variant.yaml
",
                ),
                ("variant.yaml", &variant_manifest("Imported")),
            ],
        );
        let output = dir.join("imported");

        run_import_package(vec![
            zip_path.display().to_string(),
            "--root".to_owned(),
            dir.display().to_string(),
            "--output-dir".to_owned(),
            output.display().to_string(),
        ])
        .unwrap();

        assert!(output.join("manifest.yaml").exists());
        assert!(output.join("variant.yaml").exists());

        fs::remove_dir_all(dir).unwrap();
    }

    fn descriptor_err(root: &Path) -> String {
        match load_package_descriptor(root) {
            Ok(_) => panic!("expected package descriptor to fail"),
            Err(err) => err,
        }
    }

    fn import_err(runtime_root: &Path, zip_path: &Path, output_dir: Option<&Path>) -> String {
        match import_package_zip(runtime_root, zip_path, output_dir) {
            Ok(_) => panic!("expected package import to fail"),
            Err(err) => err,
        }
    }
}
