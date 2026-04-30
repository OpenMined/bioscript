use std::{
    fs,
    path::{Path, PathBuf},
};

use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
};

use noodles::{cram, fasta};

use crate::genotype::GenotypeSourceFormat;

#[derive(Debug, Clone, Default)]
pub struct PrepareRequest {
    pub root: PathBuf,
    pub cwd: PathBuf,
    pub cache_dir: PathBuf,
    pub input_file: Option<String>,
    pub input_format: Option<GenotypeSourceFormat>,
    pub reference_file: Option<String>,
}

#[derive(Debug, Clone, Default)]
pub struct PreparedPaths {
    pub input_file: Option<PathBuf>,
    pub input_index: Option<PathBuf>,
    pub reference_file: Option<PathBuf>,
    pub reference_index: Option<PathBuf>,
    pub cache_dir: PathBuf,
}

pub fn prepare_indexes(request: &PrepareRequest) -> Result<PreparedPaths, String> {
    let root = canonical_dir(&request.root)?;
    let cache_dir = resolve_cache_dir(&request.cwd, &request.cache_dir);
    fs::create_dir_all(&cache_dir)
        .map_err(|err| format!("failed to create cache dir {}: {err}", cache_dir.display()))?;

    let input_file = request
        .input_file
        .as_deref()
        .map(|value| resolve_rooted_path(&root, value))
        .transpose()?;
    let reference_file = request
        .reference_file
        .as_deref()
        .map(|value| resolve_rooted_path(&root, value))
        .transpose()?;

    let input_index = match (&input_file, request.input_format) {
        (Some(path), Some(GenotypeSourceFormat::Cram)) => {
            Some(ensure_alignment_index(path, &cache_dir)?)
        }
        (Some(path), None) if detect_alignment_input(path) => {
            Some(ensure_alignment_index(path, &cache_dir)?)
        }
        _ => None,
    };

    let (prepared_reference_file, reference_index) = match reference_file {
        Some(path) => {
            let (reference_file, reference_index) = ensure_reference_index(&path, &cache_dir)?;
            (Some(reference_file), Some(reference_index))
        }
        None => (None, None),
    };

    Ok(PreparedPaths {
        input_file,
        input_index,
        reference_file: prepared_reference_file,
        reference_index,
        cache_dir,
    })
}

fn canonical_dir(path: &Path) -> Result<PathBuf, String> {
    path.canonicalize()
        .map_err(|err| format!("failed to canonicalize {}: {err}", path.display()))
}

fn resolve_rooted_path(root: &Path, raw: &str) -> Result<PathBuf, String> {
    let raw_path = Path::new(raw);
    let resolved = if raw_path.is_absolute() {
        raw_path
            .canonicalize()
            .map_err(|err| format!("failed to resolve {}: {err}", raw_path.display()))?
    } else {
        ensure_relative_path_safe(raw_path)?;
        root.join(raw_path)
    };
    let canonical = resolved
        .canonicalize()
        .map_err(|err| format!("failed to resolve {}: {err}", resolved.display()))?;
    ensure_path_within_root(root, &canonical)?;
    Ok(canonical)
}

fn ensure_relative_path_safe(path: &Path) -> Result<(), String> {
    for component in path.components() {
        match component {
            std::path::Component::ParentDir
            | std::path::Component::RootDir
            | std::path::Component::Prefix(_) => {
                return Err(format!("path escapes bioscript root: {}", path.display()));
            }
            std::path::Component::CurDir | std::path::Component::Normal(_) => {}
        }
    }
    Ok(())
}

fn ensure_path_within_root(root: &Path, path: &Path) -> Result<(), String> {
    let relative = path
        .strip_prefix(root)
        .map_err(|_| format!("path escapes bioscript root: {}", path.display()))?;
    ensure_relative_path_safe(relative)
}

fn resolve_cache_dir(cwd: &Path, cache_dir: &Path) -> PathBuf {
    if cache_dir.is_absolute() {
        cache_dir.to_path_buf()
    } else {
        cwd.join(cache_dir)
    }
}

fn detect_alignment_input(path: &Path) -> bool {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    lower.ends_with(".cram") || lower.ends_with(".bam")
}

fn ensure_alignment_index(path: &Path, cache_dir: &Path) -> Result<PathBuf, String> {
    if let Some(existing) = adjacent_alignment_index(path) {
        return Ok(existing);
    }

    let is_cram = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("cram"));
    if !is_cram {
        return Err(format!(
            "alignment indexing only supports CRAM in the pure-Rust backend: {}",
            path.display()
        ));
    }

    let ext = "crai";
    let out = cache_dir.join(format!("{}.{ext}", stable_stem(path)));
    if out.exists() {
        return Ok(out);
    }

    let index = cram::fs::index(path).map_err(|err| {
        format!(
            "failed to build alignment index {} for {}: {err}",
            out.display(),
            path.display()
        )
    })?;
    cram::crai::fs::write(&out, &index).map_err(|err| {
        format!(
            "failed to write alignment index {} for {}: {err}",
            out.display(),
            path.display()
        )
    })?;
    Ok(out)
}

fn adjacent_alignment_index(path: &Path) -> Option<PathBuf> {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    let candidates = if lower.ends_with(".cram") {
        vec![
            path.with_extension("cram.crai"),
            path.with_extension("crai"),
        ]
    } else {
        Vec::new()
    };

    candidates.into_iter().find(|candidate| candidate.exists())
}

fn ensure_reference_index(path: &Path, cache_dir: &Path) -> Result<(PathBuf, PathBuf), String> {
    let adjacent = adjacent_reference_index(path);
    if let Some(index) = adjacent {
        return Ok((path.to_path_buf(), index));
    }

    let cached_reference = cache_dir.join(cache_reference_name(path));
    if !cached_reference.exists() {
        create_reference_link(path, &cached_reference)?;
    }

    let cached_index = adjacent_reference_index(&cached_reference)
        .unwrap_or_else(|| cached_reference_index_path(&cached_reference));
    if !cached_index.exists() {
        let index = fasta::fs::index(&cached_reference).map_err(|err| {
            format!(
                "failed to build FASTA index {} for {}: {err}",
                cached_index.display(),
                cached_reference.display()
            )
        })?;
        fasta::fai::fs::write(&cached_index, &index).map_err(|err| {
            format!(
                "failed to write FASTA index {} for {}: {err}",
                cached_index.display(),
                cached_reference.display()
            )
        })?;
    }

    Ok((cached_reference, cached_index))
}

fn adjacent_reference_index(path: &Path) -> Option<PathBuf> {
    let candidate = cached_reference_index_path(path);
    candidate.exists().then_some(candidate)
}

fn cached_reference_index_path(path: &Path) -> PathBuf {
    if let Some(ext) = path.extension().and_then(|ext| ext.to_str()) {
        path.with_extension(format!("{ext}.fai"))
    } else {
        path.with_extension("fai")
    }
}

fn create_reference_link(source: &Path, target: &Path) -> Result<(), String> {
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent)
            .map_err(|err| format!("failed to create cache dir {}: {err}", parent.display()))?;
    }

    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(source, target).map_err(|err| {
            format!(
                "failed to create cached reference link {} -> {}: {err}",
                target.display(),
                source.display()
            )
        })?;
        Ok(())
    }

    #[cfg(not(unix))]
    {
        fs::copy(source, target).map_err(|err| {
            format!(
                "failed to copy cached reference {} -> {}: {err}",
                source.display(),
                target.display()
            )
        })?;
        Ok(())
    }
}

fn stable_stem(path: &Path) -> String {
    let mut hasher = DefaultHasher::new();
    path.to_string_lossy().hash(&mut hasher);
    let hash = hasher.finish();
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("input")
        .replace(['/', ' ', ':'], "_");
    format!("{file_name}-{hash:016x}")
}
fn cache_reference_name(path: &Path) -> String {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("reference.fa")
        .replace(['/', ' ', ':'], "_");
    format!("{}-ref", stable_stem(Path::new(&file_name)))
}

pub fn shell_flags(prepared: &PreparedPaths) -> String {
    let mut parts = Vec::new();
    if let Some(path) = &prepared.input_file {
        parts.push(format!("--input-file {}", shell_quote(path)));
    }
    if let Some(path) = &prepared.input_index {
        parts.push(format!("--input-index {}", shell_quote(path)));
    }
    if let Some(path) = &prepared.reference_file {
        parts.push(format!("--reference-file {}", shell_quote(path)));
    }
    if let Some(path) = &prepared.reference_index {
        parts.push(format!("--reference-index {}", shell_quote(path)));
    }
    parts.join(" ")
}

fn shell_quote(path: &Path) -> String {
    let text = path.to_string_lossy();
    format!("'{}'", text.replace('\'', "'\"'\"'"))
}
