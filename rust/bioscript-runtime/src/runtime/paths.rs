use std::fs;
use std::path::{Component, Path, PathBuf};

use bioscript_core::RuntimeError;

use super::{BioscriptRuntime, deepest_existing_ancestor};

impl BioscriptRuntime {
    /// Real on-disk temp directory mirroring the virtual filesystem. Native
    /// tool facades (samtools/kestrel/bcftools) can only operate on real
    /// files, so when the analysis runs under a virtual filesystem we mirror
    /// every virtual path `/X` to `<materialized_root>/X` and write virtual
    /// content there on first access. Created lazily on first resolve.
    pub(super) fn materialized_root(&self) -> Result<PathBuf, RuntimeError> {
        let mut guard = self
            .state
            .materialized_root
            .lock()
            .expect("materialized_root mutex poisoned");
        if let Some(dir) = guard.as_ref() {
            return Ok(dir.clone());
        }
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let dir =
            std::env::temp_dir().join(format!("bioscript-vfs-{}-{nanos}", std::process::id()));
        fs::create_dir_all(&dir).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to create materialized vfs root {}: {err}",
                dir.display()
            ))
        })?;
        *guard = Some(dir.clone());
        Ok(dir)
    }

    fn materialized_root_if_set(&self) -> Option<PathBuf> {
        self.state
            .materialized_root
            .lock()
            .expect("materialized_root mutex poisoned")
            .clone()
    }

    /// Canonical virtual key for a raw virtual path (e.g. `/input/genotypes`).
    fn canonical_virtual_key(raw_path: &str) -> String {
        let normalized = raw_path.replace('\\', "/");
        if normalized.starts_with('/') {
            normalized
        } else {
            format!("/{normalized}")
        }
    }

    /// Write the virtual content backing `raw_path` (script-provided config
    /// files or text written earlier in the run) to its mirrored real path so
    /// native tools can read it. No-op if the real file already exists (e.g.
    /// a prior native tool produced it).
    fn materialize_virtual_content(
        &self,
        raw_path: &str,
        real_path: &Path,
    ) -> Result<(), RuntimeError> {
        if real_path.exists() {
            return Ok(());
        }
        let key = Self::canonical_virtual_key(raw_path);
        if let Some(parent) = real_path.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to create materialized dir {}: {err}",
                    parent.display()
                ))
            })?;
        }
        if let Some(bytes) = self.config.virtual_binary_files.get(&key) {
            fs::write(real_path, bytes).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to materialize {}: {err}",
                    real_path.display()
                ))
            })?;
            return Ok(());
        }
        if let Some(text) = self.config.virtual_text_files.get(&key) {
            fs::write(real_path, text).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to materialize {}: {err}",
                    real_path.display()
                ))
            })?;
            return Ok(());
        }
        let written = self
            .state
            .virtual_written_text_files
            .lock()
            .expect("virtual file mutex poisoned");
        if let Some(text) = written.get(&key) {
            fs::write(real_path, text).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to materialize {}: {err}",
                    real_path.display()
                ))
            })?;
        }
        Ok(())
    }

    pub(super) fn resolve_user_path(&self, raw_path: &str) -> Result<PathBuf, RuntimeError> {
        let path = Path::new(raw_path);
        if path_is_rooted(path) {
            if self.uses_virtual_files() {
                for component in path.components() {
                    match component {
                        Component::ParentDir | Component::Prefix(_) => {
                            return Err(RuntimeError::InvalidArguments(format!(
                                "path escapes bioscript root: {raw_path}"
                            )));
                        }
                        Component::RootDir | Component::CurDir | Component::Normal(_) => {}
                    }
                }
                // Mirror the virtual path into the real materialized root so
                // native tool facades receive a real on-disk path.
                let rel = raw_path.trim_start_matches('/');
                return Ok(self.materialized_root()?.join(rel));
            }
            return Err(RuntimeError::InvalidArguments(format!(
                "absolute paths are not allowed: {raw_path}"
            )));
        }
        for component in path.components() {
            match component {
                Component::ParentDir | Component::RootDir | Component::Prefix(_) => {
                    return Err(RuntimeError::InvalidArguments(format!(
                        "path escapes bioscript root: {raw_path}"
                    )));
                }
                Component::CurDir | Component::Normal(_) => {}
            }
        }
        Ok(self.root.join(path))
    }

    pub(super) fn resolve_existing_user_path(
        &self,
        raw_path: &str,
    ) -> Result<PathBuf, RuntimeError> {
        let path = self.resolve_user_path(raw_path)?;
        if self.uses_virtual_files() {
            // `path` is the mirrored real path. Write any backing virtual
            // content (script-provided config files, or text the script
            // wrote earlier) to disk so native tools can read it. Files a
            // prior native tool already produced are left untouched.
            self.materialize_virtual_content(raw_path, &path)?;
            if path.exists() {
                return Ok(path);
            }
            return Err(RuntimeError::Io(format!(
                "virtual file does not exist: {raw_path}"
            )));
        }
        let canonical = path.canonicalize().map_err(|err| {
            RuntimeError::Io(format!("failed to resolve {}: {err}", path.display()))
        })?;
        self.ensure_under_root(&canonical, raw_path)?;
        Ok(canonical)
    }

    pub(super) fn resolve_user_write_path(&self, raw_path: &str) -> Result<PathBuf, RuntimeError> {
        let path = self.resolve_user_path(raw_path)?;
        if self.uses_virtual_files() {
            if path_is_rooted(Path::new(raw_path))
                && !(raw_path.starts_with("/output/") || raw_path.starts_with("/work/"))
            {
                return Err(RuntimeError::InvalidArguments(format!(
                    "virtual write path must be under /work or /output: {raw_path}"
                )));
            }
            // `path` is the mirrored real path under the materialized root;
            // make sure its parent exists so native tools can write there.
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent).map_err(|err| {
                    RuntimeError::Io(format!(
                        "failed to create materialized dir {}: {err}",
                        parent.display()
                    ))
                })?;
            }
            return Ok(path);
        }
        if path.exists() {
            let canonical = path.canonicalize().map_err(|err| {
                RuntimeError::Io(format!("failed to resolve {}: {err}", path.display()))
            })?;
            self.ensure_under_root(&canonical, raw_path)?;
            return Ok(canonical);
        }

        let parent = path.parent().unwrap_or(&self.root);
        let existing_parent = deepest_existing_ancestor(parent);
        let canonical_parent = existing_parent.canonicalize().map_err(|err| {
            RuntimeError::Io(format!(
                "failed to resolve parent dir {}: {err}",
                existing_parent.display()
            ))
        })?;
        self.ensure_under_root(&canonical_parent, raw_path)?;
        Ok(path)
    }

    pub(super) fn uses_virtual_files(&self) -> bool {
        !self.config.virtual_text_files.is_empty() || !self.config.virtual_binary_files.is_empty()
    }

    pub(super) fn virtual_key(&self, path: &Path) -> String {
        // A mirrored real path (under the materialized root) maps back to its
        // canonical virtual key `/X` so the in-memory virtual text store stays
        // consistent with the script-provided config keys and the report
        // runner's `/output/...` lookups.
        if let Some(mat) = self.materialized_root_if_set()
            && let Ok(rel) = path.strip_prefix(&mat)
        {
            return format!("/{}", rel.display()).replace('\\', "/");
        }
        path.strip_prefix(&self.root)
            .unwrap_or(path)
            .display()
            .to_string()
            .replace('\\', "/")
    }

    fn ensure_under_root(&self, path: &Path, raw_path: &str) -> Result<(), RuntimeError> {
        if path.starts_with(&self.root) {
            Ok(())
        } else {
            Err(RuntimeError::InvalidArguments(format!(
                "path escapes bioscript root: {raw_path}"
            )))
        }
    }
}

// Treat any path with a leading root (`/...`) as rooted/absolute. We cannot use
// `Path::is_absolute()` here: on `wasm32-unknown-unknown` (the wasm-pack build
// target) `is_absolute()` is gated to `unix`/`wasi` and always returns `false`
// for `/`-rooted paths, which would route the virtual report paths
// (`/input/...`, `/work/...`, `/output/...`) into the relative-path branch and
// reject their `RootDir` component as escaping the bioscript root.
// `Path::has_root()` parses components and is target-independent.
fn path_is_rooted(path: &Path) -> bool {
    path.has_root()
}

pub(crate) fn resolve_optional_loader_path(
    runtime: &BioscriptRuntime,
    path: Option<PathBuf>,
) -> Result<Option<PathBuf>, RuntimeError> {
    path.map(|path| {
        if path_is_rooted(&path) {
            Ok(path)
        } else {
            runtime.resolve_user_path(&path.to_string_lossy())
        }
    })
    .transpose()
}
