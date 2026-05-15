use std::path::{Component, Path, PathBuf};

use bioscript_core::RuntimeError;

use super::{BioscriptRuntime, deepest_existing_ancestor};

impl BioscriptRuntime {
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
                return Ok(path.to_path_buf());
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
        if self.virtual_file_exists(raw_path) {
            return Ok(path);
        }
        if self.uses_virtual_files() && path_is_rooted(Path::new(raw_path)) {
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

    fn virtual_file_exists(&self, raw_path: &str) -> bool {
        self.config.virtual_text_files.contains_key(raw_path)
            || self.config.virtual_binary_files.contains_key(raw_path)
            || self
                .state
                .virtual_written_text_files
                .lock()
                .expect("virtual file mutex poisoned")
                .contains_key(raw_path)
    }

    pub(super) fn virtual_key(&self, path: &Path) -> String {
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
