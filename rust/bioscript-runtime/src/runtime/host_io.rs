use std::{fs, io::Read, path::Path};

use bioscript_core::RuntimeError;
use monty::MontyObject;

use super::{BioscriptRuntime, args::expect_string_arg, args::reject_kwargs};

const MAX_HOST_TEXT_BYTES: u64 = 16 * 1024 * 1024;

pub(crate) fn host_read_text(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "read_text")?;
    let path = runtime.resolve_existing_user_path(&expect_string_arg(args, 0, "read_text")?)?;
    let content = read_text_limited(&path, MAX_HOST_TEXT_BYTES)?;
    Ok(MontyObject::String(content))
}

pub(crate) fn host_write_text(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "write_text")?;
    let path = runtime.resolve_user_write_path(&expect_string_arg(args, 0, "write_text")?)?;
    let content = expect_string_arg(args, 1, "write_text")?;
    if u64::try_from(content.len()).unwrap_or(u64::MAX) > MAX_HOST_TEXT_BYTES {
        return Err(RuntimeError::InvalidArguments(format!(
            "write_text content exceeds {MAX_HOST_TEXT_BYTES} bytes"
        )));
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to create parent dir {}: {err}",
                parent.display()
            ))
        })?;
    }
    fs::write(&path, content)
        .map_err(|err| RuntimeError::Io(format!("failed to write {}: {err}", path.display())))?;
    Ok(MontyObject::None)
}

pub(crate) fn deepest_existing_ancestor(path: &Path) -> &Path {
    let mut current = path;
    while !current.exists() {
        let Some(parent) = current.parent() else {
            break;
        };
        current = parent;
    }
    current
}

fn read_text_limited(path: &Path, max_bytes: u64) -> Result<String, RuntimeError> {
    let mut file = fs::File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to read {}: {err}", path.display())))?;
    let mut bytes = Vec::new();
    file.by_ref()
        .take(max_bytes.saturating_add(1))
        .read_to_end(&mut bytes)
        .map_err(|err| RuntimeError::Io(format!("failed to read {}: {err}", path.display())))?;
    if u64::try_from(bytes.len()).unwrap_or(u64::MAX) > max_bytes {
        return Err(RuntimeError::InvalidArguments(format!(
            "read_text input {} exceeds {} bytes",
            path.display(),
            max_bytes
        )));
    }
    String::from_utf8(bytes).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to decode {} as UTF-8: {err}",
            path.display()
        ))
    })
}
