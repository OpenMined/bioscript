use bioscript_core::RuntimeError;

pub(crate) fn rewrite_bioscript_imports(code: &str) -> Result<String, RuntimeError> {
    let mut out = Vec::new();
    for line in code.lines() {
        out.push(rewrite_line(line)?);
    }
    let rewritten = out.join("\n");
    if code.ends_with('\n') {
        Ok(rewritten + "\n")
    } else {
        Ok(rewritten)
    }
}

fn rewrite_line(line: &str) -> Result<String, RuntimeError> {
    let trimmed = line.trim_start();
    let Some(rest) = trimmed.strip_prefix("from bioscript import ") else {
        return Ok(line.to_owned());
    };
    if rest.contains(',') {
        return Err(RuntimeError::InvalidArguments(
            "BioScript currently supports one library import per line".to_owned(),
        ));
    }

    let indent_len = line.len() - trimmed.len();
    let indent = &line[..indent_len];
    let parts: Vec<&str> = rest.split_whitespace().collect();
    let (module, binding) = match parts.as_slice() {
        [module] => (*module, *module),
        [module, "as", binding] => (*module, *binding),
        _ => {
            return Err(RuntimeError::InvalidArguments(format!(
                "unsupported BioScript import syntax: {line}"
            )));
        }
    };
    validate_identifier(module, "module")?;
    validate_identifier(binding, "binding")?;
    Ok(format!(
        "{indent}{binding} = __bioscript_import__(\"{module}\")"
    ))
}

fn validate_identifier(value: &str, label: &str) -> Result<(), RuntimeError> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Err(RuntimeError::InvalidArguments(format!(
            "BioScript import {label} cannot be empty"
        )));
    };
    if !(first == '_' || first.is_ascii_alphabetic()) {
        return Err(RuntimeError::InvalidArguments(format!(
            "BioScript import {label} {value:?} is not a valid identifier"
        )));
    }
    if chars.all(|ch| ch == '_' || ch.is_ascii_alphanumeric()) {
        Ok(())
    } else {
        Err(RuntimeError::InvalidArguments(format!(
            "BioScript import {label} {value:?} is not a valid identifier"
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::rewrite_bioscript_imports;

    #[test]
    fn rewrites_single_bioscript_library_imports() {
        assert_eq!(
            rewrite_bioscript_imports("from bioscript import pysam\n").unwrap(),
            "pysam = __bioscript_import__(\"pysam\")\n"
        );
        assert_eq!(
            rewrite_bioscript_imports("    from bioscript import pyfaidx as fa\n").unwrap(),
            "    fa = __bioscript_import__(\"pyfaidx\")\n"
        );
    }

    #[test]
    fn rejects_multi_import_for_now() {
        let err = rewrite_bioscript_imports("from bioscript import pysam, pyfaidx").unwrap_err();
        assert!(err.to_string().contains("one library import per line"));
    }
}
