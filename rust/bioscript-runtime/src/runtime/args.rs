use std::collections::BTreeMap;

use bioscript_core::RuntimeError;
use monty::MontyObject;

pub(crate) fn reject_kwargs(
    kwargs: &[(MontyObject, MontyObject)],
    function_name: &str,
) -> Result<(), RuntimeError> {
    if kwargs.is_empty() {
        Ok(())
    } else {
        Err(RuntimeError::InvalidArguments(format!(
            "{function_name} does not accept keyword arguments"
        )))
    }
}

pub(crate) fn expect_string_arg(
    args: &[MontyObject],
    index: usize,
    function_name: &str,
) -> Result<String, RuntimeError> {
    let Some(value) = args.get(index) else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} missing argument at position {index}"
        )));
    };
    match value {
        MontyObject::String(text) => Ok(text.clone()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{function_name} expected str at position {index}, got {other:?}"
        ))),
    }
}

pub(crate) fn expect_rows(
    value: &MontyObject,
) -> Result<Vec<BTreeMap<String, String>>, RuntimeError> {
    let MontyObject::List(rows) = value else {
        return Err(RuntimeError::InvalidArguments(
            "write_tsv expects a list of dict rows".to_owned(),
        ));
    };

    let mut out = Vec::new();
    for row in rows {
        let MontyObject::Dict(dict) = row else {
            return Err(RuntimeError::InvalidArguments(
                "write_tsv row must be a dict".to_owned(),
            ));
        };
        let mut mapped = BTreeMap::new();
        for (key, value) in dict {
            let MontyObject::String(key) = key else {
                return Err(RuntimeError::InvalidArguments(
                    "write_tsv dict keys must be strings".to_owned(),
                ));
            };
            mapped.insert(key.clone(), stringify_value(value));
        }
        out.push(mapped);
    }
    Ok(out)
}

fn stringify_value(value: &MontyObject) -> String {
    match value {
        MontyObject::None => String::new(),
        MontyObject::String(text) => text.clone(),
        MontyObject::Int(v) => v.to_string(),
        MontyObject::Bool(v) => v.to_string(),
        other => format!("{other}"),
    }
}
