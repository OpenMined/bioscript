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

pub(crate) fn expect_int_arg(
    args: &[MontyObject],
    index: usize,
    function_name: &str,
) -> Result<i64, RuntimeError> {
    let Some(value) = args.get(index) else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} missing argument at position {index}"
        )));
    };
    match value {
        MontyObject::Int(value) => Ok(*value),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{function_name} expected int at position {index}, got {other:?}"
        ))),
    }
}

pub(crate) fn optional_string_kwarg(
    kwargs: &[(MontyObject, MontyObject)],
    name: &str,
    function_name: &str,
) -> Result<Option<String>, RuntimeError> {
    let mut found = None;
    for (key, value) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(format!(
                "{function_name} keyword names must be strings"
            )));
        };
        if key == name {
            if found.is_some() {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} got duplicate keyword argument {name}"
                )));
            }
            let MontyObject::String(value) = value else {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} expected keyword {name} to be str"
                )));
            };
            found = Some(value.clone());
        }
    }
    Ok(found)
}

pub(crate) fn optional_int_kwarg(
    kwargs: &[(MontyObject, MontyObject)],
    name: &str,
    function_name: &str,
) -> Result<Option<i64>, RuntimeError> {
    let mut found = None;
    for (key, value) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(format!(
                "{function_name} keyword names must be strings"
            )));
        };
        if key == name {
            if found.is_some() {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} got duplicate keyword argument {name}"
                )));
            }
            let MontyObject::Int(value) = value else {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} expected keyword {name} to be int"
                )));
            };
            found = Some(*value);
        }
    }
    Ok(found)
}

pub(crate) fn optional_float_kwarg(
    kwargs: &[(MontyObject, MontyObject)],
    name: &str,
    function_name: &str,
) -> Result<Option<f64>, RuntimeError> {
    let mut found = None;
    for (key, value) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(format!(
                "{function_name} keyword names must be strings"
            )));
        };
        if key == name {
            if found.is_some() {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} got duplicate keyword argument {name}"
                )));
            }
            // Accept int literals where a float is expected (e.g. 7 -> 7.0).
            let parsed = match value {
                MontyObject::Float(value) => *value,
                MontyObject::Int(value) => value
                    .to_string()
                    .parse::<f64>()
                    .expect("i64 string must parse as f64"),
                _ => {
                    return Err(RuntimeError::InvalidArguments(format!(
                        "{function_name} expected keyword {name} to be a number"
                    )));
                }
            };
            found = Some(parsed);
        }
    }
    Ok(found)
}

pub(crate) fn optional_bool_kwarg(
    kwargs: &[(MontyObject, MontyObject)],
    name: &str,
    function_name: &str,
) -> Result<Option<bool>, RuntimeError> {
    let mut found = None;
    for (key, value) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(format!(
                "{function_name} keyword names must be strings"
            )));
        };
        if key == name {
            if found.is_some() {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} got duplicate keyword argument {name}"
                )));
            }
            let MontyObject::Bool(value) = value else {
                return Err(RuntimeError::InvalidArguments(format!(
                    "{function_name} expected keyword {name} to be bool"
                )));
            };
            found = Some(*value);
        }
    }
    Ok(found)
}

pub(crate) fn reject_unknown_kwargs(
    kwargs: &[(MontyObject, MontyObject)],
    allowed: &[&str],
    function_name: &str,
) -> Result<(), RuntimeError> {
    for (key, _) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(format!(
                "{function_name} keyword names must be strings"
            )));
        };
        if !allowed.contains(&key.as_str()) {
            return Err(RuntimeError::InvalidArguments(format!(
                "{function_name} got unexpected keyword argument {key}"
            )));
        }
    }
    Ok(())
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
