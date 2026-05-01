use bioscript_core::{GenomicLocus, RuntimeError, VariantKind, VariantSpec};
use monty::MontyObject;

pub(crate) fn dataclass_handle_id(
    obj: &MontyObject,
    expected_name: &str,
) -> Result<u64, RuntimeError> {
    match obj {
        MontyObject::Dataclass { name, attrs, .. } if name == expected_name => {
            for (key, value) in attrs {
                if matches!(key, MontyObject::String(text) if text == "handle_id")
                    && let MontyObject::Int(id) = value
                {
                    return Ok(*id as u64);
                }
            }
            Err(RuntimeError::InvalidArguments(format!(
                "{expected_name} missing handle_id"
            )))
        }
        _ => Err(RuntimeError::InvalidArguments(format!(
            "expected {expected_name} object"
        ))),
    }
}

pub(crate) fn dataclass_to_variant_spec(obj: &MontyObject) -> Result<VariantSpec, RuntimeError> {
    let MontyObject::Dataclass { name, attrs, .. } = obj else {
        return Err(RuntimeError::InvalidArguments(
            "expected Variant object".to_owned(),
        ));
    };
    if name != "Variant" {
        return Err(RuntimeError::InvalidArguments(format!(
            "expected Variant object, got {name}"
        )));
    }

    let mut spec = VariantSpec::default();
    for (key, value) in attrs {
        let MontyObject::String(key) = key else {
            continue;
        };
        match key.as_str() {
            "rsids" => spec.rsids = string_list_from_object(value)?,
            "grch37" => {
                spec.grch37 = string_from_optional(value)?
                    .map(|v| parse_locus_string(&v))
                    .transpose()?
            }
            "grch38" => {
                spec.grch38 = string_from_optional(value)?
                    .map(|v| parse_locus_string(&v))
                    .transpose()?
            }
            "reference" => spec.reference = string_from_optional(value)?,
            "alternate" => spec.alternate = string_from_optional(value)?,
            "kind" => {
                spec.kind = string_from_optional(value)?
                    .as_deref()
                    .map(parse_variant_kind)
                    .transpose()?
            }
            "deletion_length" => {
                spec.deletion_length = int_from_optional(value)?.map(|v| v as usize)
            }
            "motifs" => spec.motifs = string_list_from_object(value)?,
            _ => {}
        }
    }
    Ok(spec)
}

pub(crate) fn variant_specs_from_plan(obj: &MontyObject) -> Result<Vec<VariantSpec>, RuntimeError> {
    match obj {
        MontyObject::List(items) => items.iter().map(dataclass_to_variant_spec).collect(),
        MontyObject::Dataclass { name, attrs, .. } if name == "VariantPlan" => {
            for (key, value) in attrs {
                if matches!(key, MontyObject::String(text) if text == "variants") {
                    return variant_specs_from_plan(value);
                }
            }
            Err(RuntimeError::InvalidArguments(
                "VariantPlan missing variants".to_owned(),
            ))
        }
        _ => Err(RuntimeError::InvalidArguments(
            "expected a list of Variant objects or a VariantPlan".to_owned(),
        )),
    }
}

pub(crate) fn variant_spec_from_kwargs(
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<VariantSpec, RuntimeError> {
    let mut spec = VariantSpec::default();
    for (key, value) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.variant keyword names must be strings".to_owned(),
            ));
        };
        match key.as_str() {
            "rsid" | "rsids" => spec.rsids = string_or_list(value)?,
            "grch37" => {
                spec.grch37 = string_from_optional(value)?
                    .map(|v| parse_locus_string(&v))
                    .transpose()?
            }
            "grch38" => {
                spec.grch38 = string_from_optional(value)?
                    .map(|v| parse_locus_string(&v))
                    .transpose()?
            }
            "ref" | "reference" => spec.reference = string_from_optional(value)?,
            "alt" | "alternate" => spec.alternate = string_from_optional(value)?,
            "kind" => {
                spec.kind = string_from_optional(value)?
                    .as_deref()
                    .map(parse_variant_kind)
                    .transpose()?
            }
            "deletion_length" => {
                spec.deletion_length = int_from_optional(value)?.map(|v| v as usize)
            }
            "motifs" => spec.motifs = string_or_list(value)?,
            other => {
                return Err(RuntimeError::InvalidArguments(format!(
                    "bioscript.variant does not accept keyword '{other}'"
                )));
            }
        }
    }
    Ok(spec)
}

fn parse_locus_string(value: &str) -> Result<GenomicLocus, RuntimeError> {
    let normalized = value.trim().strip_prefix("chr").unwrap_or(value.trim());
    let Some((chrom, rest)) = normalized.split_once(':') else {
        return Err(RuntimeError::InvalidArguments(format!(
            "invalid locus string: {value}"
        )));
    };
    let (start, end) = if let Some((start, end)) = rest.split_once('-') {
        (start, end)
    } else {
        (rest, rest)
    };
    let start = start.parse::<i64>().map_err(|err| {
        RuntimeError::InvalidArguments(format!("invalid locus start {value}: {err}"))
    })?;
    let end = end.parse::<i64>().map_err(|err| {
        RuntimeError::InvalidArguments(format!("invalid locus end {value}: {err}"))
    })?;
    Ok(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    })
}

fn parse_variant_kind(value: &str) -> Result<VariantKind, RuntimeError> {
    match value.trim().to_ascii_lowercase().as_str() {
        "snp" => Ok(VariantKind::Snp),
        "insertion" | "ins" => Ok(VariantKind::Insertion),
        "deletion" | "del" => Ok(VariantKind::Deletion),
        "indel" => Ok(VariantKind::Indel),
        "other" => Ok(VariantKind::Other),
        other => Err(RuntimeError::InvalidArguments(format!(
            "invalid variant kind: {other}"
        ))),
    }
}

pub(crate) fn string_or_list(value: &MontyObject) -> Result<Vec<String>, RuntimeError> {
    match value {
        MontyObject::String(text) => Ok(vec![text.clone()]),
        MontyObject::List(_) => string_list_from_object(value),
        MontyObject::None => Ok(Vec::new()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected string or list of strings, got {other:?}"
        ))),
    }
}

pub(crate) fn string_list_from_object(value: &MontyObject) -> Result<Vec<String>, RuntimeError> {
    match value {
        MontyObject::List(items) => {
            let mut out = Vec::new();
            for item in items {
                let MontyObject::String(text) = item else {
                    return Err(RuntimeError::InvalidArguments(
                        "expected list of strings".to_owned(),
                    ));
                };
                out.push(text.clone());
            }
            Ok(out)
        }
        MontyObject::None => Ok(Vec::new()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected list of strings, got {other:?}"
        ))),
    }
}

pub(crate) fn string_from_optional(value: &MontyObject) -> Result<Option<String>, RuntimeError> {
    match value {
        MontyObject::None => Ok(None),
        MontyObject::String(text) => Ok(Some(text.clone())),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected optional string, got {other:?}"
        ))),
    }
}

pub(crate) fn int_from_optional(value: &MontyObject) -> Result<Option<i64>, RuntimeError> {
    match value {
        MontyObject::None => Ok(None),
        MontyObject::Int(v) => Ok(Some(*v)),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected optional int, got {other:?}"
        ))),
    }
}
