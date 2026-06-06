use bioscript_core::RuntimeError;
use bioscript_libs::vcf;
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{expect_string_arg, reject_kwargs},
};

impl BioscriptRuntime {
    pub(super) fn method_vcf_read_vntyper_kestrel(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "vcf.read_vntyper_kestrel")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "vcf.read_vntyper_kestrel expects path".to_owned(),
            ));
        }
        let raw_path = expect_string_arg(args, 1, "vcf.read_vntyper_kestrel")?;
        let raw_path_buf = std::path::PathBuf::from(&raw_path);
        let records = if let Some(contents) = self.read_virtual_text_file(&raw_path_buf) {
            vcf::read_vntyper_kestrel_rows_from_str(&contents)
        } else {
            let path = self.resolve_existing_user_path(&raw_path)?;
            vcf::read_vntyper_kestrel_rows(&path)
        }
        .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        Ok(MontyObject::List(
            records
                .into_iter()
                .map(|record| {
                    MontyObject::Dict(
                        record
                            .into_iter()
                            .map(|(key, value)| {
                                (MontyObject::String(key), MontyObject::String(value))
                            })
                            .collect(),
                    )
                })
                .collect(),
        ))
    }

    pub(super) fn method_vcf_build_vntyper_report_json(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "vcf.build_vntyper_report_json")?;
        if !(4..=6).contains(&args.len()) {
            return Err(RuntimeError::InvalidArguments(
                "vcf.build_vntyper_report_json expects sample_name, input_files, rows, optional metadata, and optional coverage".to_owned(),
            ));
        }
        let sample_name = expect_string_arg(args, 1, "vcf.build_vntyper_report_json")?;
        let input_files = string_dict(&args[2], "vcf.build_vntyper_report_json input_files")?;
        let rows = row_dicts(&args[3], "vcf.build_vntyper_report_json rows")?;
        let metadata = optional_string_dict(args, 4, "vcf.build_vntyper_report_json metadata")?;
        let coverage = optional_string_dict(args, 5, "vcf.build_vntyper_report_json coverage")?;
        let report = vcf::vntyper_report_json_with_context(
            &sample_name,
            &input_files,
            &rows,
            &metadata,
            &coverage,
        )
        .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        Ok(MontyObject::String(report))
    }
}

fn optional_string_dict(
    args: &[MontyObject],
    idx: usize,
    context: &str,
) -> Result<vcf::VcfRecord, RuntimeError> {
    match args.get(idx) {
        None | Some(MontyObject::None) => Ok(vcf::VcfRecord::new()),
        Some(value) => string_dict(value, context),
    }
}

fn row_dicts(value: &MontyObject, context: &str) -> Result<Vec<vcf::VcfRecord>, RuntimeError> {
    let MontyObject::List(rows) = value else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{context} expects list"
        )));
    };
    rows.iter()
        .map(|row| string_dict(row, context))
        .collect::<Result<Vec<_>, _>>()
}

fn string_dict(value: &MontyObject, context: &str) -> Result<vcf::VcfRecord, RuntimeError> {
    let MontyObject::Dict(items) = value else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{context} expects dict"
        )));
    };
    let mut out = vcf::VcfRecord::new();
    for (key, value) in items {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(format!(
                "{context} dict keys must be strings"
            )));
        };
        out.insert(key.clone(), monty_value_string(value));
    }
    Ok(out)
}

fn monty_value_string(value: &MontyObject) -> String {
    match value {
        MontyObject::None => String::new(),
        MontyObject::Bool(value) => if *value { "True" } else { "False" }.to_owned(),
        MontyObject::Int(value) => value.to_string(),
        MontyObject::Float(value) => value.to_string(),
        MontyObject::String(value) => value.clone(),
        other => format!("{other:?}"),
    }
}
