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
        let path = self.resolve_existing_user_path(&raw_path)?;
        let records = vcf::read_vntyper_kestrel_rows(&path)
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
}
