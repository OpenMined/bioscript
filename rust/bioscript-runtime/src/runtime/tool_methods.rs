use std::path::PathBuf;

use bioscript_core::RuntimeError;
use bioscript_libs::{bcftools, kestrel::KestrelRunConfig, samtools, vcf};
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{expect_string_arg, reject_kwargs},
    timing::RuntimeInstant,
};

impl BioscriptRuntime {
    pub(super) fn method_bcftools_sort(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.sort")?;
        if args.len() != 3 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.sort expects input_vcf and output_vcf_gz".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "bcftools.sort",
            started,
            bcftools::sort(
                PathBuf::from(expect_string_arg(args, 1, "bcftools.sort")?).as_path(),
                PathBuf::from(expect_string_arg(args, 2, "bcftools.sort")?).as_path(),
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_bcftools_index(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.index")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.index expects vcf_gz".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "bcftools.index",
            started,
            bcftools::index(PathBuf::from(expect_string_arg(args, 1, "bcftools.index")?).as_path())
                .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
                .argv(),
        )
    }

    pub(super) fn method_bcftools_view_filter(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.view_filter")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.view_filter expects input_vcf, output_vcf_gz, and include_expr"
                    .to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "bcftools.view_filter",
            started,
            bcftools::view_filter(
                PathBuf::from(expect_string_arg(args, 1, "bcftools.view_filter")?).as_path(),
                PathBuf::from(expect_string_arg(args, 2, "bcftools.view_filter")?).as_path(),
                &expect_string_arg(args, 3, "bcftools.view_filter")?,
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_bcftools_norm(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.norm")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.norm expects input_vcf, reference_fasta, and output_vcf_gz".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "bcftools.norm",
            started,
            bcftools::norm(
                PathBuf::from(expect_string_arg(args, 1, "bcftools.norm")?).as_path(),
                PathBuf::from(expect_string_arg(args, 2, "bcftools.norm")?).as_path(),
                PathBuf::from(expect_string_arg(args, 3, "bcftools.norm")?).as_path(),
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_bcftools_view_header_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.view_header_native")?;
        if args.len() != 3 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.view_header_native expects input_vcf and output_vcf".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let input = self.resolve_existing_user_path(&expect_string_arg(
            args,
            1,
            "bcftools.view_header_native",
        )?)?;
        let output = self.resolve_user_write_path(&expect_string_arg(
            args,
            2,
            "bcftools.view_header_native",
        )?)?;
        bcftools::view_header_native(&input, &output)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        native_tool_none(self, "bcftools.view_header_native", started)
    }

    pub(super) fn method_bcftools_view_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.view_native")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.view_native expects input_vcf, output_vcf, and output_type".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let input =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "bcftools.view_native")?)?;
        let output =
            self.resolve_user_write_path(&expect_string_arg(args, 2, "bcftools.view_native")?)?;
        let output_type = expect_string_arg(args, 3, "bcftools.view_native")?;
        bcftools::view_native(&input, &output, &output_type)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        native_tool_none(self, "bcftools.view_native", started)
    }

    pub(super) fn method_bcftools_sort_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.sort_native")?;
        if args.len() != 5 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.sort_native expects input_vcf, output_vcf, output_type, and write_index"
                    .to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let input =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "bcftools.sort_native")?)?;
        let output =
            self.resolve_user_write_path(&expect_string_arg(args, 2, "bcftools.sort_native")?)?;
        let output_type = expect_string_arg(args, 3, "bcftools.sort_native")?;
        let write_index = expect_bool_arg(args, 4, "bcftools.sort_native")?;
        bcftools::sort_native(&input, &output, &output_type, write_index)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        native_tool_none(self, "bcftools.sort_native", started)
    }

    pub(super) fn method_bcftools_index_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bcftools.index_native")?;
        if args.len() != 5 {
            return Err(RuntimeError::InvalidArguments(
                "bcftools.index_native expects input_vcf, output_index, tbi, and force".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let input =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "bcftools.index_native")?)?;
        let output =
            self.resolve_user_write_path(&expect_string_arg(args, 2, "bcftools.index_native")?)?;
        let tbi = expect_bool_arg(args, 3, "bcftools.index_native")?;
        let force = expect_bool_arg(args, 4, "bcftools.index_native")?;
        bcftools::index_native(&input, Some(&output), tbi, force)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        native_tool_none(self, "bcftools.index_native", started)
    }

    pub(super) fn method_kestrel_build_command(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "kestrel.build_command")?;
        if args.len() != 9 {
            return Err(RuntimeError::InvalidArguments(
                "kestrel.build_command expects jar_path, reference_vntr, output_vcf, output_sam, temp_dir, sample_name, fastq_1, and fastq_2".to_owned(),
            ));
        }
        let config = KestrelRunConfig::vntyper(
            expect_string_arg(args, 1, "kestrel.build_command")?,
            expect_string_arg(args, 2, "kestrel.build_command")?,
            expect_string_arg(args, 3, "kestrel.build_command")?,
            expect_string_arg(args, 4, "kestrel.build_command")?,
            expect_string_arg(args, 5, "kestrel.build_command")?,
            expect_string_arg(args, 6, "kestrel.build_command")?,
            expect_string_arg(args, 7, "kestrel.build_command")?,
            expect_string_arg(args, 8, "kestrel.build_command")?,
        );
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "kestrel.build_command",
            started,
            config
                .command()
                .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
                .argv(),
        )
    }

    pub(super) fn method_samtools_view_region(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.view_region")?;
        if args.len() != 5 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.view_region expects bam, region, output_bam, and include_unmapped"
                    .to_owned(),
            ));
        }
        let include_unmapped = expect_bool_arg(args, 4, "samtools.view_region")?;
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.view_region",
            started,
            samtools::view_region(
                PathBuf::from(expect_string_arg(args, 1, "samtools.view_region")?).as_path(),
                &expect_string_arg(args, 2, "samtools.view_region")?,
                PathBuf::from(expect_string_arg(args, 3, "samtools.view_region")?).as_path(),
                include_unmapped,
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_samtools_fastq(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.fastq")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.fastq expects bam, fastq_1, and fastq_2".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.fastq",
            started,
            samtools::fastq(
                PathBuf::from(expect_string_arg(args, 1, "samtools.fastq")?).as_path(),
                PathBuf::from(expect_string_arg(args, 2, "samtools.fastq")?).as_path(),
                PathBuf::from(expect_string_arg(args, 3, "samtools.fastq")?).as_path(),
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_samtools_depth(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.depth")?;
        if args.len() != 3 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.depth expects bam and region".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.depth",
            started,
            samtools::depth(
                PathBuf::from(expect_string_arg(args, 1, "samtools.depth")?).as_path(),
                &expect_string_arg(args, 2, "samtools.depth")?,
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_samtools_index(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.index")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.index expects bam".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.index",
            started,
            samtools::index(PathBuf::from(expect_string_arg(args, 1, "samtools.index")?).as_path())
                .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
                .argv(),
        )
    }

    pub(super) fn method_vcf_variant_file(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "vcf.VariantFile")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "vcf.VariantFile expects path".to_owned(),
            ));
        }
        vcf::open_variant_file().map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        Ok(MontyObject::None)
    }

    pub(super) fn method_vcf_read_kestrel(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "vcf.read_kestrel")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "vcf.read_kestrel expects path".to_owned(),
            ));
        }
        let raw_path = expect_string_arg(args, 1, "vcf.read_kestrel")?;
        let path = self.resolve_existing_user_path(&raw_path)?;
        let records = vcf::read_kestrel_vcf(&path)
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

fn command_argv_object(
    runtime: &BioscriptRuntime,
    method: &str,
    started: RuntimeInstant,
    argv: Vec<String>,
) -> Result<MontyObject, RuntimeError> {
    runtime.record_timing(
        "tool_command_plan",
        started.elapsed(),
        format!("method={method} argv={}", argv.join(" ")),
    );
    Ok(MontyObject::List(
        argv.into_iter().map(MontyObject::String).collect(),
    ))
}

fn native_tool_none(
    runtime: &BioscriptRuntime,
    method: &str,
    started: RuntimeInstant,
) -> Result<MontyObject, RuntimeError> {
    runtime.record_timing(
        "native_tool_call",
        started.elapsed(),
        format!("method={method}"),
    );
    Ok(MontyObject::None)
}

fn expect_bool_arg(
    args: &[MontyObject],
    index: usize,
    function_name: &str,
) -> Result<bool, RuntimeError> {
    let Some(value) = args.get(index) else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} missing argument at position {index}"
        )));
    };
    match value {
        MontyObject::Bool(value) => Ok(*value),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{function_name} expected bool at position {index}, got {other:?}"
        ))),
    }
}
