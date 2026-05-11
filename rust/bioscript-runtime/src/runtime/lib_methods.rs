use std::path::PathBuf;

use bioscript_core::RuntimeError;
use bioscript_libs::{
    ModuleName, kestrel::KestrelRunConfig, pyfaidx::Fasta, pysam::AlignmentFile, samtools, vcf,
};
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{
        expect_int_arg, expect_string_arg, optional_string_kwarg, reject_kwargs,
        reject_unknown_kwargs,
    },
    objects::{
        kestrel_module_object, pyfaidx_fasta_object, pyfaidx_module_object,
        pysam_aligned_segment_object, pysam_alignment_file_object, pysam_module_object,
        samtools_module_object, vcf_module_object,
    },
};

pub(crate) fn host_bioscript_import(
    _runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "__bioscript_import__")?;
    let module = expect_string_arg(args, 0, "__bioscript_import__")?;
    match ModuleName::parse(&module).map_err(|err| RuntimeError::Unsupported(err.to_string()))? {
        ModuleName::Kestrel => Ok(kestrel_module_object()),
        ModuleName::Pysam => Ok(pysam_module_object()),
        ModuleName::Pyfaidx => Ok(pyfaidx_module_object()),
        ModuleName::Samtools => Ok(samtools_module_object()),
        ModuleName::Vcf => Ok(vcf_module_object()),
    }
}

impl BioscriptRuntime {
    pub(super) fn method_pysam_alignment_file(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_unknown_kwargs(
            kwargs,
            &["reference_filename", "index_filename"],
            "pysam.AlignmentFile",
        )?;
        if !(2..=3).contains(&args.len()) {
            return Err(RuntimeError::InvalidArguments(
                "pysam.AlignmentFile expects path and optional mode".to_owned(),
            ));
        }
        let path = expect_string_arg(args, 1, "pysam.AlignmentFile")?;
        let mode = if args.len() == 3 {
            expect_string_arg(args, 2, "pysam.AlignmentFile")?
        } else {
            "r".to_owned()
        };
        let reference_filename =
            optional_string_kwarg(kwargs, "reference_filename", "pysam.AlignmentFile")?
                .map(PathBuf::from);
        let index_filename =
            optional_string_kwarg(kwargs, "index_filename", "pysam.AlignmentFile")?
                .map(PathBuf::from);
        AlignmentFile::open(&path, &mode, reference_filename, index_filename)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        let reference_filename =
            optional_string_kwarg(kwargs, "reference_filename", "pysam.AlignmentFile")?;
        let index_filename =
            optional_string_kwarg(kwargs, "index_filename", "pysam.AlignmentFile")?;
        Ok(pysam_alignment_file_object(
            &path,
            &mode,
            reference_filename.as_deref(),
            index_filename.as_deref(),
        ))
    }

    pub(super) fn method_pysam_alignment_file_fetch(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "pysam.AlignmentFile.fetch")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "pysam.AlignmentFile.fetch expects contig, start, and stop".to_owned(),
            ));
        }
        let path = dataclass_string_attr(&args[0], "PysamAlignmentFile", "path")?;
        let mode = dataclass_string_attr(&args[0], "PysamAlignmentFile", "mode")?;
        let reference_filename =
            dataclass_optional_string_attr(&args[0], "PysamAlignmentFile", "reference_filename")?;
        let index_filename =
            dataclass_optional_string_attr(&args[0], "PysamAlignmentFile", "index_filename")?;
        let contig = expect_string_arg(args, 1, "pysam.AlignmentFile.fetch")?;
        let start = u64::try_from(expect_int_arg(args, 2, "pysam.AlignmentFile.fetch")?)
            .map_err(|_| RuntimeError::InvalidArguments("fetch start must be >= 0".to_owned()))?;
        let stop = u64::try_from(expect_int_arg(args, 3, "pysam.AlignmentFile.fetch")?)
            .map_err(|_| RuntimeError::InvalidArguments("fetch stop must be >= 0".to_owned()))?;
        let file = AlignmentFile::open(
            self.resolve_existing_user_path(&path)?,
            &mode,
            reference_filename
                .map(|path| self.resolve_existing_user_path(&path))
                .transpose()?,
            index_filename
                .map(|path| self.resolve_existing_user_path(&path))
                .transpose()?,
        )
        .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        let fetched = file
            .fetch(&contig, Some(start), Some(stop))
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        Ok(MontyObject::List(
            fetched
                .records
                .iter()
                .map(pysam_aligned_segment_object)
                .collect(),
        ))
    }

    pub(super) fn method_pyfaidx_fasta(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "pyfaidx.Fasta")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "pyfaidx.Fasta expects path".to_owned(),
            ));
        }
        let raw_path = expect_string_arg(args, 1, "pyfaidx.Fasta")?;
        let path = self.resolve_existing_user_path(&raw_path)?;
        Fasta::from_path(&path).map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        Ok(pyfaidx_fasta_object(&raw_path))
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
        command_argv_object(
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
        command_argv_object(
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
        command_argv_object(
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
        command_argv_object(
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
        command_argv_object(
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

fn command_argv_object(argv: Vec<String>) -> Result<MontyObject, RuntimeError> {
    Ok(MontyObject::List(
        argv.into_iter().map(MontyObject::String).collect(),
    ))
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

fn dataclass_string_attr(
    object: &MontyObject,
    expected_name: &str,
    attr_name: &str,
) -> Result<String, RuntimeError> {
    let Some(value) = dataclass_attr(object, expected_name, attr_name)? else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{expected_name}.{attr_name} is missing"
        )));
    };
    match value {
        MontyObject::String(value) => Ok(value.clone()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{expected_name}.{attr_name} expected str, got {other:?}"
        ))),
    }
}

fn dataclass_optional_string_attr(
    object: &MontyObject,
    expected_name: &str,
    attr_name: &str,
) -> Result<Option<String>, RuntimeError> {
    let Some(value) = dataclass_attr(object, expected_name, attr_name)? else {
        return Ok(None);
    };
    match value {
        MontyObject::None => Ok(None),
        MontyObject::String(value) => Ok(Some(value.clone())),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{expected_name}.{attr_name} expected str or None, got {other:?}"
        ))),
    }
}

fn dataclass_attr<'a>(
    object: &'a MontyObject,
    expected_name: &str,
    attr_name: &str,
) -> Result<Option<&'a MontyObject>, RuntimeError> {
    let MontyObject::Dataclass { name, attrs, .. } = object else {
        return Err(RuntimeError::InvalidArguments(format!(
            "expected {expected_name} object"
        )));
    };
    if name != expected_name {
        return Err(RuntimeError::InvalidArguments(format!(
            "expected {expected_name} object, got {name}"
        )));
    }
    Ok(attrs.into_iter().find_map(|(key, value)| {
        matches!(key, MontyObject::String(key) if key == attr_name).then_some(value)
    }))
}
