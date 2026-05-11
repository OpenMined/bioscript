use std::path::{Path, PathBuf};

use crate::{
    LibResult,
    tools::{CommandSpec, path_arg},
};

pub mod native;

pub const MODULE: &str = "kestrel";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KestrelRunConfig {
    pub java_program: String,
    pub java_memory: String,
    pub jar_path: PathBuf,
    pub kmer_size: u16,
    pub max_align_states: u32,
    pub max_hap_states: u32,
    pub reference_vntr: PathBuf,
    pub output_vcf: PathBuf,
    pub output_sam: PathBuf,
    pub temp_dir: PathBuf,
    pub sample_name: String,
    pub fastq_1: PathBuf,
    pub fastq_2: PathBuf,
    pub log_level: String,
    pub additional_args: Vec<String>,
}

impl KestrelRunConfig {
    pub fn vntyper(
        jar_path: impl Into<PathBuf>,
        reference_vntr: impl Into<PathBuf>,
        output_vcf: impl Into<PathBuf>,
        output_sam: impl Into<PathBuf>,
        temp_dir: impl Into<PathBuf>,
        sample_name: impl Into<String>,
        fastq_1: impl Into<PathBuf>,
        fastq_2: impl Into<PathBuf>,
    ) -> Self {
        Self {
            java_program: "java".to_owned(),
            java_memory: "12g".to_owned(),
            jar_path: jar_path.into(),
            kmer_size: 20,
            max_align_states: 40,
            max_hap_states: 40,
            reference_vntr: reference_vntr.into(),
            output_vcf: output_vcf.into(),
            output_sam: output_sam.into(),
            temp_dir: temp_dir.into(),
            sample_name: sample_name.into(),
            fastq_1: fastq_1.into(),
            fastq_2: fastq_2.into(),
            log_level: "INFO".to_owned(),
            additional_args: Vec::new(),
        }
    }

    pub fn command(&self) -> LibResult<CommandSpec> {
        let mut args = vec![
            format!("-Xmx{}", self.java_memory),
            "-jar".to_owned(),
            path_arg(&self.jar_path)?,
            "-k".to_owned(),
            self.kmer_size.to_string(),
            "--maxalignstates".to_owned(),
            self.max_align_states.to_string(),
            "--maxhapstates".to_owned(),
            self.max_hap_states.to_string(),
            "-r".to_owned(),
            path_arg(&self.reference_vntr)?,
            "-o".to_owned(),
            path_arg(&self.output_vcf)?,
            format!("-s{}", self.sample_name),
            path_arg(&self.fastq_1)?,
            path_arg(&self.fastq_2)?,
            "--hapfmt".to_owned(),
            "sam".to_owned(),
            "-p".to_owned(),
            path_arg(&self.output_sam)?,
            "--logstderr".to_owned(),
            "--logstdout".to_owned(),
            "--loglevel".to_owned(),
            self.log_level.to_ascii_uppercase(),
            "--temploc".to_owned(),
            path_arg(&self.temp_dir)?,
        ];
        args.extend(self.additional_args.clone());
        CommandSpec::new(&self.java_program, args)
    }
}

pub fn read_vcf_command(path: &Path) -> LibResult<CommandSpec> {
    CommandSpec::new("bioscript-kestrel-vcf-reader", vec![path_arg(path)?])
}
