use std::cell::RefCell;
use std::io::{Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use flate2::read::MultiGzDecoder;
use kanalyze::comp::reader::{FileSequenceSource, SequenceFormat, SequenceRecord};
use kestrel::io::{InputSample, StreamableOutput};
use kestrel::runner::KestrelRunner;
use tempfile::TempDir;

use crate::{LibError, LibResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeReferenceRegion {
    pub reference_name: String,
    pub sequence: String,
    pub md5: String,
}

impl NativeReferenceRegion {
    pub fn new(
        reference_name: impl Into<String>,
        sequence: impl Into<String>,
        md5: impl Into<String>,
    ) -> Self {
        Self {
            reference_name: reference_name.into(),
            sequence: sequence.into(),
            md5: md5.into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NativeKestrelRunOptions {
    pub sample_name: String,
    pub minimum_difference: u32,
    pub difference_quantile: f32,
    pub anchor_both_ends: bool,
    pub decay_min: f32,
    pub decay_alpha: f32,
    pub peak_scan_length: usize,
    pub scan_limit_factor: f32,
    pub call_ambiguous_regions: bool,
    pub min_kmer_count: u32,
    pub max_haplotypes: usize,
    pub max_repeat_count: usize,
    pub max_saved_states: usize,
}

impl NativeKestrelRunOptions {
    pub fn new(sample_name: impl Into<String>) -> Self {
        Self {
            sample_name: sample_name.into(),
            minimum_difference: 5,
            difference_quantile: 0.90,
            anchor_both_ends: true,
            decay_min: 0.55,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            call_ambiguous_regions: true,
            min_kmer_count: 1,
            max_haplotypes: 40,
            max_repeat_count: 0,
            max_saved_states: 40,
        }
    }
}

pub fn call_sequences_to_vcf<'a>(
    reference_name: &str,
    reference_sequence: &str,
    read_sequences: impl IntoIterator<Item = &'a str>,
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<String> {
    let temp = tempfile::tempdir().map_err(io_error)?;
    let reference_path = temp.path().join("references.fasta");
    let fastq_path = temp.path().join("reads.fastq");
    write_reference_fasta(
        &reference_path,
        &[NativeReferenceRegion::new(
            reference_name,
            reference_sequence,
            ".",
        )],
    )?;
    write_reads_fastq(&fastq_path, read_sequences)?;
    run_kestrel_to_string(&temp, &[reference_path], &[fastq_path], kmer_size, options)
}

pub fn call_fastq_paths_to_vcf<'a>(
    reference_name: &str,
    reference_sequence: &str,
    fastq_paths: impl IntoIterator<Item = &'a Path>,
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<String> {
    let temp = tempfile::tempdir().map_err(io_error)?;
    let reference_path = temp.path().join("references.fasta");
    write_reference_fasta(
        &reference_path,
        &[NativeReferenceRegion::new(
            reference_name,
            reference_sequence,
            ".",
        )],
    )?;
    let fastq_paths = prepare_fastq_paths(&temp, fastq_paths)?;
    run_kestrel_to_string(&temp, &[reference_path], &fastq_paths, kmer_size, options)
}

pub fn call_fastq_paths_to_vcf_references<'a>(
    references: &[NativeReferenceRegion],
    fastq_paths: impl IntoIterator<Item = &'a Path>,
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<String> {
    let temp = tempfile::tempdir().map_err(io_error)?;
    let reference_path = temp.path().join("references.fasta");
    write_reference_fasta(&reference_path, references)?;
    let fastq_paths = prepare_fastq_paths(&temp, fastq_paths)?;
    run_kestrel_to_string(&temp, &[reference_path], &fastq_paths, kmer_size, options)
}

pub fn call_fastq_bytes_to_vcf_references(
    references: &[NativeReferenceRegion],
    fastq_inputs: &[(String, Vec<u8>)],
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<String> {
    let reference_sources = references
        .iter()
        .enumerate()
        .map(|(index, reference)| {
            FileSequenceSource::from_records(
                format!("reference_{}.fa", index + 1),
                SequenceFormat::Fasta,
                i32::try_from(index + 1).unwrap_or(i32::MAX),
                vec![SequenceRecord {
                    name: reference.reference_name.clone(),
                    sequence: reference.sequence.as_bytes().to_vec(),
                }],
            )
            .map_err(kestrel_error)
        })
        .collect::<LibResult<Vec<_>>>()?;
    let fastq_sources = fastq_inputs
        .iter()
        .enumerate()
        .map(|(index, (name, bytes))| {
            FileSequenceSource::from_records(
                name,
                SequenceFormat::Fastq,
                i32::try_from(index + 1).unwrap_or(i32::MAX),
                parse_fastq_records(name, bytes)?,
            )
            .map_err(kestrel_error)
        })
        .collect::<LibResult<Vec<_>>>()?;
    run_kestrel_sources_to_string(reference_sources, fastq_sources, kmer_size, options)
}

pub fn load_reference_regions(path: &Path) -> LibResult<Vec<NativeReferenceRegion>> {
    let content = std::fs::read_to_string(path).map_err(io_error)?;
    load_reference_regions_from_str(&content)
}

pub fn load_reference_regions_from_str(content: &str) -> LibResult<Vec<NativeReferenceRegion>> {
    let mut records = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_sequence = String::new();

    for raw_line in content.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) =
                current_name.replace(header.split_whitespace().next().unwrap_or("").to_owned())
            {
                records.push(NativeReferenceRegion::new(
                    name,
                    std::mem::take(&mut current_sequence),
                    ".",
                ));
            }
        } else {
            if current_name.is_none() {
                return Err(LibError::InvalidArguments(
                    "FASTA sequence appeared before a record header".to_owned(),
                ));
            }
            current_sequence.push_str(line);
        }
    }

    if let Some(name) = current_name {
        records.push(NativeReferenceRegion::new(name, current_sequence, "."));
    }
    if records.is_empty() {
        return Err(LibError::InvalidArguments(
            "FASTA content contains no records".to_owned(),
        ));
    }
    for record in &records {
        validate_name(&record.reference_name)?;
        validate_sequence(&record.sequence)?;
    }
    Ok(records)
}

fn parse_fastq_records(name: &str, bytes: &[u8]) -> LibResult<Vec<SequenceRecord>> {
    let data = if is_gzip_name(name) {
        let mut decoded = Vec::new();
        MultiGzDecoder::new(Cursor::new(bytes))
            .read_to_end(&mut decoded)
            .map_err(io_error)?;
        decoded
    } else {
        bytes.to_vec()
    };
    let mut records = Vec::new();
    let text = String::from_utf8(data).map_err(|err| {
        LibError::InvalidArguments(format!("FASTQ input {name} is not UTF-8: {err}"))
    })?;
    let mut lines = text.lines();
    while let Some(header) = lines.next() {
        let Some(sequence) = lines.next() else {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ input {name} has a truncated record"
            )));
        };
        let Some(plus) = lines.next() else {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ input {name} has a truncated record"
            )));
        };
        let Some(_qualities) = lines.next() else {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ input {name} has a truncated record"
            )));
        };
        if !header.starts_with('@') || !plus.starts_with('+') {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ input {name} has an invalid record"
            )));
        }
        records.push(SequenceRecord {
            name: header.trim_start_matches('@').to_owned(),
            sequence: sequence.as_bytes().to_vec(),
        });
    }
    Ok(records)
}

fn run_kestrel_sources_to_string(
    reference_sources: Vec<FileSequenceSource>,
    fastq_sources: Vec<FileSequenceSource>,
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<String> {
    let buffer = Rc::new(RefCell::new(Vec::new()));
    let mut runner = configured_runner_inner(None, Path::new("calls.vcf"), kmer_size, options)?;
    runner.set_output_file(Some(StreamableOutput::from_buffer(
        buffer.clone(),
        Some("calls.vcf"),
    )));

    for reference_source in reference_sources {
        runner.add_reference(reference_source);
    }
    runner.add_sample(
        InputSample::new(Some(&options.sample_name), fastq_sources).map_err(kestrel_error)?,
    );

    runner.run().map_err(kestrel_error)?;
    String::from_utf8(buffer.borrow().clone()).map_err(|err| {
        LibError::InvalidArguments(format!("Kestrel VCF output is not UTF-8: {err}"))
    })
}

fn run_kestrel_to_string(
    temp: &TempDir,
    reference_paths: &[PathBuf],
    fastq_paths: &[PathBuf],
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<String> {
    let output_path = temp.path().join("calls.vcf");
    let mut runner = configured_runner(temp, &output_path, kmer_size, options)?;

    for (index, reference_path) in reference_paths.iter().enumerate() {
        runner.add_reference(sequence_source(reference_path, index + 1)?);
    }

    let sources = fastq_paths
        .iter()
        .enumerate()
        .map(|(index, path)| sequence_source(path, index + 1))
        .collect::<LibResult<Vec<_>>>()?;
    runner
        .add_sample(InputSample::new(Some(&options.sample_name), sources).map_err(kestrel_error)?);

    runner.run().map_err(kestrel_error)?;
    std::fs::read_to_string(output_path).map_err(io_error)
}

fn configured_runner(
    temp: &TempDir,
    output_path: &Path,
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<KestrelRunner> {
    configured_runner_inner(
        Some(&temp.path().display().to_string()),
        output_path,
        kmer_size,
        options,
    )
}

fn configured_runner_inner(
    temp_dir_name: Option<&str>,
    output_path: &Path,
    kmer_size: usize,
    options: &NativeKestrelRunOptions,
) -> LibResult<KestrelRunner> {
    let mut runner = KestrelRunner::new();
    runner.set_k_size(kmer_size).map_err(kestrel_error)?;
    runner.set_output_path(output_path);
    runner.set_output_format("vcf").map_err(kestrel_error)?;
    runner.set_log_file(Some(StreamableOutput::stderr()));
    if let Some(temp_dir_name) = temp_dir_name {
        runner.set_temp_dir_name(Some(temp_dir_name));
    }
    runner.set_kmer_count_in_memory(true);
    runner.set_count_reverse_kmers(true);
    runner
        .set_minimum_difference(i32::try_from(options.minimum_difference).unwrap_or(i32::MAX))
        .map_err(kestrel_error)?;
    runner
        .set_difference_quantile(f64::from(options.difference_quantile))
        .map_err(kestrel_error)?;
    runner.set_anchor_both_ends(options.anchor_both_ends);
    runner
        .set_decay_minimum(f64::from(options.decay_min))
        .map_err(kestrel_error)?;
    runner
        .set_decay_alpha(f64::from(options.decay_alpha))
        .map_err(kestrel_error)?;
    runner
        .set_peak_scan_length(i32::try_from(options.peak_scan_length).unwrap_or(i32::MAX))
        .map_err(kestrel_error)?;
    runner
        .set_scan_limit_factor(f64::from(options.scan_limit_factor))
        .map_err(kestrel_error)?;
    runner.set_call_ambiguous_regions(options.call_ambiguous_regions);
    runner
        .set_min_kmer_count(i32::try_from(options.min_kmer_count).unwrap_or(i32::MAX))
        .map_err(kestrel_error)?;
    runner
        .set_max_haplotypes(i32::try_from(options.max_haplotypes).unwrap_or(i32::MAX))
        .map_err(kestrel_error)?;
    runner
        .set_max_repeat_count(i32::try_from(options.max_repeat_count).unwrap_or(i32::MAX))
        .map_err(kestrel_error)?;
    runner
        .set_max_aligner_state(i32::try_from(options.max_saved_states).unwrap_or(i32::MAX))
        .map_err(kestrel_error)?;
    Ok(runner)
}

fn write_reference_fasta(path: &Path, references: &[NativeReferenceRegion]) -> LibResult<()> {
    let mut file = std::fs::File::create(path).map_err(io_error)?;
    for reference in references {
        validate_name(&reference.reference_name)?;
        validate_sequence(&reference.sequence)?;
        writeln!(file, ">{}", reference.reference_name).map_err(io_error)?;
        writeln!(file, "{}", reference.sequence).map_err(io_error)?;
    }
    Ok(())
}

fn write_reads_fastq<'a>(
    path: &Path,
    read_sequences: impl IntoIterator<Item = &'a str>,
) -> LibResult<()> {
    let mut file = std::fs::File::create(path).map_err(io_error)?;
    for (index, sequence) in read_sequences.into_iter().enumerate() {
        validate_sequence(sequence)?;
        writeln!(file, "@read_{index}").map_err(io_error)?;
        writeln!(file, "{sequence}").map_err(io_error)?;
        writeln!(file, "+").map_err(io_error)?;
        writeln!(file, "{}", "I".repeat(sequence.len())).map_err(io_error)?;
    }
    Ok(())
}

fn sequence_source(path: &Path, source_id: usize) -> LibResult<FileSequenceSource> {
    FileSequenceSource::from_path(path, i32::try_from(source_id).unwrap_or(i32::MAX))
        .map_err(kestrel_error)
}

fn prepare_fastq_paths<'a>(
    temp: &TempDir,
    fastq_paths: impl IntoIterator<Item = &'a Path>,
) -> LibResult<Vec<PathBuf>> {
    fastq_paths
        .into_iter()
        .enumerate()
        .map(|(index, path)| {
            if is_gzip_path(path) {
                let output = temp.path().join(format!("input_{index}.fastq"));
                decompress_gzip(path, &output)?;
                Ok(output)
            } else {
                Ok(path.to_path_buf())
            }
        })
        .collect()
}

fn decompress_gzip(input: &Path, output: &Path) -> LibResult<()> {
    let input_file = std::fs::File::open(input).map_err(io_error)?;
    let mut reader = MultiGzDecoder::new(input_file);
    let mut writer = std::fs::File::create(output).map_err(io_error)?;
    std::io::copy(&mut reader, &mut writer).map_err(io_error)?;
    Ok(())
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
}

fn is_gzip_name(name: &str) -> bool {
    name.to_ascii_lowercase().ends_with(".gz")
}

fn validate_name(name: &str) -> LibResult<()> {
    if name.trim().is_empty() {
        return Err(LibError::InvalidArguments(
            "Kestrel reference name cannot be empty".to_owned(),
        ));
    }
    Ok(())
}

fn validate_sequence(sequence: &str) -> LibResult<()> {
    if sequence.is_empty() {
        return Err(LibError::InvalidArguments(
            "Kestrel sequence cannot be empty".to_owned(),
        ));
    }
    if !sequence.bytes().all(|base| base.is_ascii_alphabetic()) {
        return Err(LibError::InvalidArguments(
            "Kestrel sequence must contain only alphabetic bases".to_owned(),
        ));
    }
    Ok(())
}

fn kestrel_error(error: impl std::fmt::Display) -> LibError {
    LibError::InvalidArguments(format!("Kestrel error: {error}"))
}

fn io_error(error: impl std::fmt::Display) -> LibError {
    LibError::InvalidArguments(format!("Kestrel IO error: {error}"))
}
