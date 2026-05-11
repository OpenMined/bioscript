use std::{collections::HashMap, io::BufRead};

use bioscript_core::RuntimeError;

use super::{
    COMMENT_PREFIXES, GenotypeSourceFormat, GenotypeStore, QueryBackend, RowParser, RsidMapBackend,
    detect_delimiter, vcf_tokens::genotype_from_vcf_gt,
};

pub(crate) fn from_vcf_reader<R: BufRead>(
    mut reader: R,
    label: &str,
) -> Result<GenotypeStore, RuntimeError> {
    let mut values = HashMap::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        let bytes = reader
            .read_line(&mut buf)
            .map_err(|err| RuntimeError::Io(format!("failed to read {label}: {err}")))?;
        if bytes == 0 {
            break;
        }
        read_vcf_rsid_line(buf.trim_end_matches(['\n', '\r']), &mut values);
    }

    Ok(from_rsid_map(
        GenotypeSourceFormat::Vcf,
        values,
        HashMap::new(),
    ))
}

pub(crate) fn from_delimited_reader<R: BufRead>(
    format: GenotypeSourceFormat,
    mut reader: R,
    label: &str,
) -> Result<GenotypeStore, RuntimeError> {
    // Buffer lines up to the first non-empty/non-comment line so delimiter
    // detection sees representative input, then stream the rest directly.
    let mut prelude: Vec<String> = Vec::new();
    let mut buf = String::new();
    let mut delimiter = None;
    loop {
        buf.clear();
        let bytes = reader
            .read_line(&mut buf)
            .map_err(|err| RuntimeError::Io(format!("failed to read {label}: {err}")))?;
        if bytes == 0 {
            break;
        }
        let line = buf.trim_end_matches(['\n', '\r']).to_owned();
        let trimmed = line.trim();
        let is_data = !trimmed.is_empty()
            && !COMMENT_PREFIXES
                .iter()
                .any(|prefix| trimmed.starts_with(prefix));
        prelude.push(line);
        if is_data {
            delimiter = Some(detect_delimiter(&prelude));
            break;
        }
    }

    let mut parser = RowParser::new(delimiter.unwrap_or(super::Delimiter::Tab));
    let mut values = HashMap::new();
    let mut source_lines = HashMap::new();
    for line in prelude {
        consume_delimited_line(&mut parser, &line, &mut values, &mut source_lines)?;
    }
    loop {
        buf.clear();
        let bytes = reader
            .read_line(&mut buf)
            .map_err(|err| RuntimeError::Io(format!("failed to read {label}: {err}")))?;
        if bytes == 0 {
            break;
        }
        consume_delimited_line(
            &mut parser,
            buf.trim_end_matches(['\n', '\r']),
            &mut values,
            &mut source_lines,
        )?;
    }

    Ok(from_rsid_map(format, values, source_lines))
}

pub(crate) fn from_vcf_lines(lines: Vec<String>) -> Result<GenotypeStore, RuntimeError> {
    let mut values = HashMap::new();
    for line in lines {
        read_vcf_rsid_line(line.trim(), &mut values);
    }
    Ok(from_rsid_map(
        GenotypeSourceFormat::Vcf,
        values,
        HashMap::new(),
    ))
}

fn read_vcf_rsid_line(line: &str, values: &mut HashMap<String, String>) {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with("##") || trimmed.starts_with("#CHROM") {
        return;
    }

    let fields: Vec<&str> = trimmed.split('\t').collect();
    if fields.len() < 10 {
        return;
    }

    let rsid = fields[2].trim();
    if rsid.is_empty() || rsid == "." {
        return;
    }

    let reference = fields[3].trim();
    let alternates: Vec<&str> = fields[4]
        .split(',')
        .map(str::trim)
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .collect();
    if reference.is_empty() || alternates.is_empty() {
        return;
    }

    let sample_gt = fields[9].split(':').next().unwrap_or(".");
    if let Some(genotype) = genotype_from_vcf_gt(sample_gt, reference, &alternates) {
        values.insert(rsid.to_owned(), genotype);
    }
}

fn consume_delimited_line(
    parser: &mut RowParser,
    line: &str,
    values: &mut HashMap<String, String>,
    source_lines: &mut HashMap<String, String>,
) -> Result<(), RuntimeError> {
    let trimmed = line.trim().to_owned();
    if let Some((rsid, genotype)) = parser.consume_line(line)? {
        values.insert(rsid.clone(), genotype);
        source_lines.insert(rsid, trimmed);
    }
    Ok(())
}

fn from_rsid_map(
    format: GenotypeSourceFormat,
    values: HashMap<String, String>,
    source_lines: HashMap<String, String>,
) -> GenotypeStore {
    GenotypeStore {
        backend: QueryBackend::RsidMap(RsidMapBackend {
            format,
            values,
            source_lines,
        }),
    }
}
