use std::io::Read;

use bioscript_core::RuntimeError;
use noodles::{
    bcf, vcf,
    vcf::variant::record::{
        AlternateBases as _, Ids as _, ReferenceBases as _,
        samples::Sample as _,
        samples::series::{Value as SampleValue, value::Array},
    },
};

use super::{super::vcf::ParsedVcfRow, chrom_from_bcf_shard_name};

pub(super) fn read_bcf_header_lenient<R: Read>(
    reader: &mut bcf::io::Reader<R>,
    label: &str,
) -> Result<vcf::Header, RuntimeError> {
    let mut header_reader = reader.header_reader();
    header_reader
        .read_magic_number()
        .map_err(|err| RuntimeError::Io(format!("failed to read BCF magic {label}: {err}")))?;
    header_reader
        .read_format_version()
        .map_err(|err| RuntimeError::Io(format!("failed to read BCF version {label}: {err}")))?;
    let mut raw_reader = header_reader
        .raw_vcf_header_reader()
        .map_err(|err| RuntimeError::Io(format!("failed to open BCF VCF header {label}: {err}")))?;
    let mut raw = String::new();
    raw_reader
        .read_to_string(&mut raw)
        .map_err(|err| RuntimeError::Io(format!("failed to read BCF VCF header {label}: {err}")))?;
    raw_reader.discard_to_end().map_err(|err| {
        RuntimeError::Io(format!(
            "failed to discard BCF VCF header padding {label}: {err}"
        ))
    })?;

    let sanitized = sanitize_bcf_vcf_header(&raw);
    sanitized
        .parse::<vcf::Header>()
        .map_err(|err| RuntimeError::Io(format!("failed to parse BCF VCF header {label}: {err}")))
}

pub(super) fn bcf_record_to_vcf_row(
    header: &vcf::Header,
    string_maps: &vcf::header::StringMaps,
    record: &bcf::Record,
    label: &str,
) -> Result<ParsedVcfRow, RuntimeError> {
    let chrom = record
        .reference_sequence_name(string_maps)
        .map(str::to_owned)
        .or_else(|_| {
            chrom_from_bcf_shard_name(label)
                .map(|chrom| {
                    if chrom == "mt" {
                        "chrM".to_owned()
                    } else {
                        format!("chr{chrom}")
                    }
                })
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "missing reference sequence name in contig string map",
                    )
                })
        })
        .map_err(|err| {
            RuntimeError::Io(format!("{label}: failed to read BCF chromosome: {err}"))
        })?;
    let position = i64::try_from(
        record
            .variant_start()
            .transpose()
            .map_err(|err| {
                RuntimeError::Io(format!("{label}: failed to read BCF position: {err}"))
            })?
            .ok_or_else(|| RuntimeError::Io(format!("{label}: BCF record missing position")))?
            .get(),
    )
    .map_err(|_| RuntimeError::Io(format!("{label}: BCF position exceeds i64 range")))?;
    let ids_buf = record.ids();
    let ids: Vec<&str> = ids_buf.iter().collect();
    let rsid = ids
        .iter()
        .find(|id| !id.is_empty() && **id != ".")
        .map(|id| (*id).to_owned());
    let reference = String::from_utf8(
        record
            .reference_bases()
            .iter()
            .collect::<Result<Vec<_>, _>>()
            .map_err(|err| {
                RuntimeError::Io(format!(
                    "{label}: failed to read BCF reference bases: {err}"
                ))
            })?,
    )
    .map_err(|err| RuntimeError::Io(format!("{label}: invalid BCF reference bases: {err}")))?;
    let alternates = record
        .alternate_bases()
        .iter()
        .map(|result| result.map(str::to_owned))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|err| {
            RuntimeError::Io(format!(
                "{label}: failed to read BCF alternate bases: {err}"
            ))
        })?;
    let hds = extract_hds(header, record, label)?;
    let genotype = genotype_from_hds(&reference, alternates.first().map(String::as_str), &hds)
        .unwrap_or_else(|| "--".to_owned());
    let raw_line = format!(
        "{}  {}  {}  {}  {}  FORMAT=HDS  HDS={}",
        chrom,
        position,
        rsid.as_deref().unwrap_or("."),
        reference,
        alternates.join(","),
        hds.iter()
            .map(|value| format_hds_value(*value))
            .collect::<Vec<_>>()
            .join(",")
    );

    Ok(ParsedVcfRow {
        rsid,
        chrom,
        position,
        reference,
        alternates,
        genotype,
        raw_line,
    })
}

fn sanitize_bcf_vcf_header(raw: &str) -> String {
    let is_23andme_imputed = raw.contains("##DISCLAIMER=") && raw.contains("23andMe");
    let mut inserted_23andme_contigs = false;
    raw.lines()
        .filter_map(|line| {
            if is_23andme_imputed && line.starts_with("##contig=<") {
                if inserted_23andme_contigs {
                    return None;
                }
                inserted_23andme_contigs = true;
                return Some(
                    (1..=22)
                        .map(|chrom| format!("##contig=<ID=chr{chrom}>"))
                        .chain(std::iter::once("##contig=<ID=chrX>".to_owned()))
                        .collect::<Vec<_>>()
                        .join("\n"),
                );
            }
            if line.starts_with("##FORMAT=<")
                && line.contains("ID=HDS")
                && !line.contains("Description=")
            {
                Some(
                    "##FORMAT=<ID=HDS,Number=.,Type=Float,Description=\"haploid alternate dosage\">"
                        .to_owned(),
                )
            } else {
                Some(line.to_owned())
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn extract_hds(
    header: &vcf::Header,
    record: &bcf::Record,
    label: &str,
) -> Result<Vec<f32>, RuntimeError> {
    let samples = record
        .samples()
        .map_err(|err| RuntimeError::Io(format!("{label}: failed to read BCF samples: {err}")))?;
    if let Some(values) = parse_first_bcf_float_series(samples.as_ref()) {
        return Ok(values);
    }
    let Some(sample) = samples.get_index(0) else {
        return Ok(Vec::new());
    };
    let Some(value) = sample
        .get(header, "HDS")
        .transpose()
        .map_err(|err| RuntimeError::Io(format!("{label}: failed to read BCF HDS: {err}")))?
        .flatten()
    else {
        return Ok(Vec::new());
    };
    match value {
        SampleValue::Array(Array::Float(values)) => values
            .iter()
            .filter_map(|result| {
                result
                    .map_err(|err| {
                        RuntimeError::Io(format!("{label}: failed to read HDS value: {err}"))
                    })
                    .transpose()
            })
            .collect(),
        SampleValue::Float(value) => Ok(vec![value]),
        _ => Ok(Vec::new()),
    }
}

fn parse_first_bcf_float_series(src: &[u8]) -> Option<Vec<f32>> {
    let mut offset = 0usize;
    skip_bcf_typed_value(src, &mut offset)?;
    let descriptor = *src.get(offset)?;
    offset += 1;
    let (len, ty) = bcf_descriptor_len_ty(src, descriptor, &mut offset)?;
    if ty != 5 {
        return None;
    }
    let mut values = Vec::with_capacity(len);
    for _ in 0..len {
        let bytes: [u8; 4] = src.get(offset..offset + 4)?.try_into().ok()?;
        offset += 4;
        let value = f32::from_le_bytes(bytes);
        if !value.is_nan() {
            values.push(value);
        }
    }
    Some(values)
}

fn skip_bcf_typed_value(src: &[u8], offset: &mut usize) -> Option<()> {
    let descriptor = *src.get(*offset)?;
    *offset += 1;
    let (len, ty) = bcf_descriptor_len_ty(src, descriptor, offset)?;
    let width = match ty {
        1 | 7 => 1,
        2 => 2,
        3 | 5 => 4,
        _ => return None,
    };
    *offset = offset.checked_add(len.checked_mul(width)?)?;
    (*offset <= src.len()).then_some(())
}

fn bcf_descriptor_len_ty(src: &[u8], descriptor: u8, offset: &mut usize) -> Option<(usize, u8)> {
    let ty = descriptor & 0x0f;
    let mut len = usize::from(descriptor >> 4);
    if len == 15 {
        let len_descriptor = *src.get(*offset)?;
        *offset += 1;
        let (len_len, len_ty) = bcf_descriptor_len_ty(src, len_descriptor, offset)?;
        if len_len != 1 {
            return None;
        }
        len = match len_ty {
            1 => usize::from(*src.get(*offset)?),
            2 => {
                let bytes: [u8; 2] = src.get(*offset..*offset + 2)?.try_into().ok()?;
                usize::try_from(i16::from_le_bytes(bytes)).ok()?
            }
            3 => {
                let bytes: [u8; 4] = src.get(*offset..*offset + 4)?.try_into().ok()?;
                usize::try_from(i32::from_le_bytes(bytes)).ok()?
            }
            _ => return None,
        };
        *offset += match len_ty {
            1 => 1,
            2 => 2,
            3 => 4,
            _ => return None,
        };
    }
    Some((len, ty))
}

fn genotype_from_hds(reference: &str, alternate: Option<&str>, hds: &[f32]) -> Option<String> {
    let alternate = alternate?;
    if reference.is_empty() || alternate.is_empty() || hds.is_empty() {
        return None;
    }
    let mut alleles = Vec::with_capacity(hds.len().max(2));
    for dosage in hds {
        if *dosage >= 0.5 {
            alleles.push(alternate.to_owned());
        } else {
            alleles.push(reference.to_owned());
        }
    }
    if alleles.len() == 1 {
        alleles.push(reference.to_owned());
    }
    if alleles.iter().any(|allele| allele.len() > 1) {
        return Some(alleles.join("/"));
    }
    Some(super::super::normalize_genotype(&alleles.join("/")))
}

fn format_hds_value(value: f32) -> String {
    if (value.fract()).abs() < f32::EPSILON {
        format!("{value:.0}")
    } else {
        format!("{value:.3}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned()
    }
}
