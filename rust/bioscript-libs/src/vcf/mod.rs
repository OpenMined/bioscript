use std::{collections::BTreeMap, fs, path::Path};

use crate::{LibError, LibResult};

pub const MODULE: &str = "vcf";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfDirection {
    PysamVariantFile,
}

pub fn chosen_initial_surface() -> VcfDirection {
    VcfDirection::PysamVariantFile
}

pub fn open_variant_file() -> LibResult<()> {
    Err(LibError::unsupported_feature(
        MODULE,
        "VariantFile; planned as bioscript.pysam.VariantFile first",
    ))
}

pub type VcfRecord = BTreeMap<String, String>;

pub fn read_kestrel_vcf(path: &Path) -> LibResult<Vec<VcfRecord>> {
    let contents = fs::read_to_string(path).map_err(|err| {
        LibError::InvalidArguments(format!("failed to read VCF {}: {err}", path.display()))
    })?;
    parse_kestrel_vcf(&contents)
}

pub fn parse_kestrel_vcf(contents: &str) -> LibResult<Vec<VcfRecord>> {
    let mut header: Option<Vec<String>> = None;
    let mut records = Vec::new();
    for line in contents.lines() {
        if line.trim().is_empty() || line.starts_with("##") {
            continue;
        }
        if let Some(header_line) = line.strip_prefix("#CHROM") {
            let mut names = vec!["CHROM".to_owned()];
            names.extend(
                header_line
                    .trim_start_matches('\t')
                    .split('\t')
                    .map(str::to_owned),
            );
            header = Some(names);
            continue;
        }
        let Some(header) = header.as_ref() else {
            continue;
        };
        let values = line.split('\t').collect::<Vec<_>>();
        let mut record = VcfRecord::new();
        for (idx, key) in header.iter().enumerate() {
            record.insert(
                key.clone(),
                values
                    .get(idx)
                    .map_or_else(String::new, |value| (*value).to_owned()),
            );
        }
        if let Some(sample) = record.get("SAMPLE").cloned()
            && !record.contains_key("Sample")
        {
            record.insert("Sample".to_owned(), sample);
        }
        records.push(record);
    }
    Ok(records)
}
