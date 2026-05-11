use std::collections::{BTreeMap, HashMap};

use crate::{LibError, LibResult};

use super::variant::{NativeVariantCall, ReferenceRegion};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSequence {
    pub name: String,
    pub length: usize,
    pub md5: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantCall {
    pub sample_name: String,
    pub chrom: String,
    pub pos: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_depth: u32,
    pub locus_depth: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KestrelVcfWriter {
    source_version: String,
    references: Vec<ReferenceSequence>,
    sample_names: Vec<String>,
    records: BTreeMap<VcfRecordKey, HashMap<String, SampleDepth>>,
}

impl KestrelVcfWriter {
    pub fn new(source_version: impl Into<String>, references: Vec<ReferenceSequence>) -> Self {
        Self {
            source_version: source_version.into(),
            references,
            sample_names: Vec::new(),
            records: BTreeMap::new(),
        }
    }

    pub fn add_sample(&mut self, sample_name: impl Into<String>) -> LibResult<()> {
        let sample_name = sample_name.into();
        validate_sample_name(&sample_name)?;
        if self.sample_names.contains(&sample_name) {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel VCF sample already exists: {sample_name}"
            )));
        }
        self.sample_names.push(sample_name);
        Ok(())
    }

    pub fn add_variant(&mut self, variant: VariantCall) -> LibResult<()> {
        if !self.sample_names.contains(&variant.sample_name) {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel VCF variant references unknown sample: {}",
                variant.sample_name
            )));
        }
        validate_variant(&variant)?;
        let key = VcfRecordKey {
            chrom: variant.chrom,
            pos: variant.pos,
            ref_allele: variant.ref_allele,
            alt_allele: variant.alt_allele,
        };
        self.records.entry(key).or_default().insert(
            variant.sample_name,
            SampleDepth {
                variant_depth: variant.variant_depth,
                locus_depth: variant.locus_depth,
            },
        );
        Ok(())
    }

    pub fn add_native_variant(
        &mut self,
        variant: &NativeVariantCall,
        region: &ReferenceRegion,
    ) -> LibResult<()> {
        self.add_variant(variant.to_vcf_call(region)?)
    }

    pub fn to_vcf_string(&self) -> String {
        let mut out = String::new();
        out.push_str("##fileformat=VCF4.2\n");
        out.push_str(&format!("##source=Kestrel{}\n", self.source_version));
        for reference in &self.references {
            out.push_str(&format!(
                "##contig=<ID={},length={},md5={}>\n",
                reference.name, reference.length, reference.md5
            ));
        }
        out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        out.push_str("##FORMAT=<ID=GDP,Number=A,Type=Integer,Description=\"Estimated depth of all haplotypes supporting the alternate variant\">\n");
        out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Estimated depth of all haplotypes in the variant active region\">\n");
        out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for sample_name in &self.sample_names {
            out.push('\t');
            out.push_str(sample_name);
        }
        out.push('\n');
        for (key, sample_depths) in &self.records {
            out.push_str(&key.vcf_prefix());
            for sample_name in &self.sample_names {
                out.push('\t');
                if let Some(depth) = sample_depths.get(sample_name) {
                    out.push_str(&format!("1:{}:{}", depth.variant_depth, depth.locus_depth));
                } else {
                    out.push_str("0:.:.");
                }
            }
            out.push('\n');
        }
        out
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct VcfRecordKey {
    chrom: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
}

impl VcfRecordKey {
    fn vcf_prefix(&self) -> String {
        format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT:GDP:DP",
            self.chrom, self.pos, self.ref_allele, self.alt_allele
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SampleDepth {
    variant_depth: u32,
    locus_depth: u32,
}

fn validate_sample_name(sample_name: &str) -> LibResult<()> {
    if sample_name.is_empty() {
        return Err(LibError::InvalidArguments(
            "Kestrel VCF sample name cannot be empty".to_owned(),
        ));
    }
    if sample_name.chars().any(char::is_whitespace) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel VCF sample name cannot contain whitespace: {sample_name:?}"
        )));
    }
    Ok(())
}

fn validate_variant(variant: &VariantCall) -> LibResult<()> {
    if variant.chrom.is_empty() {
        return Err(LibError::InvalidArguments(
            "Kestrel VCF variant chromosome cannot be empty".to_owned(),
        ));
    }
    if variant.pos == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel VCF variant position must be >= 1".to_owned(),
        ));
    }
    if variant.ref_allele.is_empty() || variant.alt_allele.is_empty() {
        return Err(LibError::InvalidArguments(
            "Kestrel VCF variant REF and ALT cannot be empty".to_owned(),
        ));
    }
    if variant.locus_depth < variant.variant_depth {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel VCF locus depth {} is less than variant depth {}",
            variant.locus_depth, variant.variant_depth
        )));
    }
    Ok(())
}
