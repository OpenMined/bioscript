use crate::{LibError, LibResult};

use super::vcf::VariantCall;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceRegion {
    pub reference_name: String,
    pub sequence: String,
}

impl ReferenceRegion {
    pub fn base_at(&self, position: u32) -> LibResult<char> {
        if position == 0 {
            return Err(LibError::InvalidArguments(
                "Kestrel reference-region positions are 1-based".to_owned(),
            ));
        }
        self.sequence
            .chars()
            .nth(usize::try_from(position - 1).unwrap_or(usize::MAX))
            .ok_or_else(|| {
                LibError::InvalidArguments(format!(
                    "Kestrel reference position {position} is outside {}",
                    self.reference_name
                ))
            })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKind {
    Snp,
    Insertion,
    Deletion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeVariantCall {
    pub sample_name: String,
    pub kind: VariantKind,
    pub start: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub variant_depth: u32,
    pub locus_depth: u32,
}

impl NativeVariantCall {
    pub fn snp(
        sample_name: impl Into<String>,
        start: u32,
        ref_base: impl Into<String>,
        alt_base: impl Into<String>,
        variant_depth: u32,
        locus_depth: u32,
    ) -> Self {
        Self {
            sample_name: sample_name.into(),
            kind: VariantKind::Snp,
            start,
            ref_allele: ref_base.into(),
            alt_allele: alt_base.into(),
            variant_depth,
            locus_depth,
        }
    }

    pub fn insertion(
        sample_name: impl Into<String>,
        start: u32,
        inserted_bases: impl Into<String>,
        variant_depth: u32,
        locus_depth: u32,
    ) -> Self {
        Self {
            sample_name: sample_name.into(),
            kind: VariantKind::Insertion,
            start,
            ref_allele: String::new(),
            alt_allele: inserted_bases.into(),
            variant_depth,
            locus_depth,
        }
    }

    pub fn deletion(
        sample_name: impl Into<String>,
        start: u32,
        deleted_bases: impl Into<String>,
        variant_depth: u32,
        locus_depth: u32,
    ) -> Self {
        Self {
            sample_name: sample_name.into(),
            kind: VariantKind::Deletion,
            start,
            ref_allele: deleted_bases.into(),
            alt_allele: String::new(),
            variant_depth,
            locus_depth,
        }
    }

    pub fn to_vcf_call(&self, region: &ReferenceRegion) -> LibResult<VariantCall> {
        let (pos, ref_allele, alt_allele) = match self.kind {
            VariantKind::Snp => self.snp_vcf_fields()?,
            VariantKind::Insertion => self.insertion_vcf_fields(region)?,
            VariantKind::Deletion => self.deletion_vcf_fields(region)?,
        };
        Ok(VariantCall {
            sample_name: self.sample_name.clone(),
            chrom: region.reference_name.clone(),
            pos,
            ref_allele,
            alt_allele,
            variant_depth: self.variant_depth,
            locus_depth: self.locus_depth,
        })
    }

    fn snp_vcf_fields(&self) -> LibResult<(u32, String, String)> {
        if self.ref_allele.chars().count() != 1 || self.alt_allele.chars().count() != 1 {
            return Err(LibError::InvalidArguments(
                "Kestrel SNP REF and ALT must each be one base".to_owned(),
            ));
        }
        Ok((self.start, self.ref_allele.clone(), self.alt_allele.clone()))
    }

    fn insertion_vcf_fields(&self, region: &ReferenceRegion) -> LibResult<(u32, String, String)> {
        if self.alt_allele.is_empty() {
            return Err(LibError::InvalidArguments(
                "Kestrel insertion ALT cannot be empty".to_owned(),
            ));
        }
        if self.start == 0 {
            return Err(LibError::InvalidArguments(
                "Kestrel insertion start must be >= 1".to_owned(),
            ));
        }
        let anchor_pos = self.start.saturating_sub(1).max(1);
        let anchor = region.base_at(anchor_pos)?;
        let pos = if self.start == 1 { 1 } else { self.start - 1 };
        let alt = if self.start == 1 {
            format!("{}{anchor}", self.alt_allele)
        } else {
            format!("{anchor}{}", self.alt_allele)
        };
        Ok((pos, anchor.to_string(), alt))
    }

    fn deletion_vcf_fields(&self, region: &ReferenceRegion) -> LibResult<(u32, String, String)> {
        if self.ref_allele.is_empty() {
            return Err(LibError::InvalidArguments(
                "Kestrel deletion REF cannot be empty".to_owned(),
            ));
        }
        if self.start == 0 {
            return Err(LibError::InvalidArguments(
                "Kestrel deletion start must be >= 1".to_owned(),
            ));
        }
        if self.start == 1 {
            let anchor = region.base_at(self.reference_end() + 1)?;
            return Ok((
                1,
                format!("{}{anchor}", self.ref_allele),
                anchor.to_string(),
            ));
        }
        let anchor = region.base_at(self.start - 1)?;
        Ok((
            self.start - 1,
            format!("{anchor}{}", self.ref_allele),
            anchor.to_string(),
        ))
    }

    fn reference_end(&self) -> u32 {
        self.start + u32::try_from(self.ref_allele.chars().count()).unwrap_or(u32::MAX) - 1
    }
}
