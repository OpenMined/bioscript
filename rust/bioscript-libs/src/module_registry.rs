use crate::{LibError, LibResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModuleName {
    Pysam,
    Pyfaidx,
    Vcf,
}

impl ModuleName {
    pub fn parse(name: &str) -> LibResult<Self> {
        match name {
            "pysam" => Ok(Self::Pysam),
            "pyfaidx" => Ok(Self::Pyfaidx),
            "vcf" => Ok(Self::Vcf),
            other => Err(LibError::UnknownModule(other.to_owned())),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pysam => "pysam",
            Self::Pyfaidx => "pyfaidx",
            Self::Vcf => "vcf",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModuleDescriptor {
    pub name: ModuleName,
    pub import_path: &'static str,
    pub summary: &'static str,
}

pub fn supported_modules() -> &'static [ModuleDescriptor] {
    &[
        ModuleDescriptor {
            name: ModuleName::Pysam,
            import_path: "from bioscript import pysam",
            summary: "pysam-compatible alignment and variant IO subset",
        },
        ModuleDescriptor {
            name: ModuleName::Pyfaidx,
            import_path: "from bioscript import pyfaidx",
            summary: "pyfaidx-compatible indexed FASTA subset",
        },
        ModuleDescriptor {
            name: ModuleName::Vcf,
            import_path: "from bioscript import vcf",
            summary: "BioScript VCF compatibility namespace; may become pysam.VariantFile",
        },
    ]
}
