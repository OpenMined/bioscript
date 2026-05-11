use crate::{LibError, LibResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModuleName {
    Kestrel,
    Pysam,
    Pyfaidx,
    Samtools,
    Vcf,
}

impl ModuleName {
    pub fn parse(name: &str) -> LibResult<Self> {
        match name {
            "kestrel" => Ok(Self::Kestrel),
            "pysam" => Ok(Self::Pysam),
            "pyfaidx" => Ok(Self::Pyfaidx),
            "samtools" => Ok(Self::Samtools),
            "vcf" => Ok(Self::Vcf),
            other => Err(LibError::UnknownModule(other.to_owned())),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Kestrel => "kestrel",
            Self::Pysam => "pysam",
            Self::Pyfaidx => "pyfaidx",
            Self::Samtools => "samtools",
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
            name: ModuleName::Kestrel,
            import_path: "from bioscript import kestrel",
            summary: "structured Kestrel mapping-free variant-caller wrapper",
        },
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
            name: ModuleName::Samtools,
            import_path: "from bioscript import samtools",
            summary: "structured samtools command wrapper for allowed VNtyper verbs",
        },
        ModuleDescriptor {
            name: ModuleName::Vcf,
            import_path: "from bioscript import vcf",
            summary: "BioScript VCF compatibility namespace; may become pysam.VariantFile",
        },
    ]
}
