#![allow(
    clippy::missing_errors_doc,
    clippy::module_name_repetitions,
    clippy::must_use_candidate
)]

mod errors;
pub mod kestrel;
mod module_registry;
pub mod pyfaidx;
pub mod pysam;
pub mod samtools;
pub mod tools;
mod value;
pub mod vcf;

pub use errors::{LibError, LibResult};
pub use module_registry::{ModuleDescriptor, ModuleName, supported_modules};
pub use value::{LibValue, ObjectKind};
