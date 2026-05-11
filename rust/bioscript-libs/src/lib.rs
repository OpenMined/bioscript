#![allow(
    clippy::missing_errors_doc,
    clippy::module_name_repetitions,
    clippy::must_use_candidate
)]

mod errors;
mod module_registry;
pub mod pyfaidx;
pub mod pysam;
mod value;
pub mod vcf;

pub use errors::{LibError, LibResult};
pub use module_registry::{ModuleDescriptor, ModuleName, supported_modules};
pub use value::{LibValue, ObjectKind};
