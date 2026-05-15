//! Browser-facing bindings around the existing bioscript Rust code.
//! See docs/architecture/bioscript-is-source-of-truth.md — the app layer
//! must not reimplement file parsing or lookups in TS/JS. It goes through here.
//!
//! Current surface:
//! - `inspectBytes(name, bytes, options)` — file classification / vendor sniff
//! - `lookupCramVariants(cramReadAt, cramLen, craiBytes, fastaReadAt, fastaLen,
//!   faiBytes, variantsJson)` — SNP lookups against an indexed CRAM + FASTA
//!   through JS-supplied random-read callbacks.
//!
//! Pending (see migration checklist in the architecture doc):
//! - `loadGenotypesBytes(name, bytes)` / `lookupVariants(storeId, planJson)`
//! - `compileVariantYaml(yamlText)`
//! - Index-less fallback (linear scan or on-the-fly index build).
//! - Indel / deletion observations on CRAM.

mod index_api;
mod inspect_api;
mod js_reader;
mod lookup_api;
mod package_api;
mod report_api;
mod variant_yaml;

pub use index_api::{
    generate_bam_bai, generate_bam_bai_from_reader, generate_cram_crai,
    generate_cram_crai_from_reader, generate_fasta_fai, generate_fasta_fai_from_reader,
    generate_vcf_tbi,
};
pub use inspect_api::{inspect_bytes, resolve_remote_resource_text};
pub use lookup_api::{
    lookup_cram_variants, lookup_genotype_bytes_rsids, lookup_genotype_bytes_variants,
    lookup_vcf_variants,
};
pub use package_api::{
    resolve_package_release_text, resolve_package_zip_bytes, verify_package_artifact_sha256,
};
pub use report_api::run_package_report_bytes;
pub use variant_yaml::compile_variant_yaml_text;

#[wasm_bindgen::prelude::wasm_bindgen(start)]
pub fn start() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}
