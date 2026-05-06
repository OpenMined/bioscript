//! Browser-facing bindings around the existing bioscript Rust code.
//! See docs/architecture/bioscript-is-source-of-truth.md — the app layer
//! must not reimplement file parsing or lookups in TS/JS. It goes through here.
//!
//! Current surface:
//! - `inspectBytes(name, bytes, options)` — file classification / vendor sniff
//! - `lookupCramVariants(cramReadAt, cramLen, craiBytes, fastaReadAt, fastaLen,
//!    faiBytes, variantsJson)` — SNP lookups against an indexed CRAM + FASTA
//!    through JS-supplied random-read callbacks.
//!
//! Pending (see migration checklist in the architecture doc):
//! - `loadGenotypesBytes(name, bytes)` / `lookupVariants(storeId, planJson)`
//! - `compileVariantYaml(yamlText)`
//! - Index-less fallback (linear scan or on-the-fly index build).
//! - Indel / deletion observations on CRAM.

mod inspect_api;
mod js_reader;
mod lookup_api;
mod variant_yaml;

pub use inspect_api::{inspect_bytes, resolve_remote_resource_text};
pub use lookup_api::{
    lookup_cram_variants, lookup_genotype_bytes_rsids, lookup_genotype_bytes_variants,
    lookup_vcf_variants,
};
pub use variant_yaml::compile_variant_yaml_text;

#[wasm_bindgen::prelude::wasm_bindgen(start)]
pub fn start() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}
