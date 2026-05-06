// Keep validator source files small and grouped by schema responsibility.
// If a file approaches 500 lines, split it by validation domain rather than
// creating arbitrary numbered chunks.
include!("validator_types.rs");
include!("validator_load.rs");
include!("validator_roots.rs");
include!("validator_alleles_findings.rs");
include!("validator_panel.rs");
include!("validator_parse.rs");
include!("validator_helpers.rs");
