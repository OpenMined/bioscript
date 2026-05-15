// Keep included source files small and named by responsibility.
// If a file approaches 500 lines, split it by domain behavior rather than
// creating arbitrary numbered chunks.
include!("cli_bootstrap.rs");
include!("cli_commands.rs");
include!("report_options.rs");
include!("package.rs");
include!("report_review.rs");
include!("report_execution.rs");
include!("report_output.rs");
include!("manifest_runner.rs");
