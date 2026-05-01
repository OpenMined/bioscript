use std::{env, fs, path::PathBuf};

use bioscript_formats::{
    GenotypeSourceFormat, InspectOptions, PrepareRequest, inspect_file, prepare_indexes,
    shell_flags,
};
use bioscript_schema::{validate_panels_path, validate_variants_path};

pub(crate) fn run_prepare(args: Vec<String>) -> Result<(), String> {
    let mut root: Option<PathBuf> = None;
    let mut input_file: Option<String> = None;
    let mut reference_file: Option<String> = None;
    let mut input_format: Option<GenotypeSourceFormat> = None;
    let mut cache_dir: Option<PathBuf> = None;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--root" => {
                root = Some(PathBuf::from(
                    iter.next().ok_or("--root requires a directory")?,
                ));
            }
            "--input-file" => {
                input_file = Some(iter.next().ok_or("--input-file requires a path")?);
            }
            "--reference-file" => {
                reference_file = Some(iter.next().ok_or("--reference-file requires a path")?);
            }
            "--input-format" => {
                let value = iter.next().ok_or("--input-format requires a value")?;
                if !value.eq_ignore_ascii_case("auto") {
                    input_format = Some(
                        value
                            .parse::<GenotypeSourceFormat>()
                            .map_err(|err| format!("invalid --input-format: {err}"))?,
                    );
                }
            }
            "--cache-dir" => {
                cache_dir = Some(PathBuf::from(
                    iter.next().ok_or("--cache-dir requires a path")?,
                ));
            }
            other => {
                return Err(format!("unexpected argument: {other}"));
            }
        }
    }

    let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
    let effective_root = root.unwrap_or_else(|| cwd.clone());
    let effective_cache = cache_dir.unwrap_or_else(|| cwd.join(".bioscript-cache"));

    let request = PrepareRequest {
        root: effective_root,
        cwd,
        cache_dir: effective_cache,
        input_file,
        input_format,
        reference_file,
    };

    let prepared = prepare_indexes(&request)?;
    let flags = shell_flags(&prepared);
    if flags.is_empty() {
        eprintln!("bioscript prepare: nothing to index");
    } else {
        println!("{flags}");
    }

    Ok(())
}

pub(crate) fn run_inspect(args: Vec<String>) -> Result<(), String> {
    let mut path: Option<PathBuf> = None;
    let mut options = InspectOptions::default();

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--input-index" => {
                options.input_index = Some(PathBuf::from(
                    iter.next().ok_or("--input-index requires a path")?,
                ));
            }
            "--reference-file" => {
                options.reference_file = Some(PathBuf::from(
                    iter.next().ok_or("--reference-file requires a path")?,
                ));
            }
            "--reference-index" => {
                options.reference_index = Some(PathBuf::from(
                    iter.next().ok_or("--reference-index requires a path")?,
                ));
            }
            other if path.is_none() => {
                path = Some(PathBuf::from(other));
            }
            other => {
                return Err(format!("unexpected argument: {other}"));
            }
        }
    }

    let Some(path) = path else {
        return Err(
            "usage: bioscript inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>]"
                .to_owned(),
        );
    };

    let inspection = inspect_file(&path, &options).map_err(|err| err.to_string())?;
    println!("{}", inspection.render_text());
    Ok(())
}

pub(crate) fn run_validate_variants(args: Vec<String>) -> Result<(), String> {
    run_validation_command(
        args,
        "usage: bioscript validate-variants <path> [--report <file>]",
        validate_variants_path,
    )
}

pub(crate) fn run_validate_panels(args: Vec<String>) -> Result<(), String> {
    run_validation_command(
        args,
        "usage: bioscript validate-panels <path> [--report <file>]",
        validate_panels_path,
    )
}

fn run_validation_command<F>(args: Vec<String>, usage: &str, validate: F) -> Result<(), String>
where
    F: FnOnce(&std::path::Path) -> Result<bioscript_schema::ValidationReport, String>,
{
    let mut path: Option<PathBuf> = None;
    let mut report_path: Option<PathBuf> = None;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        if arg == "--report" {
            let Some(value) = iter.next() else {
                return Err("--report requires a path".to_owned());
            };
            report_path = Some(PathBuf::from(value));
        } else if path.is_none() {
            path = Some(PathBuf::from(arg));
        } else {
            return Err(format!("unexpected argument: {arg}"));
        }
    }

    let Some(path) = path else {
        return Err(usage.to_owned());
    };

    let report = validate(&path)?;
    let text = report.render_text();
    print!("{text}");

    if let Some(report_path) = report_path {
        if let Some(parent) = report_path.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create report dir {}: {err}", parent.display())
            })?;
        }
        fs::write(&report_path, text)
            .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    }

    if report.has_errors() {
        return Err(format!(
            "validation found {} errors and {} warnings",
            report.total_errors(),
            report.total_warnings()
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioscript_schema::ValidationReport;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock drift")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "bioscript-cli-commands-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn empty_report() -> ValidationReport {
        ValidationReport {
            files_scanned: 1,
            reports: Vec::new(),
        }
    }

    #[test]
    fn command_parsers_report_prepare_and_inspect_argument_errors() {
        for (args, expected) in [
            (vec!["--root"], "--root requires a directory"),
            (vec!["--input-file"], "--input-file requires a path"),
            (vec!["--reference-file"], "--reference-file requires a path"),
            (vec!["--input-format"], "--input-format requires a value"),
            (vec!["--input-format", "bad"], "invalid --input-format"),
            (vec!["--cache-dir"], "--cache-dir requires a path"),
            (vec!["--unexpected"], "unexpected argument"),
        ] {
            let err = run_prepare(args.into_iter().map(str::to_owned).collect()).unwrap_err();
            assert!(err.contains(expected), "{err}");
        }

        for (args, expected) in [
            (Vec::<&str>::new(), "usage: bioscript inspect"),
            (vec!["--input-index"], "--input-index requires a path"),
            (vec!["--reference-file"], "--reference-file requires a path"),
            (
                vec!["--reference-index"],
                "--reference-index requires a path",
            ),
            (vec!["input.txt", "extra"], "unexpected argument"),
        ] {
            let err = run_inspect(args.into_iter().map(str::to_owned).collect()).unwrap_err();
            assert!(err.contains(expected), "{err}");
        }
    }

    #[test]
    fn validation_command_covers_report_success_and_error_paths() {
        let dir = temp_dir("validation");
        let input = dir.join("input.yaml");
        fs::write(&input, "schema: bioscript:variant:1.0\n").unwrap();
        let report = dir.join("reports/report.txt");

        run_validation_command(
            vec![
                input.display().to_string(),
                "--report".to_owned(),
                report.display().to_string(),
            ],
            "usage",
            |_| Ok(empty_report()),
        )
        .unwrap();
        assert!(
            fs::read_to_string(&report)
                .unwrap()
                .contains("files_scanned")
        );

        let err =
            run_validation_command(Vec::new(), "usage text", |_| Ok(empty_report())).unwrap_err();
        assert_eq!(err, "usage text");

        let err =
            run_validation_command(vec!["--report".to_owned()], "usage", |_| Ok(empty_report()))
                .unwrap_err();
        assert!(err.contains("--report requires a path"));

        let err = run_validation_command(vec!["one".to_owned(), "two".to_owned()], "usage", |_| {
            Ok(empty_report())
        })
        .unwrap_err();
        assert!(err.contains("unexpected argument"));

        let err = run_validation_command(vec!["input".to_owned()], "usage", |_| {
            Err("validator failed".to_owned())
        })
        .unwrap_err();
        assert_eq!(err, "validator failed");
    }

    #[test]
    fn public_validation_and_inspect_commands_cover_successful_argument_branches() {
        let dir = temp_dir("public-commands");
        let vcf = dir.join("sample.vcf");
        fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let index = dir.join("sample.vcf.tbi");
        let reference = dir.join("ref.fa");
        let reference_index = dir.join("ref.fa.fai");
        fs::write(&index, b"index").unwrap();
        fs::write(&reference, b">chr1\nA\n").unwrap();
        fs::write(&reference_index, b"chr1\t1\t6\t1\t2\n").unwrap();

        run_inspect(vec![
            vcf.display().to_string(),
            "--input-index".to_owned(),
            index.display().to_string(),
            "--reference-file".to_owned(),
            reference.display().to_string(),
            "--reference-index".to_owned(),
            reference_index.display().to_string(),
        ])
        .unwrap();

        let invalid_variant = dir.join("invalid-variant.yaml");
        fs::write(&invalid_variant, "schema: bioscript:variant:1.0\n").unwrap();
        let variant_report = dir.join("variant-report.txt");
        let err = run_validate_variants(vec![
            invalid_variant.display().to_string(),
            "--report".to_owned(),
            variant_report.display().to_string(),
        ])
        .unwrap_err();
        assert!(err.contains("validation found"));
        assert!(
            fs::read_to_string(&variant_report)
                .unwrap()
                .contains("errors:")
        );

        let invalid_panel = dir.join("invalid-panel.yaml");
        fs::write(&invalid_panel, "schema: bioscript:panel:1.0\n").unwrap();
        let panel_report = dir.join("panel-report.txt");
        let err = run_validate_panels(vec![
            invalid_panel.display().to_string(),
            "--report".to_owned(),
            panel_report.display().to_string(),
        ])
        .unwrap_err();
        assert!(err.contains("validation found"));
        assert!(
            fs::read_to_string(&panel_report)
                .unwrap()
                .contains("errors:")
        );
    }
}
