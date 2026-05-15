use super::ReportOptionsInput;

/// Pull a head sample from the JS-backed reader and run it through the same
/// `inspect_bytes_rs` path the text/zip flow uses. This populates the
/// "Input" metadata block in the rust HTML report (Format / Source /
/// Assembly / Inferred sex / Evidence) for CRAM and VCF inputs, which
/// otherwise had no inspection data.
pub(super) fn inspect_head_via_js_reader(
    read_at: &js_sys::Function,
    total_len: u64,
    input_name: &str,
    options: &bioscript_formats::InspectOptions,
    detect_sex: bool,
) -> bioscript_formats::FileInspection {
    use crate::js_reader::JsReader;
    use std::io::Read;

    // For VCF/CRAM, sex detection and full assembly inference require
    // scanning many records. 8 MiB is enough for hundreds of decompressed
    // VCF records and several CRAM containers.
    let head_len = total_len.min(8 * 1024 * 1024);
    let mut reader = JsReader::new(read_at.clone(), total_len, "inspect");
    let mut buf = vec![0u8; head_len as usize];
    let mut filled = 0usize;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(_) => break,
        }
    }
    buf.truncate(filled);
    let mut opts = options.clone();
    opts.detect_sex = detect_sex;
    match bioscript_formats::inspect_bytes(input_name, &buf, &opts) {
        Ok(inspection) => inspection,
        Err(err) => bioscript_formats::FileInspection {
            path: std::path::PathBuf::from(input_name),
            container: bioscript_formats::FileContainer::Plain,
            detected_kind: bioscript_formats::DetectedKind::Unknown,
            confidence: bioscript_formats::DetectionConfidence::Unknown,
            source: None,
            assembly: None,
            phased: None,
            selected_entry: None,
            has_index: None,
            index_path: None,
            reference_matches: None,
            inferred_sex: None,
            evidence: vec![format!("inspect_bytes failed: {err:?}")],
            warnings: Vec::new(),
            duration_ms: 0,
        },
    }
}

/// Build the `SexInference` the CLI produces when `--sample-sex` is passed,
/// without dragging in the bioscript-cli crate. Mirrors
/// `bioscript_cli::report_options::explicit_sample_sex_inference`.
pub(super) fn explicit_sex_from_options(
    options: &ReportOptionsInput,
) -> Option<bioscript_formats::SexInference> {
    let raw = options.sample_sex.as_deref()?.trim().to_ascii_lowercase();
    let sex = match raw.as_str() {
        "male" | "m" => bioscript_formats::InferredSex::Male,
        "female" | "f" => bioscript_formats::InferredSex::Female,
        "unknown" | "u" | "" => bioscript_formats::InferredSex::Unknown,
        _ => return None,
    };
    Some(bioscript_formats::SexInference {
        sex,
        confidence: bioscript_formats::SexDetectionConfidence::High,
        method: "explicit_sample_sex".to_owned(),
        evidence: vec!["source=sample_sex_option".to_owned()],
    })
}

pub(super) fn vcf_sex_via_js_reader(
    read_at: js_sys::Function,
    total_len: u64,
    input_name: &str,
) -> Option<bioscript_formats::SexInference> {
    let reader = crate::js_reader::JsReader::new(read_at, total_len, "vcf-sex");
    bioscript_formats::infer_sex_from_named_reader(
        input_name,
        reader,
        bioscript_formats::DetectedKind::Vcf,
    )
    .ok()
}
