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
    let opts = bioscript_formats::InspectOptions {
        input_index: None,
        reference_file: None,
        reference_index: None,
        detect_sex,
    };
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

const VCF_X_NON_PAR_WINDOWS_GRCH38: &[(i64, i64)] = &[
    (10_000_000, 11_000_000),
    (40_000_000, 41_000_000),
    (70_000_000, 71_000_000),
    (100_000_000, 101_000_000),
    (130_000_000, 131_000_000),
];
const VCF_Y_WINDOWS_GRCH38: &[(i64, i64)] = &[
    (3_500_000, 4_500_000),
    (10_000_000, 11_000_000),
    (15_000_000, 16_000_000),
];

/// Sex inference for indexed VCFs that streams only X non-PAR + Y windows
/// instead of scanning the whole file. Reuses the shared
/// `infer_sex_from_text_lines` so classification rules match the CLI.
pub(super) fn vcf_sex_via_tabix<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::csi::io::IndexedReader<
        noodles::bgzf::io::Reader<R>,
        noodles::tabix::Index,
    >,
    head_lines: &[String],
) -> Option<bioscript_formats::SexInference> {
    let mut lines = head_lines.to_vec();
    for chrom_label in ["X", "chrX"] {
        for (start, end) in VCF_X_NON_PAR_WINDOWS_GRCH38 {
            let Some(region) = build_region(chrom_label, *start, *end) else {
                continue;
            };
            if let Ok(query) = reader.query(&region) {
                for record_result in query {
                    let Ok(record) = record_result else {
                        continue;
                    };
                    let line: &str = record.as_ref();
                    lines.push(line.to_owned());
                }
            }
        }
        for (start, end) in VCF_Y_WINDOWS_GRCH38 {
            let y_label = if chrom_label == "X" { "Y" } else { "chrY" };
            let Some(region) = build_region(y_label, *start, *end) else {
                continue;
            };
            if let Ok(query) = reader.query(&region) {
                for record_result in query {
                    let Ok(record) = record_result else {
                        continue;
                    };
                    let line: &str = record.as_ref();
                    lines.push(line.to_owned());
                }
            }
        }
    }
    bioscript_formats::infer_sex_from_text_lines(&lines, bioscript_formats::DetectedKind::Vcf).ok()
}

fn build_region(chrom: &str, start: i64, end: i64) -> Option<noodles::core::Region> {
    use noodles::core::{Position, Region};
    let s = Position::try_from(usize::try_from(start.max(1)).ok()?).ok()?;
    let e = Position::try_from(usize::try_from(end.max(start)).ok()?).ok()?;
    Some(Region::new(chrom, s..=e))
}

/// Pull the VCF header from the bgzf head so `infer_sex_from_text_lines` can
/// resolve delimiter and column indexes for the X/Y records added via tabix.
pub(super) fn decompress_vcf_head_lines(read_at: &js_sys::Function, total_len: u64) -> Vec<String> {
    use crate::js_reader::JsReader;
    use std::io::{BufRead, BufReader, Read};

    let head_len = total_len.min(2 * 1024 * 1024);
    let mut reader = JsReader::new(read_at.clone(), total_len, "vcf-head");
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
    let cursor = std::io::Cursor::new(buf);
    let mut bgzf_reader = BufReader::new(noodles::bgzf::io::Reader::new(cursor));
    let mut lines = Vec::new();
    let mut line = String::new();
    for _ in 0..1024 {
        line.clear();
        match bgzf_reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => lines.push(line.trim_end_matches(['\n', '\r']).to_owned()),
            Err(_) => break,
        }
    }
    lines
}
