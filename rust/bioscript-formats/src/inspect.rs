use std::path::{Path, PathBuf};

#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

// std::time::Instant::now() panics on wasm32-unknown-unknown ("time not
// implemented on this platform"). duration_ms is diagnostic-only, so on wasm
// we stub the timer instead of pulling in a perf/date shim crate.
#[cfg(target_arch = "wasm32")]
struct Instant;

#[cfg(target_arch = "wasm32")]
impl Instant {
    fn now() -> Self {
        Self
    }
    fn elapsed(&self) -> StubDuration {
        StubDuration
    }
}

#[cfg(target_arch = "wasm32")]
struct StubDuration;

#[cfg(target_arch = "wasm32")]
impl StubDuration {
    fn as_millis(&self) -> u128 {
        0
    }
}

use bioscript_core::{Assembly, RuntimeError};

mod heuristics;
mod io;
mod render;

pub(crate) use heuristics::*;
pub(crate) use io::*;
#[cfg(test)]
pub(crate) use render::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileContainer {
    Plain,
    Zip,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedKind {
    GenotypeText,
    Vcf,
    AlignmentCram,
    AlignmentBam,
    ReferenceFasta,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectionConfidence {
    Authoritative,
    StrongHeuristic,
    WeakHeuristic,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceMetadata {
    pub vendor: Option<String>,
    pub platform_version: Option<String>,
    pub confidence: DetectionConfidence,
    pub evidence: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct InspectOptions {
    pub input_index: Option<PathBuf>,
    pub reference_file: Option<PathBuf>,
    pub reference_index: Option<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileInspection {
    pub path: PathBuf,
    pub container: FileContainer,
    pub detected_kind: DetectedKind,
    pub confidence: DetectionConfidence,
    pub source: Option<SourceMetadata>,
    pub assembly: Option<Assembly>,
    pub phased: Option<bool>,
    pub selected_entry: Option<String>,
    pub has_index: Option<bool>,
    pub index_path: Option<PathBuf>,
    pub reference_matches: Option<bool>,
    pub evidence: Vec<String>,
    pub warnings: Vec<String>,
    pub duration_ms: u128,
}

/// Classify a file from in-memory bytes. Mirrors `inspect_file` but sources
/// its sample lines / zip entries from a byte buffer instead of the
/// filesystem. Needed by wasm targets where `std::fs` isn't available.
///
/// `name` is used for extension-based detection (.cram / .bam / .fa / .zip /
/// .vcf.gz) and vendor sniffing from the filename. `bytes` is the file
/// contents; for zips we scan the central directory out of these bytes.
pub fn inspect_bytes(
    name: &str,
    bytes: &[u8],
    options: &InspectOptions,
) -> Result<FileInspection, RuntimeError> {
    let started = Instant::now();
    let lower = name.to_ascii_lowercase();
    let mut evidence = Vec::new();
    let mut warnings = Vec::new();
    let path = Path::new(name);

    if lower.ends_with(".zip") {
        let selected_entry = select_zip_entry_from_bytes(bytes)?;
        let sample_lines = read_zip_sample_lines_from_bytes(bytes, &selected_entry)?;
        let mut inspection = inspect_from_textual_sample(
            path,
            FileContainer::Zip,
            &selected_entry,
            &sample_lines,
            options,
        );
        inspection.duration_ms = started.elapsed().as_millis();
        return Ok(inspection);
    }

    let detected_kind = if lower.ends_with(".cram") {
        evidence.push("extension .cram".to_owned());
        DetectedKind::AlignmentCram
    } else if lower.ends_with(".bam") {
        evidence.push("extension .bam".to_owned());
        DetectedKind::AlignmentBam
    } else if is_reference_path(path) {
        evidence.push("reference fasta extension".to_owned());
        DetectedKind::ReferenceFasta
    } else {
        let sample_lines = read_plain_sample_lines_from_bytes(&lower, bytes)?;
        let sample_lower = sample_lines.join("\n").to_ascii_lowercase();
        if looks_like_vcf_lines(&sample_lines) {
            evidence.push("vcf header markers".to_owned());
            DetectedKind::Vcf
        } else if looks_like_genotype_text(&sample_lines) {
            if sample_lower.contains("rsid") || sample_lower.contains("allele1") {
                evidence.push("genotype-like sampled rows and headers".to_owned());
            } else {
                evidence.push("genotype-like sampled rows".to_owned());
            }
            DetectedKind::GenotypeText
        } else {
            warnings.push("file did not match known textual heuristics".to_owned());
            DetectedKind::Unknown
        }
    };

    let sample_lines = match detected_kind {
        DetectedKind::AlignmentCram | DetectedKind::AlignmentBam | DetectedKind::ReferenceFasta => {
            Vec::new()
        }
        _ => read_plain_sample_lines_from_bytes(&lower, bytes)?,
    };
    let source = detect_source(&lower, &sample_lines, detected_kind);
    let assembly = detect_assembly(&lower, &sample_lines);
    let phased = (detected_kind == DetectedKind::Vcf)
        .then(|| detect_vcf_phasing(&sample_lines))
        .flatten();
    // Index discovery is filesystem-only; wasm callers pass indexes separately
    // through `InspectOptions.input_index` / `reference_index` and we surface
    // whichever is provided.
    let has_index = options
        .input_index
        .as_ref()
        .or(options.reference_index.as_ref())
        .map(|_| true);
    let index_path = options
        .input_index
        .clone()
        .or_else(|| options.reference_index.clone());
    let confidence = classify_confidence(detected_kind, &sample_lines, source.as_ref());

    Ok(FileInspection {
        path: path.to_path_buf(),
        container: FileContainer::Plain,
        detected_kind,
        confidence,
        source,
        assembly,
        phased,
        selected_entry: None,
        has_index,
        index_path,
        reference_matches: None,
        evidence,
        warnings,
        duration_ms: started.elapsed().as_millis(),
    })
}

pub fn inspect_file(path: &Path, options: &InspectOptions) -> Result<FileInspection, RuntimeError> {
    let started = Instant::now();
    let lower = path.to_string_lossy().to_ascii_lowercase();
    let mut evidence = Vec::new();
    let mut warnings = Vec::new();

    if lower.ends_with(".zip") {
        let selected_entry = select_zip_entry(path)?;
        let sample_lines = read_zip_sample_lines(path, &selected_entry)?;
        let mut inspection = inspect_from_textual_sample(
            path,
            FileContainer::Zip,
            &selected_entry,
            &sample_lines,
            options,
        );
        inspection.duration_ms = started.elapsed().as_millis();
        return Ok(inspection);
    }

    let detected_kind = if lower.ends_with(".cram") {
        evidence.push("extension .cram".to_owned());
        DetectedKind::AlignmentCram
    } else if lower.ends_with(".bam") {
        evidence.push("extension .bam".to_owned());
        DetectedKind::AlignmentBam
    } else if is_reference_path(path) {
        evidence.push("reference fasta extension".to_owned());
        DetectedKind::ReferenceFasta
    } else {
        let sample_lines = read_plain_sample_lines(path)?;
        let sample_lower = sample_lines.join("\n").to_ascii_lowercase();
        if looks_like_vcf_lines(&sample_lines) {
            evidence.push("vcf header markers".to_owned());
            DetectedKind::Vcf
        } else if looks_like_genotype_text(&sample_lines) {
            if sample_lower.contains("rsid") || sample_lower.contains("allele1") {
                evidence.push("genotype-like sampled rows and headers".to_owned());
            } else {
                evidence.push("genotype-like sampled rows".to_owned());
            }
            DetectedKind::GenotypeText
        } else {
            warnings.push("file did not match known textual heuristics".to_owned());
            DetectedKind::Unknown
        }
    };

    let sample_lines = match detected_kind {
        DetectedKind::AlignmentCram | DetectedKind::AlignmentBam | DetectedKind::ReferenceFasta => {
            Vec::new()
        }
        _ => read_plain_sample_lines(path)?,
    };
    let source = detect_source(&lower, &sample_lines, detected_kind);
    let assembly = detect_assembly(&lower, &sample_lines);
    let phased = (detected_kind == DetectedKind::Vcf)
        .then(|| detect_vcf_phasing(&sample_lines))
        .flatten();
    let (has_index, index_path) = detect_index(path, detected_kind, options);
    let confidence = classify_confidence(detected_kind, &sample_lines, source.as_ref());

    Ok(FileInspection {
        path: path.to_path_buf(),
        container: FileContainer::Plain,
        detected_kind,
        confidence,
        source,
        assembly,
        phased,
        selected_entry: None,
        has_index,
        index_path,
        reference_matches: None,
        evidence,
        warnings,
        duration_ms: started.elapsed().as_millis(),
    })
}

fn inspect_from_textual_sample(
    path: &Path,
    container: FileContainer,
    selected_entry: &str,
    sample_lines: &[String],
    options: &InspectOptions,
) -> FileInspection {
    let lower = selected_entry.to_ascii_lowercase();
    let detected_kind = if lower.ends_with(".vcf")
        || lower.ends_with(".vcf.gz")
        || looks_like_vcf_lines(sample_lines)
    {
        DetectedKind::Vcf
    } else if looks_like_genotype_text(sample_lines) {
        DetectedKind::GenotypeText
    } else {
        DetectedKind::Unknown
    };
    let path_lower = path.to_string_lossy().to_ascii_lowercase();
    let combined_name = format!("{path_lower}\n{lower}");
    let source = detect_source(&combined_name, sample_lines, detected_kind);
    let assembly = detect_assembly(&combined_name, sample_lines);
    let phased = (detected_kind == DetectedKind::Vcf)
        .then(|| detect_vcf_phasing(sample_lines))
        .flatten();
    let confidence = classify_confidence(detected_kind, sample_lines, source.as_ref());
    let mut evidence = vec![format!("selected zip entry {selected_entry}")];
    if detected_kind == DetectedKind::Vcf {
        evidence.push("vcf entry markers".to_owned());
    } else if detected_kind == DetectedKind::GenotypeText {
        evidence.push("genotype-like sampled rows".to_owned());
    }
    let (has_index, index_path) = detect_index(path, detected_kind, options);

    FileInspection {
        path: path.to_path_buf(),
        container,
        detected_kind,
        confidence,
        source,
        assembly,
        phased,
        selected_entry: Some(selected_entry.to_owned()),
        has_index,
        index_path,
        reference_matches: None,
        evidence,
        warnings: Vec::new(),
        duration_ms: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::bgzf;
    use std::io::{Cursor, Write as _};
    use std::path::PathBuf;

    #[test]
    fn inspect_zip_entry_limited_reader_rejects_oversized_output() {
        let mut reader = Cursor::new(b"abcdef".to_vec());
        let err = read_entry_limited(&mut reader, 5, "inspect zip entry").unwrap_err();
        assert!(
            err.to_string()
                .contains("inspect zip entry exceeds decompressed limit of 5 bytes"),
            "{err}"
        );
    }

    #[test]
    fn inspect_helpers_cover_text_shape_source_and_assembly_edges() {
        assert_eq!(split_fields("rs1\t1\t2\tAA"), vec!["rs1", "1", "2", "AA"]);
        assert_eq!(
            split_fields("\"rs1\", 1, 2, \"AG\""),
            vec!["rs1", "1", "2", "AG"]
        );
        assert!(looks_like_genotype_text(&[
            "// header".to_owned(),
            "i12345 XY 10 A G".to_owned(),
            "rs2 chr26 20 DD".to_owned(),
        ]));
        assert!(!looks_like_genotype_text(&["not enough fields".to_owned()]));
        assert!(!matches_genotype_shape(&[
            "bad".to_owned(),
            "1".to_owned(),
            "2".to_owned(),
            "AA".to_owned()
        ]));
        assert!(!matches_genotype_shape(&[
            "rs1".to_owned(),
            "badchr".to_owned(),
            "2".to_owned(),
            "AA".to_owned()
        ]));
        assert!(!is_valid_genotype(""));
        assert!(!is_valid_genotype("ACGTI"));
        assert!(!is_valid_allele("N"));

        let gfg = detect_source(
            "genesforgood.txt",
            &["# Genes for Good v1 export".to_owned()],
            DetectedKind::GenotypeText,
        )
        .unwrap();
        assert_eq!(gfg.vendor.as_deref(), Some("Genes for Good"));
        assert_eq!(gfg.platform_version.as_deref(), Some("v1"));

        let twenty_three =
            detect_source("/tmp/v5/23andme.txt", &[], DetectedKind::GenotypeText).unwrap();
        assert_eq!(twenty_three.vendor.as_deref(), Some("23andMe"));
        assert_eq!(twenty_three.platform_version.as_deref(), Some("v5"));
        assert_eq!(
            detect_source(
                "sequencing.com.vcf",
                &["##source=sequencing.com".to_owned()],
                DetectedKind::Vcf
            )
            .unwrap()
            .confidence,
            DetectionConfidence::WeakHeuristic
        );
        assert_eq!(
            detect_source("cari-genetics.txt", &[], DetectedKind::GenotypeText)
                .unwrap()
                .vendor
                .as_deref(),
            Some("CariGenetics")
        );
        assert_eq!(canonicalize_ancestry_version("v2.0"), "V2.0");

        assert_eq!(
            detect_assembly("sample", &["##reference=human_g1k_v37".to_owned()]),
            Some(Assembly::Grch37)
        );
        assert_eq!(
            detect_assembly(
                "sample",
                &["##contig=<ID=chr1,length=248956422>".to_owned()]
            ),
            Some(Assembly::Grch38)
        );
        assert_eq!(detect_assembly("sample", &[]), None);
    }

    #[test]
    fn inspect_helpers_cover_index_and_render_edges() {
        let explicit = PathBuf::from("/tmp/explicit.idx");
        let options = InspectOptions {
            input_index: Some(explicit.clone()),
            ..InspectOptions::default()
        };
        assert_eq!(
            detect_index(
                Path::new("sample.txt"),
                DetectedKind::GenotypeText,
                &options
            ),
            (Some(false), Some(explicit))
        );

        let no_ext_ref = Path::new("reference");
        assert_eq!(
            detect_index(
                no_ext_ref,
                DetectedKind::ReferenceFasta,
                &InspectOptions::default()
            )
            .1,
            Some(PathBuf::from("reference.fai"))
        );
        assert_eq!(
            detect_index(
                Path::new("sample.dat"),
                DetectedKind::AlignmentCram,
                &InspectOptions::default()
            ),
            (Some(false), None)
        );

        assert_eq!(render_container(FileContainer::Plain), "plain");
        assert_eq!(render_container(FileContainer::Zip), "zip");
        assert_eq!(render_kind(DetectedKind::AlignmentBam), "alignment_bam");
        assert_eq!(render_kind(DetectedKind::Unknown), "unknown");
        assert_eq!(render_confidence(DetectionConfidence::Unknown), "unknown");
        assert_eq!(render_assembly(None), "");
        assert_eq!(render_bool(Some(true)), "true");
        assert_eq!(render_bool(Some(false)), "false");
        assert_eq!(render_bool(None), "");
    }

    #[test]
    fn inspect_helpers_cover_bgzip_zip_and_index_edges() {
        let mut bgzf_writer = bgzf::io::Writer::new(Vec::new());
        bgzf_writer
            .write_all(
                b"##fileformat=VCFv4.3\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                  chr1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0|1\n",
            )
            .unwrap();
        let bgzf_vcf = bgzf_writer.finish().unwrap();

        let bgzip_inspection =
            inspect_bytes("sample.vcf.gz", &bgzf_vcf, &InspectOptions::default()).unwrap();
        assert_eq!(bgzip_inspection.detected_kind, DetectedKind::Vcf);
        assert_eq!(bgzip_inspection.phased, Some(true));

        let cursor = Cursor::new(Vec::new());
        let mut zip_writer = zip::ZipWriter::new(cursor);
        zip_writer
            .add_directory("nested/", zip::write::SimpleFileOptions::default())
            .unwrap();
        zip_writer
            .start_file(
                "nested/sample.vcf.gz",
                zip::write::SimpleFileOptions::default(),
            )
            .unwrap();
        zip_writer.write_all(&bgzf_vcf).unwrap();
        let zip_bytes = zip_writer.finish().unwrap().into_inner();
        let zip_inspection =
            inspect_bytes("archive.zip", &zip_bytes, &InspectOptions::default()).unwrap();
        assert_eq!(zip_inspection.container, FileContainer::Zip);
        assert_eq!(zip_inspection.detected_kind, DetectedKind::Vcf);
        assert_eq!(
            zip_inspection.selected_entry.as_deref(),
            Some("nested/sample.vcf.gz")
        );

        let missing = read_zip_sample_lines_from_bytes(&zip_bytes, "missing.vcf").unwrap_err();
        assert!(missing.to_string().contains("failed to open zip entry"));

        let dir =
            std::env::temp_dir().join(format!("bioscript-inspect-unit-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let cram_no_ext = dir.join("sample.dat");
        let bam_short = dir.join("reads.bam");
        let short_bai = dir.join("reads.bai");
        std::fs::write(&cram_no_ext, b"cram").unwrap();
        std::fs::write(&bam_short, b"bam").unwrap();
        std::fs::write(&short_bai, b"bai").unwrap();

        assert_eq!(
            detect_index(
                &cram_no_ext,
                DetectedKind::AlignmentCram,
                &InspectOptions::default()
            ),
            (Some(false), None)
        );
        assert_eq!(
            detect_index(
                &bam_short,
                DetectedKind::AlignmentBam,
                &InspectOptions::default()
            ),
            (Some(true), Some(short_bai))
        );
        assert_eq!(
            classify_confidence(DetectedKind::Vcf, &[], None),
            DetectionConfidence::StrongHeuristic
        );
    }

    #[test]
    fn inspect_helpers_cover_unheaded_text_zip_fallbacks_and_render_edges() {
        let unheaded = inspect_bytes(
            "sample.txt",
            b"rs123\t1\t12345\tAG\n",
            &InspectOptions::default(),
        )
        .unwrap();
        assert_eq!(unheaded.detected_kind, DetectedKind::GenotypeText);
        assert!(
            unheaded
                .evidence
                .iter()
                .any(|line| line == "genotype-like sampled rows")
        );

        let dir =
            std::env::temp_dir().join(format!("bioscript-inspect-more-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let unknown_path = dir.join("unknown.dat");
        std::fs::write(&unknown_path, b"not enough structure\n").unwrap();
        let unknown = inspect_file(&unknown_path, &InspectOptions::default()).unwrap();
        assert_eq!(unknown.detected_kind, DetectedKind::Unknown);
        assert!(unknown.warnings[0].contains("known textual heuristics"));

        let mut bgzf_writer = bgzf::io::Writer::new(Vec::new());
        bgzf_writer
            .write_all(
                b"##fileformat=VCFv4.3\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                  chr1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            )
            .unwrap();
        let bgzf_vcf = bgzf_writer.finish().unwrap();
        let vcf_gz_path = dir.join("sample.vcf.gz");
        std::fs::write(&vcf_gz_path, &bgzf_vcf).unwrap();
        assert_eq!(read_plain_sample_lines(&vcf_gz_path).unwrap().len(), 3);

        let zip_path = dir.join("fallback.zip");
        let cursor = Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .add_directory("__MACOSX/", zip::write::SimpleFileOptions::default())
            .unwrap();
        writer
            .start_file("notes.bin", zip::write::SimpleFileOptions::default())
            .unwrap();
        writer.write_all(b"fallback bytes\n").unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        std::fs::write(&zip_path, &bytes).unwrap();
        assert_eq!(select_zip_entry(&zip_path).unwrap(), "notes.bin");

        let zip_gz_path = dir.join("vcf-gz.zip");
        let cursor = Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .start_file(
                "nested/sample.vcf.gz",
                zip::write::SimpleFileOptions::default(),
            )
            .unwrap();
        writer.write_all(&bgzf_vcf).unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        std::fs::write(&zip_gz_path, &bytes).unwrap();
        assert_eq!(
            read_zip_sample_lines(&zip_gz_path, "nested/sample.vcf.gz")
                .unwrap()
                .len(),
            3
        );

        let empty_zip_path = dir.join("empty.zip");
        let cursor = Cursor::new(Vec::new());
        let writer = zip::ZipWriter::new(cursor);
        let bytes = writer.finish().unwrap().into_inner();
        std::fs::write(&empty_zip_path, bytes).unwrap();
        let err = select_zip_entry(&empty_zip_path).unwrap_err();
        assert!(
            err.to_string()
                .contains("does not contain a supported file")
        );

        let source = detect_source(
            "dynamicdna.txt",
            &["# Dynamic DNA GSAv3 report".to_owned()],
            DetectedKind::GenotypeText,
        )
        .unwrap();
        assert_eq!(source.platform_version.as_deref(), Some("GSAv3"));
        assert_eq!(canonicalize_ancestry_version("v2"), "V2");
        assert_eq!(render_kind(DetectedKind::AlignmentCram), "alignment_cram");
        assert_eq!(render_kind(DetectedKind::AlignmentBam), "alignment_bam");
        assert_eq!(render_assembly(Some(Assembly::Grch38)), "grch38");
        assert_eq!(render_bool(Some(true)), "true");
    }
}
