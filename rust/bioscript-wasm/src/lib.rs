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

mod js_reader;

use std::{io::BufReader, path::PathBuf};

use bioscript_core::GenomicLocus;
use bioscript_formats::{
    DetectedKind, DetectionConfidence, FileContainer, FileInspection, InspectOptions,
    SourceMetadata, alignment, inspect_bytes as inspect_bytes_rs, observe_cram_snp_with_reader,
    observe_vcf_snp_with_reader,
};
use noodles::csi as noodles_csi;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::js_reader::JsReader;

#[wasm_bindgen(start)]
pub fn start() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

#[derive(Default, Deserialize)]
struct InspectOptionsJs {
    input_index: Option<String>,
    reference_file: Option<String>,
    reference_index: Option<String>,
}

/// Classify bytes as a known genomic file. Mirrors `bioscript-formats::inspect::inspect_bytes`.
/// Returns JSON matching the `Inspection` shape the app already uses.
#[wasm_bindgen(js_name = inspectBytes)]
pub fn inspect_bytes(name: &str, bytes: &[u8], options_json: Option<String>) -> Result<String, JsError> {
    let options_js: InspectOptionsJs = match options_json {
        Some(text) if !text.is_empty() => serde_json::from_str(&text)
            .map_err(|err| JsError::new(&format!("invalid InspectOptions JSON: {err}")))?,
        _ => InspectOptionsJs::default(),
    };
    let options = InspectOptions {
        input_index: options_js.input_index.map(PathBuf::from),
        reference_file: options_js.reference_file.map(PathBuf::from),
        reference_index: options_js.reference_index.map(PathBuf::from),
    };

    let inspection = inspect_bytes_rs(name, bytes, &options)
        .map_err(|err| JsError::new(&format!("inspect_bytes failed: {err:?}")))?;

    let resp = InspectionJs::from(inspection);
    serde_json::to_string(&resp)
        .map_err(|err| JsError::new(&format!("failed to encode response: {err}")))
}

#[derive(Deserialize)]
struct VariantInput {
    name: String,
    chrom: String,
    // 1-based genomic position of the SNP.
    pos: i64,
    #[serde(rename = "ref")]
    ref_base: String,
    #[serde(rename = "alt")]
    alt_base: String,
    #[serde(default)]
    rsid: Option<String>,
    #[serde(default)]
    assembly: Option<String>,
}

#[derive(Serialize)]
struct VariantObservationJs {
    name: String,
    backend: String,
    #[serde(rename = "matchedRsid", skip_serializing_if = "Option::is_none")]
    matched_rsid: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    genotype: Option<String>,
    #[serde(rename = "refCount", skip_serializing_if = "Option::is_none")]
    ref_count: Option<u32>,
    #[serde(rename = "altCount", skip_serializing_if = "Option::is_none")]
    alt_count: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    depth: Option<u32>,
    #[serde(rename = "rawCounts")]
    raw_counts: std::collections::BTreeMap<String, u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    decision: Option<String>,
    evidence: Vec<String>,
}

/// Observe a list of SNP variants against an indexed CRAM + reference FASTA,
/// with the bulk bytes pulled on demand via JS-supplied `readAt(offset, len)`
/// callbacks. The small index payloads (`.crai`, `.fai`) are passed inline.
///
/// Both callbacks must return a `Uint8Array` synchronously (or via a Node
/// sync read) — wasm's `Read + Seek` contract is synchronous. Async reads are
/// a follow-up that needs buffered pre-fetch on the JS side.
#[wasm_bindgen(js_name = lookupCramVariants)]
pub fn lookup_cram_variants(
    cram_read_at: js_sys::Function,
    cram_len: f64,
    crai_bytes: &[u8],
    fasta_read_at: js_sys::Function,
    fasta_len: f64,
    fai_bytes: &[u8],
    variants_json: &str,
) -> Result<String, JsError> {
    let crai_index = alignment::parse_crai_bytes(crai_bytes)
        .map_err(|err| JsError::new(&format!("parse crai: {err:?}")))?;
    let fai_index = alignment::parse_fai_bytes(fai_bytes)
        .map_err(|err| JsError::new(&format!("parse fai: {err:?}")))?;

    let fasta_reader = BufReader::new(JsReader::new(fasta_read_at, fasta_len as u64, "fasta"));
    let repository = alignment::build_reference_repository_from_readers(fasta_reader, fai_index);

    let cram_reader = JsReader::new(cram_read_at, cram_len as u64, "cram");
    let mut indexed = alignment::build_cram_indexed_reader_from_reader(
        cram_reader,
        crai_index,
        repository,
    )
    .map_err(|err| JsError::new(&format!("build cram reader: {err:?}")))?;

    let variants: Vec<VariantInput> = serde_json::from_str(variants_json)
        .map_err(|err| JsError::new(&format!("parse variantsJson: {err}")))?;

    let mut results = Vec::with_capacity(variants.len());
    for variant in variants {
        let ref_char = variant
            .ref_base
            .chars()
            .next()
            .ok_or_else(|| JsError::new(&format!("variant {}: empty ref", variant.name)))?;
        let alt_char = variant
            .alt_base
            .chars()
            .next()
            .ok_or_else(|| JsError::new(&format!("variant {}: empty alt", variant.name)))?;
        let assembly = variant
            .assembly
            .as_deref()
            .and_then(parse_assembly_str);
        let locus = GenomicLocus {
            chrom: variant.chrom.clone(),
            start: variant.pos,
            end: variant.pos,
        };
        let observation = observe_cram_snp_with_reader(
            &mut indexed,
            &variant.name,
            &locus,
            ref_char,
            alt_char,
            variant.rsid.clone(),
            assembly,
        )
        .map_err(|err| JsError::new(&format!("lookup {}: {err:?}", variant.name)))?;
        results.push(VariantObservationJs {
            name: variant.name,
            backend: observation.backend,
            matched_rsid: observation.matched_rsid,
            assembly: observation.assembly.map(|a| render_assembly(a).to_owned()),
            genotype: observation.genotype,
            ref_count: observation.ref_count,
            alt_count: observation.alt_count,
            depth: observation.depth,
            raw_counts: observation.raw_counts,
            decision: observation.decision,
            evidence: observation.evidence,
        });
    }

    serde_json::to_string(&results).map_err(|err| JsError::new(&format!("encode results: {err}")))
}

/// Observe a list of SNP variants against a bgzipped + tabix-indexed VCF,
/// with the bulk bytes pulled on demand via a JS-supplied `readAt(offset, len)`
/// callback. The small `.tbi` payload is passed inline.
///
/// The reader must provide the VCF synchronously — on web this is a
/// `FileReaderSync`-backed callback running inside a Web Worker.
#[wasm_bindgen(js_name = lookupVcfVariants)]
pub fn lookup_vcf_variants(
    vcf_read_at: js_sys::Function,
    vcf_len: f64,
    tbi_bytes: &[u8],
    variants_json: &str,
) -> Result<String, JsError> {
    let tabix_index = alignment::parse_tbi_bytes(tbi_bytes)
        .map_err(|err| JsError::new(&format!("parse tbi: {err:?}")))?;
    let vcf_reader = JsReader::new(vcf_read_at, vcf_len as u64, "vcf");
    let mut indexed = noodles_csi::io::IndexedReader::new(vcf_reader, tabix_index);

    let variants: Vec<VariantInput> = serde_json::from_str(variants_json)
        .map_err(|err| JsError::new(&format!("parse variantsJson: {err}")))?;

    let mut results = Vec::with_capacity(variants.len());
    for variant in variants {
        let ref_char = variant
            .ref_base
            .chars()
            .next()
            .ok_or_else(|| JsError::new(&format!("variant {}: empty ref", variant.name)))?;
        let alt_char = variant
            .alt_base
            .chars()
            .next()
            .ok_or_else(|| JsError::new(&format!("variant {}: empty alt", variant.name)))?;
        let assembly = variant
            .assembly
            .as_deref()
            .and_then(parse_assembly_str);
        let locus = GenomicLocus {
            chrom: variant.chrom.clone(),
            start: variant.pos,
            end: variant.pos,
        };
        let observation = observe_vcf_snp_with_reader(
            &mut indexed,
            &variant.name,
            &locus,
            ref_char,
            alt_char,
            variant.rsid.clone(),
            assembly,
        )
        .map_err(|err| JsError::new(&format!("vcf lookup {}: {err:?}", variant.name)))?;
        results.push(VariantObservationJs {
            name: variant.name,
            backend: observation.backend,
            matched_rsid: observation.matched_rsid,
            assembly: observation.assembly.map(|a| render_assembly(a).to_owned()),
            genotype: observation.genotype,
            ref_count: observation.ref_count,
            alt_count: observation.alt_count,
            depth: observation.depth,
            raw_counts: observation.raw_counts,
            decision: observation.decision,
            evidence: observation.evidence,
        });
    }

    serde_json::to_string(&results).map_err(|err| JsError::new(&format!("encode results: {err}")))
}

fn parse_assembly_str(s: &str) -> Option<bioscript_core::Assembly> {
    match s.to_ascii_lowercase().as_str() {
        "grch37" | "hg19" | "b37" => Some(bioscript_core::Assembly::Grch37),
        "grch38" | "hg38" => Some(bioscript_core::Assembly::Grch38),
        _ => None,
    }
}

// Wire types — we flatten the Rust FileInspection into the shape the app's
// TS Inspection type already expects (matches widgets/FilePicker/types.ts).
#[derive(Serialize)]
struct InspectionJs {
    #[serde(rename = "fileName")]
    file_name: String,
    container: &'static str,
    #[serde(rename = "detectedKind")]
    detected_kind: &'static str,
    confidence: &'static str,
    assembly: Option<&'static str>,
    phased: Option<bool>,
    source: Option<SourceJs>,
    #[serde(rename = "selectedEntry", skip_serializing_if = "Option::is_none")]
    selected_entry: Option<String>,
    #[serde(rename = "hasIndex", skip_serializing_if = "Option::is_none")]
    has_index: Option<bool>,
    #[serde(rename = "referenceMatches", skip_serializing_if = "Option::is_none")]
    reference_matches: Option<bool>,
    evidence: Vec<String>,
    warnings: Vec<String>,
    #[serde(rename = "durationMs")]
    duration_ms: u128,
}

#[derive(Serialize)]
struct SourceJs {
    vendor: String,
    #[serde(rename = "platformVersion", skip_serializing_if = "Option::is_none")]
    platform_version: Option<String>,
    confidence: &'static str,
    evidence: Vec<String>,
}

impl From<FileInspection> for InspectionJs {
    fn from(i: FileInspection) -> Self {
        InspectionJs {
            file_name: i.path.display().to_string(),
            container: render_container(i.container),
            detected_kind: render_kind(i.detected_kind),
            confidence: render_confidence(i.confidence),
            assembly: i.assembly.map(render_assembly),
            phased: i.phased,
            source: i.source.map(SourceJs::from),
            selected_entry: i.selected_entry,
            has_index: i.has_index,
            reference_matches: i.reference_matches,
            evidence: i.evidence,
            warnings: i.warnings,
            duration_ms: i.duration_ms,
        }
    }
}

impl From<SourceMetadata> for SourceJs {
    fn from(s: SourceMetadata) -> Self {
        SourceJs {
            vendor: s.vendor.unwrap_or_default(),
            platform_version: s.platform_version,
            confidence: render_confidence(s.confidence),
            evidence: s.evidence,
        }
    }
}

fn render_container(c: FileContainer) -> &'static str {
    match c {
        FileContainer::Plain => "plain",
        FileContainer::Zip => "zip",
    }
}

fn render_kind(k: DetectedKind) -> &'static str {
    match k {
        DetectedKind::GenotypeText => "genotype_text",
        DetectedKind::Vcf => "vcf",
        DetectedKind::AlignmentCram => "alignment_cram",
        DetectedKind::AlignmentBam => "alignment_bam",
        DetectedKind::ReferenceFasta => "reference_fasta",
        DetectedKind::Unknown => "unknown",
    }
}

fn render_confidence(c: DetectionConfidence) -> &'static str {
    match c {
        DetectionConfidence::Authoritative => "authoritative",
        DetectionConfidence::StrongHeuristic => "strong_heuristic",
        DetectionConfidence::WeakHeuristic => "weak_heuristic",
        DetectionConfidence::Unknown => "unknown",
    }
}

fn render_assembly(a: bioscript_core::Assembly) -> &'static str {
    match a {
        bioscript_core::Assembly::Grch37 => "grch37",
        bioscript_core::Assembly::Grch38 => "grch38",
    }
}
