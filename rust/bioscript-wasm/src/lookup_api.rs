use std::io::BufReader;

use bioscript_core::{GenomicLocus, VariantKind, VariantObservation, VariantSpec};
use bioscript_formats::{
    alignment, observe_cram_indel_with_reader, observe_cram_snp_with_reader,
    observe_vcf_snp_with_reader, GenotypeStore,
};
use noodles::csi as noodles_csi;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::{inspect_api::render_assembly, js_reader::JsReader};

#[derive(Deserialize)]
struct VariantInput {
    name: String,
    chrom: String,
    // 1-based genomic interval. `pos` is accepted for older callers.
    #[serde(default)]
    pos: Option<i64>,
    #[serde(default)]
    start: Option<i64>,
    #[serde(default)]
    end: Option<i64>,
    #[serde(rename = "ref")]
    ref_base: String,
    #[serde(rename = "alt")]
    alt_base: String,
    #[serde(default)]
    rsid: Option<String>,
    #[serde(default)]
    assembly: Option<String>,
    #[serde(default)]
    kind: Option<String>,
}

#[derive(Serialize)]
struct VariantObservationJs {
    name: String,
    backend: String,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    reference: Option<String>,
    #[serde(rename = "alt", skip_serializing_if = "Option::is_none")]
    alternate: Option<String>,
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
    let mut indexed =
        alignment::build_cram_indexed_reader_from_reader(cram_reader, crai_index, repository)
            .map_err(|err| JsError::new(&format!("build cram reader: {err:?}")))?;

    let variants = parse_variants_json(variants_json)?;
    let mut results = Vec::with_capacity(variants.len());
    for variant in variants {
        let assembly = variant.assembly.as_deref().and_then(parse_assembly_str);
        let locus = variant_locus(&variant)?;
        let kind = parse_variant_kind(variant.kind.as_deref()).unwrap_or(VariantKind::Snp);
        let observation = match kind {
            VariantKind::Snp => {
                ensure_single_base_variant(&variant)?;
                let ref_char = first_allele_char(&variant.name, &variant.ref_base, "ref")?;
                let alt_char = first_allele_char(&variant.name, &variant.alt_base, "alt")?;
                observe_cram_snp_with_reader(
                    &mut indexed,
                    &variant.name,
                    &locus,
                    ref_char,
                    alt_char,
                    variant.rsid.clone(),
                    assembly,
                )
            }
            VariantKind::Insertion | VariantKind::Indel => observe_cram_indel_with_reader(
                &mut indexed,
                &variant.name,
                &locus,
                &variant.ref_base,
                &variant.alt_base,
                variant.rsid.clone(),
                assembly,
            ),
            other => {
                return Err(JsError::new(&format!(
                    "variant {} has unsupported kind {:?} for web CRAM lookup",
                    variant.name, other
                )));
            }
        }
        .map_err(|err| JsError::new(&format!("lookup {}: {err:?}", variant.name)))?;
        results.push(observation_to_js(variant, observation));
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

    let variants = parse_variants_json(variants_json)?;
    let mut results = Vec::with_capacity(variants.len());
    for variant in variants {
        ensure_single_base_variant(&variant)?;
        let ref_char = first_allele_char(&variant.name, &variant.ref_base, "ref")?;
        let alt_char = first_allele_char(&variant.name, &variant.alt_base, "alt")?;
        let assembly = variant.assembly.as_deref().and_then(parse_assembly_str);
        let locus = variant_locus(&variant)?;
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
        results.push(observation_to_js(variant, observation));
    }

    serde_json::to_string(&results).map_err(|err| JsError::new(&format!("encode results: {err}")))
}

#[wasm_bindgen(js_name = lookupGenotypeBytesVariants)]
pub fn lookup_genotype_bytes_variants(
    name: &str,
    bytes: &[u8],
    variants_json: &str,
) -> Result<String, JsError> {
    let store = GenotypeStore::from_bytes(name, bytes)
        .map_err(|err| JsError::new(&format!("load genotype bytes {name}: {err:?}")))?;
    let variants = parse_variants_json(variants_json)?;
    let specs = variants
        .iter()
        .map(variant_input_to_spec)
        .collect::<Result<Vec<_>, _>>()?;
    let observations = store
        .lookup_variants(&specs)
        .map_err(|err| JsError::new(&format!("lookup genotype bytes {name}: {err:?}")))?;
    let rows = variants
        .into_iter()
        .zip(observations)
        .map(|(variant, observation)| observation_to_js(variant, observation))
        .collect::<Vec<_>>();
    serde_json::to_string(&rows).map_err(|err| JsError::new(&format!("encode results: {err}")))
}

#[wasm_bindgen(js_name = lookupGenotypeBytesRsids)]
pub fn lookup_genotype_bytes_rsids(
    name: &str,
    bytes: &[u8],
    rsids_json: &str,
) -> Result<String, JsError> {
    let store = GenotypeStore::from_bytes(name, bytes)
        .map_err(|err| JsError::new(&format!("load genotype bytes {name}: {err:?}")))?;
    let rsids: Vec<String> = serde_json::from_str(rsids_json)
        .map_err(|err| JsError::new(&format!("parse rsidsJson: {err}")))?;
    let values = rsids
        .iter()
        .map(|rsid| {
            store
                .get(rsid)
                .map_err(|err| JsError::new(&format!("lookup genotype rsid {rsid}: {err:?}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    serde_json::to_string(&values).map_err(|err| JsError::new(&format!("encode results: {err}")))
}

fn parse_variants_json(variants_json: &str) -> Result<Vec<VariantInput>, JsError> {
    serde_json::from_str(variants_json)
        .map_err(|err| JsError::new(&format!("parse variantsJson: {err}")))
}

fn variant_locus(variant: &VariantInput) -> Result<GenomicLocus, JsError> {
    let start = variant
        .start
        .or(variant.pos)
        .ok_or_else(|| JsError::new(&format!("variant {}: start/pos missing", variant.name)))?;
    let end = variant.end.unwrap_or(start);
    Ok(GenomicLocus {
        chrom: variant.chrom.clone(),
        start,
        end,
    })
}

fn first_allele_char(variant_name: &str, allele: &str, label: &str) -> Result<char, JsError> {
    allele
        .chars()
        .next()
        .ok_or_else(|| JsError::new(&format!("variant {variant_name}: empty {label}")))
}

fn ensure_single_base_variant(variant: &VariantInput) -> Result<(), JsError> {
    let kind = variant
        .kind
        .as_deref()
        .unwrap_or("snv")
        .to_ascii_lowercase();
    let is_snp_kind = matches!(kind.as_str(), "snp" | "snv" | "variant" | "");
    if !is_snp_kind
        || variant.ref_base.chars().count() != 1
        || variant.alt_base.chars().count() != 1
    {
        return Err(JsError::new(&format!(
            "variant {} has kind/ref/alt {} {}/{}; web CRAM/VCF lookup currently supports single-base SNV observations only",
            variant.name,
            variant.kind.as_deref().unwrap_or("snv"),
            variant.ref_base,
            variant.alt_base
        )));
    }
    Ok(())
}

fn variant_input_to_spec(variant: &VariantInput) -> Result<VariantSpec, JsError> {
    let locus = variant_locus(variant)?;
    let assembly = variant.assembly.as_deref().and_then(parse_assembly_str);
    let kind = parse_variant_kind(variant.kind.as_deref());
    Ok(VariantSpec {
        rsids: variant.rsid.clone().into_iter().collect(),
        grch37: if assembly == Some(bioscript_core::Assembly::Grch37) {
            Some(locus.clone())
        } else {
            None
        },
        grch38: if assembly == Some(bioscript_core::Assembly::Grch38) || assembly.is_none() {
            Some(locus)
        } else {
            None
        },
        reference: Some(variant.ref_base.clone()),
        alternate: Some(variant.alt_base.clone()),
        kind,
        deletion_length: None,
        motifs: Vec::new(),
    })
}

fn parse_variant_kind(kind: Option<&str>) -> Option<VariantKind> {
    match kind.unwrap_or("").to_ascii_lowercase().as_str() {
        "snp" | "snv" | "variant" | "" => Some(VariantKind::Snp),
        "insertion" => Some(VariantKind::Insertion),
        "deletion" => Some(VariantKind::Deletion),
        "indel" => Some(VariantKind::Indel),
        "other" => Some(VariantKind::Other),
        _ => Some(VariantKind::Other),
    }
}

fn observation_to_js(
    variant: VariantInput,
    observation: VariantObservation,
) -> VariantObservationJs {
    VariantObservationJs {
        name: variant.name,
        backend: observation.backend,
        reference: Some(variant.ref_base),
        alternate: Some(variant.alt_base),
        matched_rsid: observation.matched_rsid,
        assembly: observation.assembly.map(|a| render_assembly(a).to_owned()),
        genotype: observation.genotype,
        ref_count: observation.ref_count,
        alt_count: observation.alt_count,
        depth: observation.depth,
        raw_counts: observation.raw_counts,
        decision: observation.decision,
        evidence: observation.evidence,
    }
}

fn parse_assembly_str(s: &str) -> Option<bioscript_core::Assembly> {
    match s.to_ascii_lowercase().as_str() {
        "grch37" | "hg19" | "b37" => Some(bioscript_core::Assembly::Grch37),
        "grch38" | "hg38" => Some(bioscript_core::Assembly::Grch38),
        _ => None,
    }
}
