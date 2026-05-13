use wasm_bindgen::{JsError, prelude::wasm_bindgen};

use crate::js_reader::JsReader;

#[wasm_bindgen(js_name = generateVcfTbi)]
pub fn generate_vcf_tbi(input_name: &str, vcf_bytes: &[u8]) -> Result<Vec<u8>, JsError> {
    bioscript_formats::alignment::generate_vcf_tbi_bytes(vcf_bytes)
        .map_err(|err| JsError::new(&format!("generate tabix index for {input_name}: {err:?}")))
}

#[wasm_bindgen(js_name = generateCramCrai)]
pub fn generate_cram_crai(input_name: &str, cram_bytes: &[u8]) -> Result<Vec<u8>, JsError> {
    bioscript_formats::alignment::generate_cram_crai_bytes(cram_bytes)
        .map_err(|err| JsError::new(&format!("generate CRAI index for {input_name}: {err:?}")))
}

#[wasm_bindgen(js_name = generateCramCraiFromReader)]
pub fn generate_cram_crai_from_reader(
    input_name: &str,
    cram_read_at: js_sys::Function,
    cram_len: f64,
) -> Result<Vec<u8>, JsError> {
    let reader = JsReader::new(cram_read_at, cram_len as u64, "cram-index");
    bioscript_formats::alignment::generate_cram_crai_reader(reader)
        .map_err(|err| JsError::new(&format!("generate CRAI index for {input_name}: {err:?}")))
}

#[wasm_bindgen(js_name = generateBamBai)]
pub fn generate_bam_bai(input_name: &str, bam_bytes: &[u8]) -> Result<Vec<u8>, JsError> {
    bioscript_formats::alignment::generate_bam_bai_bytes(bam_bytes)
        .map_err(|err| JsError::new(&format!("generate BAI index for {input_name}: {err:?}")))
}

#[wasm_bindgen(js_name = generateBamBaiFromReader)]
pub fn generate_bam_bai_from_reader(
    input_name: &str,
    bam_read_at: js_sys::Function,
    bam_len: f64,
) -> Result<Vec<u8>, JsError> {
    let reader = JsReader::new(bam_read_at, bam_len as u64, "bam-index");
    bioscript_formats::alignment::generate_bam_bai_reader(reader)
        .map_err(|err| JsError::new(&format!("generate BAI index for {input_name}: {err:?}")))
}

#[wasm_bindgen(js_name = generateFastaFai)]
pub fn generate_fasta_fai(input_name: &str, fasta_bytes: &[u8]) -> Result<Vec<u8>, JsError> {
    bioscript_formats::alignment::generate_fasta_fai_bytes(fasta_bytes)
        .map_err(|err| JsError::new(&format!("generate FAI index for {input_name}: {err:?}")))
}

#[wasm_bindgen(js_name = generateFastaFaiFromReader)]
pub fn generate_fasta_fai_from_reader(
    input_name: &str,
    fasta_read_at: js_sys::Function,
    fasta_len: f64,
) -> Result<Vec<u8>, JsError> {
    let reader = JsReader::new(fasta_read_at, fasta_len as u64, "fasta-index");
    bioscript_formats::alignment::generate_fasta_fai_reader(reader)
        .map_err(|err| JsError::new(&format!("generate FAI index for {input_name}: {err:?}")))
}
