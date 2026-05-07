struct AppObservationJson {
    allele_balance: Option<f64>,
    alt_count: Option<u32>,
    assay_id: String,
    assembly: String,
    call: ObservationCallValues,
    chrom: String,
    depth: Option<u32>,
    evidence_raw: String,
    gene: String,
    genotype: String,
    genotype_display: String,
    kind: String,
    locus: Option<bioscript_core::GenomicLocus>,
    manifest: bioscript_schema::VariantManifest,
    non_reportable_status: Option<&'static str>,
    observed_alt_alleles: Vec<String>,
    ref_allele: String,
    ref_count: Option<u32>,
    reportable_alt: String,
    row: BTreeMap<String, String>,
    row_path: String,
    source: serde_json::Value,
    weak_indel_match: bool,
    zygosity: String,
}

struct ObservationCallValues {
    outcome: &'static str,
    status: &'static str,
    reported_genotype_display: String,
}

fn observation_call_values(
    depth: Option<u32>,
    non_reportable_status: Option<&'static str>,
    genotype: &str,
    zygosity: &str,
    genotype_display: &str,
) -> ObservationCallValues {
    let outcome = if depth == Some(0) {
        "not_covered"
    } else if non_reportable_status == Some("observed_alt") {
        "observed_alt"
    } else if non_reportable_status == Some("unknown_alt") {
        "unknown_alt"
    } else if genotype == "./." {
        "no_call"
    } else if zygosity == "hom_ref" || zygosity == "hem_ref" {
        "reference"
    } else if zygosity == "het" || zygosity == "hom_alt" || zygosity == "hem_alt" {
        "variant"
    } else {
        "unknown"
    };
    let status = if matches!(outcome, "observed_alt" | "unknown_alt") {
        outcome
    } else if genotype == "./." {
        "no_call"
    } else {
        "called"
    };
    let reported_genotype_display = if matches!(zygosity, "hem_ref" | "hem_alt") {
        hemizygous_display_genotype(genotype_display)
    } else if genotype_display.is_empty() && matches!(outcome, "no_call" | "not_covered") {
        "??".to_owned()
    } else {
        genotype_display.to_owned()
    };
    ObservationCallValues {
        outcome,
        status,
        reported_genotype_display,
    }
}

fn render_app_observation_json(input: AppObservationJson) -> serde_json::Value {
    let AppObservationJson {
        allele_balance,
        alt_count,
        assay_id,
        assembly,
        call,
        chrom,
        depth,
        evidence_raw,
        gene,
        genotype,
        genotype_display,
        kind,
        locus,
        manifest,
        non_reportable_status,
        observed_alt_alleles,
        ref_allele,
        ref_count,
        reportable_alt,
        row,
        row_path,
        source,
        weak_indel_match,
        zygosity,
    } = input;
    serde_json::json!({
        "participant_id": row.get("participant_id").cloned().unwrap_or_default(),
        "assay_id": assay_id,
        "assay_version": "1.0",
        "variant_key": manifest.name,
        "variant_path": row_path,
        "rsid": row.get("matched_rsid").filter(|value| !value.is_empty()).cloned().or_else(|| manifest.spec.rsids.first().cloned()),
        "gene": gene,
        "assembly": if assembly.is_empty() { serde_json::Value::Null } else { serde_json::Value::String(assembly.to_uppercase()) },
        "chrom": chrom,
        "pos_start": locus.as_ref().map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.start)),
        "pos_end": locus.as_ref().map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.end)),
        "ref": ref_allele,
        "alt": reportable_alt,
        "kind": kind,
        "match_status": if row.get("matched_rsid").is_some_and(|value| !value.is_empty()) || !genotype_display.is_empty() { "found" } else { "not_found" },
        "coverage_status": depth.map_or("covered", |depth| if depth > 0 { "covered" } else { "not_covered" }),
        "call_status": call.status,
        "genotype": genotype,
        "genotype_display": call.reported_genotype_display,
        "zygosity": zygosity,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "depth": depth,
        "genotype_quality": serde_json::Value::Null,
        "allele_balance": allele_balance,
        "outcome": call.outcome,
        "evidence_type": if row.get("backend").is_some_and(|value| value == "cram") { "mpileup" } else { "genotype_file" },
        "evidence_raw": evidence_raw,
        "source": source,
        "match_quality": if weak_indel_match { serde_json::Value::String("weak".to_owned()) } else { serde_json::Value::Null },
        "match_notes": if weak_indel_match {
            serde_json::Value::String("consumer genotype file reported an insertion/deletion token at the marker, not sequence-resolved evidence for the exact deletion allele".to_owned())
        } else {
            serde_json::Value::Null
        },
        "facets": observation_facets(non_reportable_status, &observed_alt_alleles),
    })
}
