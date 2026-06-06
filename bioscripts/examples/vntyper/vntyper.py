# MUC1 VNTR (VNtyper) advanced assay analysis.
#
# Unlike a SNP-lookup assay, this does not use bioscript.query_plan. It treats
# the dragged-in genome (input_file) as an aligned BAM/CRAM and runs the native
# samtools -> kestrel -> bcftools pipeline that was built on the
# madhava/libs branch, then VNtyper post-processing via vcf.read_vntyper_kestrel.
#
# Runtime globals provided by the assay runner (same input/output contract as
# any assay): input_file, output_file, participant_id, asset_paths. Alignment
# package-report runners should also provide input_index for BAM/CRAM indexes
# and alignment_reference_fasta/reference_fasta plus
# alignment_reference_index/reference_index for CRAM region slicing.

from bioscript import kestrel
from bioscript import samtools
from bioscript import vcf

# MUC1 region for GRCh37/hg19 aligned inputs (the assembly the VNtyper test
# fixtures use). bam_region is the slice fed to Kestrel; vntr_region is the
# tighter span used for depth context.
MUC1_BAM_REGION = "chr1:155158000-155163000"
MUC1_VNTR_REGION = "chr1:155160500-155162000"

# VNtyper-correct Kestrel configuration. These mirror the values the proven
# native VNtyper pipeline uses (ports/vntyper/bioscript/vntyper_config.py):
# the MUC1 VNTR call only resolves with min_kmer_count=5 and the bounded
# haplotype/saved-state caps.
KMER_SIZE = 20
MIN_KMER_COUNT = 5
MAX_HAPLOTYPES = 2
MAX_SAVED_STATES = 2


def best_passing_row(rows):
    best = None
    best_score = -1.0
    for row in rows:
        if str(row.get("passes_vntyper_filters")) != "True":
            continue
        score = 0.0
        raw = row.get("Depth_Score")
        if raw is not None and raw != "" and raw != "None":
            score = float(raw)
        if score > best_score:
            best_score = score
            best = row
    return best


def text_value(value):
    if value is None:
        return ""
    return str(value)


def html_escape(value):
    text = text_value(value)
    text = text.replace("&", "&amp;")
    text = text.replace("<", "&lt;")
    text = text.replace(">", "&gt;")
    text = text.replace('"', "&quot;")
    return text


def json_escape(value):
    text = text_value(value)
    text = text.replace("\\", "\\\\")
    text = text.replace('"', '\\"')
    text = text.replace("\n", "\\n")
    return text


def write_json_text(path, entries):
    parts = ["{"]
    first = True
    for key in entries:
        if not first:
            parts.append(",")
        first = False
        parts.append('\n  "' + json_escape(key) + '": "' + json_escape(entries[key]) + '"')
    parts.append("\n}\n")
    bioscript.write_text(path, "".join(parts))


def write_summary_html(path, sample, status, outcome, confidence, variant, alt_depth, notes, rows):
    body = [
        "<!doctype html><html><head><meta charset=\"utf-8\">",
        "<title>VNtyper MUC1 report</title>",
        "<style>",
        "body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;margin:32px;color:#111827}",
        "h1{font-size:28px;margin:0 0 8px} h2{margin-top:28px}",
        ".badge{display:inline-block;border-radius:999px;padding:4px 10px;font-weight:700}",
        ".positive{background:#fee2e2;color:#991b1b}.negative{background:#dcfce7;color:#166534}",
        "table{border-collapse:collapse;width:100%;font-size:13px}td,th{border:1px solid #d1d5db;padding:6px;text-align:left}",
        "code{background:#f3f4f6;padding:2px 4px;border-radius:4px}",
        "</style></head><body>",
        "<h1>VNtyper MUC1 report</h1>",
        "<p>Sample <code>" + html_escape(sample) + "</code></p>",
    ]
    badge_class = "positive"
    if status != "positive":
        badge_class = "negative"
    body.append("<p><span class=\"badge " + badge_class + "\">" + html_escape(status) + "</span></p>")
    body.append("<table><tbody>")
    body.append("<tr><th>Outcome</th><td>" + html_escape(outcome) + "</td></tr>")
    body.append("<tr><th>Confidence</th><td>" + html_escape(confidence) + "</td></tr>")
    body.append("<tr><th>Variant</th><td>" + html_escape(variant) + "</td></tr>")
    body.append("<tr><th>Alternate depth</th><td>" + html_escape(alt_depth) + "</td></tr>")
    body.append("<tr><th>Notes</th><td>" + html_escape(notes) + "</td></tr>")
    body.append("</tbody></table>")
    body.append("<h2>Kestrel result</h2>")
    if len(rows) == 0:
        body.append("<p>No Kestrel rows passed parsing.</p>")
    else:
        display_rows = rows
        if len(display_rows) > 20:
            display_rows = rows[:20]
        headers = []
        for key in display_rows[0]:
            headers.append(key)
        body.append("<table><thead><tr>")
        for header in headers:
            body.append("<th>" + html_escape(header) + "</th>")
        body.append("</tr></thead><tbody>")
        for row in display_rows:
            body.append("<tr>")
            for header in headers:
                body.append("<td>" + html_escape(row.get(header)) + "</td>")
            body.append("</tr>")
        body.append("</tbody></table>")
        if len(rows) > len(display_rows):
            body.append("<p>Showing " + str(len(display_rows)) + " of " + str(len(rows)) + " parsed Kestrel rows. See <code>kestrel/kestrel_result.tsv</code> for all rows.</p>")
    body.append("</body></html>\n")
    bioscript.write_text(path, "".join(body))


def main():
    sample = participant_id
    reference_fasta = asset_paths["muc1_reference"]
    bam_region = MUC1_BAM_REGION
    if "vntyper_bam_region" in bioscript.context:
        bam_region = bioscript.context["vntyper_bam_region"]
    alignment_index = None
    if "input_index" in bioscript.context:
        alignment_index = bioscript.context["input_index"]
    if alignment_index is None and "input_bai" in bioscript.context:
        alignment_index = bioscript.context["input_bai"]
    alignment_reference_fasta = None
    if "alignment_reference_fasta" in bioscript.context:
        alignment_reference_fasta = bioscript.context["alignment_reference_fasta"]
    if alignment_reference_fasta is None and "reference_fasta" in bioscript.context:
        alignment_reference_fasta = bioscript.context["reference_fasta"]
    alignment_reference_index = None
    if "alignment_reference_index" in bioscript.context:
        alignment_reference_index = bioscript.context["alignment_reference_index"]
    if alignment_reference_index is None and "reference_index" in bioscript.context:
        alignment_reference_index = bioscript.context["reference_index"]

    sliced_bam = "/work/sliced.bam"
    sliced_bai = "/work/sliced.bam.bai"
    fastq_1 = "/work/reads_R1.fastq.gz"
    fastq_2 = "/work/reads_R2.fastq.gz"
    kestrel_vcf = "/work/kestrel.vcf"

    if alignment_index is None:
        alignment_index = "/work/input.bam.bai"
        samtools.index(input_file, alignment_index)
    samtools.view_region_native(
        input_file,
        bam_region,
        sliced_bam,
        alignment_index,
        alignment_reference_fasta,
        alignment_reference_index,
    )
    samtools.fastq_native(
        input_file,
        bam_region,
        fastq_1,
        fastq_2,
        alignment_index,
        alignment_reference_fasta,
        alignment_reference_index,
    )
    samtools.index(sliced_bam, sliced_bai)

    kestrel.run_native(
        reference_fasta,
        [fastq_1, fastq_2],
        kestrel_vcf,
        kmer_size=KMER_SIZE,
        sample_name=sample,
        min_kmer_count=MIN_KMER_COUNT,
        max_haplotypes=MAX_HAPLOTYPES,
        max_saved_states=MAX_SAVED_STATES,
    )
    rows = vcf.read_vntyper_kestrel(kestrel_vcf)
    called = best_passing_row(rows)
    depth_summary = samtools.depth_native(
        input_file,
        MUC1_VNTR_REGION,
        alignment_index,
    )

    if called is None:
        outcome = "normal"
        status = "negative"
        confidence = "Negative"
        variant = "none"
        alt_depth = "0"
        notes = (
            "No MUC1 VNTR frameshift passed VNtyper filters. This is the "
            "expected (normal) result. Consult a licensed doctor for advice."
        )
    else:
        outcome = "possible frameshift"
        status = "positive"
        confidence = str(called.get("Confidence"))
        variant = (
            str(called.get("CHROM"))
            + ":"
            + str(called.get("POS"))
            + " "
            + str(called.get("REF"))
            + ">"
            + str(called.get("ALT"))
        )
        alt_depth = str(called.get("Estimated_Depth_AlternateVariant"))
        notes = (
            'A possible MUC1 VNTR frameshift was detected ("'
            + confidence
            + '" confidence, variant '
            + variant
            + "). This is consistent with ADTKD-MUC1 and should be confirmed "
            + "with an orthogonal method (e.g. SNaPshot or long-read "
            + "sequencing). Consult a licensed doctor for advice."
        )

    rows_out = [
        {
            "participant_id": participant_id,
            "vntyper_outcome": outcome,
            "vntyper_status": status,
            "vntyper_confidence": confidence,
            "vntyper_variant": variant,
            "vntyper_alt_depth": alt_depth,
            "notes": notes,
        }
    ]
    bioscript.write_tsv(output_file, rows_out)
    bioscript.write_tsv("/output/kestrel/kestrel_result.tsv", rows)
    bioscript.write_tsv("/output/kestrel/kestrel_pre_result.tsv", rows)
    bioscript.write_text("/output/kestrel/output.vcf", bioscript.read_text(kestrel_vcf))
    bioscript.write_text("/output/kestrel/output_insertion.vcf", bioscript.read_text(kestrel_vcf))
    bioscript.write_text("/output/kestrel/output_deletion.vcf", "##fileformat=VCFv4.2\n")
    bioscript.write_text("/output/kestrel/output_indel.vcf", bioscript.read_text(kestrel_vcf))
    bioscript.write_text(
        "/output/kestrel/output.bed",
        "chr1\t155158000\t155163000\tMUC1\n",
    )
    bioscript.write_tsv(
        "/output/coverage/coverage_summary.tsv",
        [
            {
                "file": input_file,
                "region": MUC1_VNTR_REGION,
                "mean": depth_summary["mean"],
                "median": depth_summary["median"],
                "stdev": depth_summary["stdev"],
                "min": depth_summary["min"],
                "max": depth_summary["max"],
                "region_length": depth_summary["region_length"],
                "uncovered_bases": depth_summary["uncovered_bases"],
                "percent_uncovered": depth_summary["percent_uncovered"],
            }
        ],
    )
    bioscript.write_text(
        "/output/coverage/coverage_vntr_coverage.txt",
        "region\tmean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\n"
        + MUC1_VNTR_REGION
        + "\t"
        + str(depth_summary["mean"])
        + "\t"
        + str(depth_summary["median"])
        + "\t"
        + str(depth_summary["stdev"])
        + "\t"
        + str(depth_summary["min"])
        + "\t"
        + str(depth_summary["max"])
        + "\t"
        + str(depth_summary["region_length"])
        + "\t"
        + str(depth_summary["uncovered_bases"])
        + "\t"
        + str(depth_summary["percent_uncovered"])
        + "\n",
    )
    bioscript.write_text("/output/coverage/coverage_vntr_coverage.depth.log", "bioscript samtools.depth_native\n")
    bioscript.copy_file(fastq_1, "/output/fastq_bam_processing/output_R1.fastq.gz")
    bioscript.copy_file(fastq_2, "/output/fastq_bam_processing/output_R2.fastq.gz")
    bioscript.copy_file(sliced_bam, "/output/fastq_bam_processing/output_sliced.bam")
    bioscript.copy_file(sliced_bai, "/output/fastq_bam_processing/output_sliced.bam.bai")
    bioscript.write_text("/output/fastq_bam_processing/output_slice.log", "bioscript samtools.view_region_native\n")
    bioscript.write_text("/output/fastq_bam_processing/output_sort_fastq.log", "bioscript samtools.fastq_native\n")
    bioscript.write_text("/output/fastq_bam_processing/output_index.log", "bioscript samtools.index\n")
    bioscript.write_text("/output/fastq_bam_processing/output_merge.log", "not generated by current BioScript VNtyper path\n")
    write_json_text(
        "/output/fastq_bam_processing/pipeline_info.json",
        {
            "bam_region": bam_region,
            "fastq_1": "fastq_bam_processing/output_R1.fastq.gz",
            "fastq_2": "fastq_bam_processing/output_R2.fastq.gz",
            "sliced_bam": "fastq_bam_processing/output_sliced.bam",
        },
    )
    bioscript.write_text(
        "/output/predefined_regions_hg19.bed",
        "chr1\t155158000\t155163000\tMUC1\n",
    )
    bioscript.write_text(
        "/output/pipeline.log",
        "BioScript VNtyper MUC1 pipeline\n"
        + "input: "
        + input_file
        + "\nregion: "
        + bam_region
        + "\nstatus: "
        + status
        + "\n",
    )
    write_json_text(
        "/output/pipeline_summary.json",
        {
            "sample": sample,
            "status": status,
            "outcome": outcome,
            "confidence": confidence,
            "variant": variant,
            "alt_depth": alt_depth,
            "bam_region": bam_region,
            "vntr_region": MUC1_VNTR_REGION,
        },
    )
    write_summary_html(
        "/output/summary_report.html",
        sample,
        status,
        outcome,
        confidence,
        variant,
        alt_depth,
        notes,
        rows,
    )
    bioscript.write_text(
        "/output/igv_report.html",
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>VNtyper IGV report</title></head>"
        + "<body><h1>VNtyper IGV report</h1><p>IGV.js view is not yet generated by the BioScript port.</p></body></html>\n",
    )
    print(outcome)


if __name__ == "__main__":
    main()
