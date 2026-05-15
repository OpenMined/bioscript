# MUC1 VNTR (VNtyper) advanced assay analysis.
#
# Unlike a SNP-lookup assay, this does not use bioscript.query_plan. It treats
# the dragged-in genome (input_file) as an aligned BAM/CRAM and runs the native
# samtools -> kestrel -> bcftools pipeline that was built on the
# madhava/libs branch, then VNtyper post-processing via vcf.read_vntyper_kestrel.
#
# Runtime globals provided by the assay runner (same input/output contract as
# any assay): input_file, output_file, participant_id, asset_paths.

from bioscript import bcftools
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


def main():
    sample = participant_id
    reference_fasta = asset_paths["muc1_reference"]

    bai = "/work/input.bam.bai"
    sliced_bam = "/work/sliced.bam"
    fastq_1 = "/work/reads_R1.fastq.gz"
    fastq_2 = "/work/reads_R2.fastq.gz"
    kestrel_vcf = "/work/kestrel.vcf"
    sorted_vcf = "/work/kestrel.sorted.vcf.gz"

    samtools.index(input_file, bai)
    samtools.view_region_native(input_file, MUC1_BAM_REGION, sliced_bam, bai)
    samtools.fastq_native(input_file, MUC1_BAM_REGION, fastq_1, fastq_2, bai)

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
    bcftools.sort(kestrel_vcf, sorted_vcf)
    bcftools.index(sorted_vcf)

    rows = vcf.read_vntyper_kestrel(kestrel_vcf)
    called = best_passing_row(rows)

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
    print(outcome)


if __name__ == "__main__":
    main()
