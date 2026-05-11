from bioscript import bcftools
from bioscript import kestrel
from bioscript import samtools
from bioscript import vcf


def main():
    sample = participant_id
    work_dir = "vntyper"

    fastq_1 = work_dir + "/" + sample + "_R1.fastq.gz"
    fastq_2 = work_dir + "/" + sample + "_R2.fastq.gz"

    fastq_command = samtools.fastq(
        input_file,
        fastq_1,
        fastq_2,
    )

    kestrel_command = kestrel.build_command(
        "ports/vntyper/kestrel/kestrel.jar",
        "ports/vntyper/vntyper/reference/All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa",
        work_dir + "/kestrel/output.vcf",
        work_dir + "/kestrel/output.sam",
        work_dir + "/kestrel/tmp",
        sample,
        fastq_1,
        fastq_2,
    )
    sorted_vcf = work_dir + "/kestrel/output.sorted.vcf.gz"
    bcftools_sort_command = bcftools.sort(work_dir + "/kestrel/output.vcf", sorted_vcf)
    bcftools_index_command = bcftools.index(sorted_vcf)

    report = {
        "participant_id": sample,
        "samtools_fastq_command": fastq_command,
        "kestrel_command": kestrel_command,
        "bcftools_sort_command": bcftools_sort_command,
        "bcftools_index_command": bcftools_index_command,
    }
    bioscript.write_tsv(output_file, [report])


if __name__ == "__main__":
    main()
