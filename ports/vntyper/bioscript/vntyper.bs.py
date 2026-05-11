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

    report = {
        "participant_id": sample,
        "samtools_fastq_command": fastq_command,
        "kestrel_command": kestrel_command,
    }
    bioscript.write_tsv(output_file, [report])


if __name__ == "__main__":
    main()
