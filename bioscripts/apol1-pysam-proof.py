from bioscript import pysam


APOL1_SITES = [
    {
        "key": "G1_SITE_1",
        "chrom": "22",
        "start": 36265859,
        "stop": 36265860,
        "ref": "A",
        "alt": "G",
    },
    {
        "key": "G1_SITE_2",
        "chrom": "22",
        "start": 36265987,
        "stop": 36265988,
        "ref": "T",
        "alt": "G",
    },
    {
        "key": "G2_SITE",
        "chrom": "22",
        "start": 36265999,
        "stop": 36266005,
        "ref": "TTATAA",
        "alt": "<DEL:6>",
    },
]


def count_region_reads(bam, site):
    total = 0
    for read in bam.fetch(site["chrom"], site["start"], site["stop"]):
        if not read.is_unmapped:
            total = total + 1
    return total


def main():
    bam = pysam.AlignmentFile(
        input_file,
        "rc",
        reference_filename=reference_file,
        index_filename=input_index,
    )
    rows = []
    for site in APOL1_SITES:
        rows.append(
            {
                "participant_id": participant_id,
                "variant_key": site["key"],
                "chrom": site["chrom"],
                "start": str(site["start"]),
                "stop": str(site["stop"]),
                "depth": str(count_region_reads(bam, site)),
                "proof_status": "region_fetch_only",
            }
        )
    bioscript.write_tsv(output_file, rows)


if __name__ == "__main__":
    main()
