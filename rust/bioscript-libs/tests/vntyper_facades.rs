use std::io::{Read, Write};

use bioscript_libs::{
    bcftools,
    kestrel::native::{
        NativeKestrelRunOptions, NativeReferenceRegion, call_fastq_paths_to_vcf_references,
    },
    samtools,
    vcf::parse_kestrel_vcf,
};

#[test]
fn native_vntyper_facades_can_extract_fastq_call_and_sort_vcf() {
    let temp = tempfile::tempdir().unwrap();
    let sam = temp.path().join("reads.sam");
    let bam = temp.path().join("reads.bam");
    let fastq_1 = temp.path().join("reads_R1.fastq.gz");
    let fastq_2 = temp.path().join("reads_R2.fastq.gz");
    let calls = temp.path().join("calls.vcf");
    let sorted = temp.path().join("calls.sorted.vcf.gz");
    let sorted_index = temp.path().join("calls.sorted.vcf.gz.csi");

    write_variant_pair_sam(&sam);
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();
    samtools_rs::native::index(&bam, Option::<&std::path::Path>::None, Some(1)).unwrap();

    let fastq = samtools::fastq_native(&bam, None, "chr1:1-16", &fastq_1, &fastq_2).unwrap();
    assert_eq!(fastq.read1_records, 5);
    assert_eq!(fastq.read2_records, 5);

    let mut options = NativeKestrelRunOptions::new("sample1");
    options.minimum_difference = 1;
    options.max_haplotypes = 4;
    options.max_saved_states = 4;

    let vcf = call_fastq_paths_to_vcf_references(
        &[NativeReferenceRegion::new(
            "chr1",
            "AAAACCCCGGGGTTTT",
            "2a9fd43653a81f9ec44e34c7ec038636",
        )],
        [fastq_1.as_path(), fastq_2.as_path()],
        4,
        &options,
    )
    .unwrap();
    assert!(vcf.contains("chr1\t5\t.\tC\tT"), "{vcf}");
    // kestrel-rs is now bug-compatible with Java Kestrel (vendor PR #3),
    // which emits the full motif-reference-equivalent record set rather
    // than a single collapsed row. The canonical C>T call is asserted
    // above; the parsed set is the Java-parity output, not 1.
    assert_eq!(parse_kestrel_vcf(&vcf).unwrap().len(), 7);

    std::fs::write(&calls, vcf).unwrap();
    bcftools::sort_native(&calls, &sorted, "z", true).unwrap();
    assert!(std::fs::metadata(&sorted_index).unwrap().len() > 0);

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(sorted).unwrap());
    let mut sorted_vcf = String::new();
    decoder.read_to_string(&mut sorted_vcf).unwrap();
    assert!(sorted_vcf.contains("chr1\t5\t.\tC\tT"), "{sorted_vcf}");
}

fn write_variant_pair_sam(path: &std::path::Path) {
    let mut file = std::fs::File::create(path).unwrap();
    writeln!(file, "@HD\tVN:1.6\tSO:coordinate").unwrap();
    writeln!(file, "@SQ\tSN:chr1\tLN:16").unwrap();
    for index in 0..5 {
        writeln!(
            file,
            "r{index}\t65\tchr1\t1\t60\t16M\t=\t1\t0\tAAAATCCCGGGGTTTT\tIIIIIIIIIIIIIIII"
        )
        .unwrap();
        writeln!(
            file,
            "r{index}\t129\tchr1\t1\t60\t16M\t=\t1\t0\tAAAATCCCGGGGTTTT\tIIIIIIIIIIIIIIII"
        )
        .unwrap();
    }
}
