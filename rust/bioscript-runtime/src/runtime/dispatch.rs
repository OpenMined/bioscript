use bioscript_core::RuntimeError;
use monty::MontyObject;

use super::BioscriptRuntime;

impl BioscriptRuntime {
    pub(super) fn dispatch_method_call(
        &self,
        method_name: &str,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let class_name = match args.first() {
            Some(MontyObject::Dataclass { name, .. }) => name.as_str(),
            _ => "<unknown>",
        };

        match (class_name, method_name) {
            ("Bioscript", "load_genotypes") => self.method_load_genotypes(args, kwargs),
            ("Bioscript", "variant") => self.method_variant(args, kwargs),
            ("Bioscript", "query_plan") => self.method_query_plan(args, kwargs),
            ("Bioscript", "write_tsv") => self.method_write_tsv(args, kwargs),
            ("Bioscript", "read_text") => self.method_read_text(args, kwargs),
            ("Bioscript", "write_text") => self.method_write_text(args, kwargs),
            ("Bioscript", "exists") => self.method_exists(args, kwargs),
            ("PysamModule", "AlignmentFile") => self.method_pysam_alignment_file(args, kwargs),
            ("PysamAlignmentFile", "fetch") => self.method_pysam_alignment_file_fetch(args, kwargs),
            ("PyfaidxModule", "Fasta") => self.method_pyfaidx_fasta(args, kwargs),
            ("BcftoolsModule", "sort") => self.method_bcftools_sort_native(args, kwargs),
            ("BcftoolsModule", "plan_sort") => self.method_bcftools_sort(args, kwargs),
            ("BcftoolsModule", "index") => self.method_bcftools_index_native(args, kwargs),
            ("BcftoolsModule", "plan_index") => self.method_bcftools_index(args, kwargs),
            ("BcftoolsModule", "view") => self.method_bcftools_view_native(args, kwargs),
            ("BcftoolsModule", "plan_view") => self.method_bcftools_view(args, kwargs),
            ("BcftoolsModule", "view_filter") => self.method_bcftools_view_filter(args, kwargs),
            ("BcftoolsModule", "plan_view_filter") => {
                self.method_bcftools_view_filter(args, kwargs)
            }
            ("BcftoolsModule", "norm") => self.method_bcftools_norm(args, kwargs),
            ("BcftoolsModule", "plan_norm") => self.method_bcftools_norm(args, kwargs),
            ("BcftoolsModule", "view_header_native") => {
                self.method_bcftools_view_header_native(args, kwargs)
            }
            ("BcftoolsModule", "view_native") => self.method_bcftools_view_native(args, kwargs),
            ("BcftoolsModule", "sort_native") => self.method_bcftools_sort_native(args, kwargs),
            ("BcftoolsModule", "index_native") => self.method_bcftools_index_native(args, kwargs),
            ("VcfModule", "VariantFile") => self.method_vcf_variant_file(args, kwargs),
            ("VcfModule", "read_kestrel") => self.method_vcf_read_kestrel(args, kwargs),
            ("VcfModule", "read_vntyper_kestrel") => {
                self.method_vcf_read_vntyper_kestrel(args, kwargs)
            }
            ("KestrelModule", "build_command") => self.method_kestrel_build_command(args, kwargs),
            ("KestrelModule", "plan_command") => self.method_kestrel_build_command(args, kwargs),
            ("KestrelModule", "run_native") => self.method_kestrel_run_native(args, kwargs),
            ("SamtoolsModule", "view") => self.method_samtools_view_region_native(args, kwargs),
            ("SamtoolsModule", "plan_view") => self.method_samtools_view(args, kwargs),
            ("SamtoolsModule", "view_region") => {
                self.method_samtools_view_region_default_native(args, kwargs)
            }
            ("SamtoolsModule", "plan_view_region") => {
                self.method_samtools_view_region(args, kwargs)
            }
            ("SamtoolsModule", "fastq") => self.method_samtools_fastq_all_native(args, kwargs),
            ("SamtoolsModule", "plan_fastq") => self.method_samtools_fastq(args, kwargs),
            ("SamtoolsModule", "sort") => self.method_samtools_sort_native(args, kwargs),
            ("SamtoolsModule", "plan_sort") => self.method_samtools_sort(args, kwargs),
            ("SamtoolsModule", "depth") => self.method_samtools_depth_native(args, kwargs),
            ("SamtoolsModule", "plan_depth") => self.method_samtools_depth(args, kwargs),
            ("SamtoolsModule", "index") => self.method_samtools_index_native(args, kwargs),
            ("SamtoolsModule", "plan_index") => self.method_samtools_index(args, kwargs),
            ("SamtoolsModule", "faidx") => self.method_samtools_faidx(args, kwargs),
            ("SamtoolsModule", "plan_faidx") => self.method_samtools_faidx(args, kwargs),
            ("SamtoolsModule", "view_region_native") => {
                self.method_samtools_view_region_native(args, kwargs)
            }
            ("SamtoolsModule", "fastq_native") => self.method_samtools_fastq_native(args, kwargs),
            ("SamtoolsModule", "depth_native") => self.method_samtools_depth_native(args, kwargs),
            ("GenotypeFile", "get") => self.method_genotype_get(args, kwargs),
            ("GenotypeFile", "lookup_variant") => self.method_genotype_lookup_variant(args, kwargs),
            ("GenotypeFile", "lookup_variant_details") => {
                self.method_genotype_lookup_variant_details(args, kwargs)
            }
            ("GenotypeFile", "lookup_variants") => {
                self.method_genotype_lookup_variants(args, kwargs)
            }
            ("GenotypeFile", "lookup_variants_details") => {
                self.method_genotype_lookup_variants_details(args, kwargs)
            }
            _ => Err(RuntimeError::Unsupported(format!(
                "'{class_name}' object has no attribute '{method_name}'"
            ))),
        }
    }
}
