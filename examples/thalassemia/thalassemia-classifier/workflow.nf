// BioVault workflow export v0.1.1

nextflow.enable.dsl=2

workflow USER {
    take:
        context
        participants  // Channel emitting GenotypeRecord maps

    main:
        def assetsDir = context.assets_dir
        if (!assetsDir) {
            throw new IllegalStateException("Missing assets directory in context")
        }
        def assetsDirPath = file(assetsDir)

        // Pair the assets directory with each (participant_id, genotype_file) tuple
        def participant_work_items = participants.map { record ->
            tuple(
                assetsDirPath,
                record.participant_id,
                file(record.genotype_file)
            )
        }

        // Process each participant
        def per_participant_results = thalassemia_classifier(
            participant_work_items
        )

        // Aggregate all results into single file
        def aggregated = aggregate_results(
            per_participant_results.collect()
        )

    emit:
        classification_result = aggregated
}

process thalassemia_classifier {
    container 'ghcr.io/openmined/bioscript:0.1.5'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'result_THALASSEMIA_*.tsv'
    tag { participant_id }
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple path(assets_dir), val(participant_id), path(genotype_file)

    output:
        path "result_THALASSEMIA_${participant_id}.tsv"

    script:
    def genoFileName = genotype_file.getName()
    """
    GENO_FILE=\$(printf '%q' "${genoFileName}")
    bioscript classify "${assets_dir}/classify_thalassemia.py" --file \$GENO_FILE --participant_id "${participant_id}"
    """
}

process aggregate_results {
    container 'ghcr.io/openmined/bioscript:0.1.5'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path individual_results

    output:
        path "result_THALASSEMIA.tsv"

    script:
    def manifestContent = individual_results.collect { it.toString() }.join('\n') + '\n'
    """
    cat <<'EOF' > results.list\n${manifestContent}EOF
    bioscript combine --list results.list --output result_THALASSEMIA.tsv
    """
}
