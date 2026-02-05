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
        def per_participant_results = herc2_classifier(
            participant_work_items
        )

        // Aggregate all results into single file
        def aggregated = aggregate_results(
            per_participant_results.collect()
        )

        // Aggregate population statistics
        def population_stats_ch = aggregate_population_stats(
            Channel.value(assetsDirPath),
            aggregated
        )

    emit:
        classification_result = aggregated
        population_stats = population_stats_ch
}

process herc2_classifier {
    container 'ghcr.io/openmined/bioscript:0.1.6'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'result_HERC2_*.tsv'
    tag { participant_id }
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }

    input:
        tuple path(assets_dir), val(participant_id), path(genotype_file)

    output:
        path "result_HERC2_${participant_id}.tsv"

    script:
    def genoFileName = genotype_file.getName()
    """
    bioscript classify "${assets_dir}/classify_herc2.py" --file "${genoFileName}" --participant_id "${participant_id}"
    """
}

process aggregate_results {
    container 'ghcr.io/openmined/bioscript:0.1.6'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path individual_results

    output:
        path "result_HERC2.tsv"

    script:
    def manifestContent = individual_results.collect { it.toString() }.join('\n') + '\n'
    """
    cat <<'EOF' > results.list\n${manifestContent}EOF
    bioscript combine --list results.list --output result_HERC2.tsv
    """
}

process aggregate_population_stats {
    container 'ghcr.io/openmined/bioscript:0.1.6'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path assets_dir
        path aggregated_results

    output:
        path "result_HERC2_stats.tsv"

    script:
    """
    python3 "${assets_dir}/aggregate_population_stats.py"       --input "${aggregated_results}"       --output result_HERC2_stats.tsv
    """
}
