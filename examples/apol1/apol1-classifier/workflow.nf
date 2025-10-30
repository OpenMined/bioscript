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
        def per_participant_results = apol1_classifier(
            participant_work_items
        )

        // Aggregate all results into single file
        def aggregated = aggregate_results(
            per_participant_results.collect()
        )

    emit:
        classification_result = aggregated
}

process apol1_classifier {
    container 'ghcr.io/openmined/bioscript:0.1.2'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'result_APOL1_*.tsv'
    tag { participant_id }

    input:
        tuple path(assets_dir), val(participant_id), path(genotype_file)

    output:
        path "result_APOL1_${participant_id}.tsv"

    script:
    """
    bioscript classify "${assets_dir}/classify_apol1.py" --file "${genotype_file}" --participant_id "${participant_id}"
    """
}

process aggregate_results {
    container 'ghcr.io/openmined/bioscript:0.1.2'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path individual_results

    output:
        path "result_APOL1.tsv"

    script:
    """
    # Extract header from first file
    head -n 1 ${individual_results[0]} > result_APOL1.tsv

    # Append all data rows (skip headers)
    for file in ${individual_results}; do
        tail -n +2 "\$file" >> result_APOL1.tsv
    done
    """
}
