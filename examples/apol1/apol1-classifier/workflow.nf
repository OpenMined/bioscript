nextflow.enable.dsl=2

workflow USER {
    take:
        context
        participants  // Channel emitting GenotypeRecord maps

    main:
        def assetsDir = file(context.params.assets_dir)
        def workflowScript = file("${assetsDir}/classify_apol1.py")

        // Extract (participant_id, genotype_file) tuples from the records channel
        def participant_tuples = participants.map { record ->
            tuple(
                record.participant_id,
                file(record.genotype_file)
            )
        }

        // Process each participant
        def per_participant_results = apol1_classifier(
            workflowScript,
            participant_tuples
        )

        // Aggregate all results into single file
        def aggregated = aggregate_results(
            per_participant_results.collect()
        )

    emit:
        classification_result = aggregated
}

process apol1_classifier {
    container 'ghcr.io/openmined/bioscript:0.1.1'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'result_APOL1_*.tsv'
    tag { participant_id }

    input:
        path script
        tuple val(participant_id), path(genotype_file)

    output:
        path "result_APOL1_${participant_id}.tsv"

    script:
    """
    bioscript classify "${script}" --file "${genotype_file}" --participant_id "${participant_id}"
    """
}

process aggregate_results {
    container 'ghcr.io/openmined/bioscript:0.1.1'
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
