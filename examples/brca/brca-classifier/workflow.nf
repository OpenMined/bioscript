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
        def per_participant_results = brca_classifier(
            participant_work_items
        )

        // Aggregate all results into single file
        def aggregated = aggregate_results(
            per_participant_results.collect()
        )

    emit:
        classification_result = aggregated
}

process brca_classifier {
    container 'ghcr.io/openmined/bioscript:0.1.1'
    publishDir params.results_dir, mode: 'copy', overwrite: true, pattern: 'result_BRCA_*.tsv'
    tag { participant_id }

    input:
        tuple path(assets_dir), val(participant_id), path(genotype_file)

    output:
        path "result_BRCA_${participant_id}.tsv"

    script:
    def genoFileName = genotype_file.getName()
    """
    GENO_FILE=$(printf '%q' "${{genoFileName}}")
    bioscript classify "${{assets_dir}}/classify_brca.py" --file $GENO_FILE --participant_id "${{participant_id}}"
    """
}

process aggregate_results {
    container 'ghcr.io/openmined/bioscript:0.1.1'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path individual_results

    output:
        path "result_BRCA.tsv"

    script:
    """
    # Extract header from first file
    head -n 1 ${individual_results[0]} > result_BRCA.tsv

    # Append all data rows (skip headers)
    for file in ${individual_results}; do
        tail -n +2 "\$file" >> result_BRCA.tsv
    done
    """
}
