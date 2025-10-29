nextflow.enable.dsl=2

workflow USER {
    take:
        context
        genotype_file
        data_dir

    main:
        def assetsDir = file(context.params.assets_dir)
        def workflowScript = file("${assetsDir}/classify_herc2.py")
        def script_ch = Channel.value(workflowScript)
        def classification_result_ch = herc2_classifier(
            script_ch,
            genotype_file,
            data_dir
        )

    emit:
        classification_result = classification_result_ch
}

process herc2_classifier {
    container 'ghcr.io/openmined/bioscript:0.1.1'
    publishDir params.results_dir, mode: 'copy', overwrite: true

    input:
        path script
        path genotype_file
        path data_dir

    output:
        path 'result_HERC2_{participant_id}.tsv', emit: classification_result

    script:
    """
    python3 ${script} \n        --input "${genotype_file}"
        --data-dir "${data_dir}"
        --output "result_HERC2_{participant_id}.tsv"
    """
}
