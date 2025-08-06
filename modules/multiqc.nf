process MULTIQC {
    label 'process_single'

    container "https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0"
    publishDir "results/multiqc", mode: 'copy'

    input:
    path '*/'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}
