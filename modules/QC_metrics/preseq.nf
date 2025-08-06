process PRESEQ {
    label 'process_low'

    container 'https://depot.galaxyproject.org/singularity/preseq:3.2.0--hdcf5f25_6'
    publishDir "results/QC/preseq/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.lc_extrap.txt")        , emit: lc_extrap
    tuple val(sample_id), path("${sample_id}.log")                  , optional: true, emit: log
    tuple val(sample_id), path("${sample_id}_ccurve.txt")           , emit: ccurve
    tuple val(sample_id), path("${sample_id}_yield.txt")            , optional: true, emit: yield

    script:
    """
    preseq c_curve -P -o ${sample_id}_ccurve.txt ${bed}
    preseq lc_extrap -P -D -v -output ${sample_id}.lc_extrap.txt ${bed}
    """
}