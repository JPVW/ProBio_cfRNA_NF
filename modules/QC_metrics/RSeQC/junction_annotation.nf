process JUNCTION_ANNOTATION {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
    publishDir "results/QC/RSeQC/junction_annotation/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(bed)
    path(bed12)

    output:
    tuple val(sample_id), path("*.xls")         , emit: xls
    tuple val(sample_id), path("*.r")           , emit: rscript
    tuple val(sample_id), path("*.log")         , emit: log
    tuple val(sample_id), path("*.junction.bed"), emit: bed
    tuple val(sample_id), path("*.Interact.bed"), emit: interact_bed
    tuple val(sample_id), path("*junction.pdf") , optional: true, emit: pdf
    tuple val(sample_id), path("*events.pdf")   , optional: true, emit: events_pdf

    script:
    """
    junction_annotation.py -i ${dedup_bam} -r $bed -o ${sample_id} \\
    2> >(grep -v 'E::idx_find_and_load' | tee ${sample_id}.junction_annotation.log >&2)
    """
}