process INNER_DISTANCE {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
    publishDir "results/QC/RSeQC/Inner_distance/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(bed)
    path(bed12)

    output:
    tuple val(sample_id), path("*distance.txt"), optional:true, emit: distance
    tuple val(sample_id), path("*freq.txt")    , optional:true, emit: freq
    tuple val(sample_id), path("*mean.txt")    , optional:true, emit: mean
    tuple val(sample_id), path("*.pdf")        , optional:true, emit: pdf
    tuple val(sample_id), path("*.r")          , optional:true, emit: rscript

    script:
    """
    inner_distance.py -i ${dedup_bam} -o ${sample_id} -r $bed -u 5000
    """
}