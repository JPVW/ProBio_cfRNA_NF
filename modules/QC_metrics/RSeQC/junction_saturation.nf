process JUNCTION_SATURATION {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
    publishDir "results/QC/RSeQC/junction_saturation/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(bed)
    path(bed12)

    output:
    tuple val(sample_id), path("*distance.txt"), optional:true, emit: distance

    script:
    """
    junction_saturation.py -i ${dedup_bam} -r $bed12 -o ${sample_id}
    """
}