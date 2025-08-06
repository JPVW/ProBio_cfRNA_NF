process INFER_EXPERIMENT {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0'
    publishDir "results/QC/RSeQC/infer_experiment/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(bed)
    path(bed12)

    output:
    tuple val(sample_id), path("${sample_id}.infer_experiment.txt"), emit: txt

    script:
    """
    infer_experiment.py -i ${dedup_bam} -r $bed > ${sample_id}.infer_experiment.txt
    """
}