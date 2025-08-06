process HSMETRICS {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0'
    publishDir "results/QC/HSmetrics/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(fasta)
    path(fai)
    path(target)

    output:
    tuple val(sample_id), path("${sample_id}_output_HSmetrics.txt")                     , emit: metrics
    
    script:
    """    
    picard -Xmx12g CollectHsMetrics \\
      -I ${dedup_bam} \\
      -O ${sample_id}_output_HSmetrics.txt \\
      -R $fasta \\
      -BAIT_INTERVALS $target \\
      -TARGET_INTERVALS $target
    """
}
