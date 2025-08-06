process RNASEQMETRICS {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0'
    publishDir "results/QC/RnaSeqMetrics/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(refFlat)

    output:
    tuple val(sample_id), path("${sample_id}.RNA_metrics")                     , emit: RNA_metrics
    
    script:
    """
    picard CollectRnaSeqMetrics \\
        I=${dedup_bam} \\
        O=${sample_id}.RNA_metrics \\
        REF_FLAT= $refFlat \\
        STRAND=SECOND_READ_TRANSCRIPTION_STRAND
    """
}