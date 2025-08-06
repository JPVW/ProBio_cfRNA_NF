process COLLECTALIGNMENTMETRICS {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0'
    publishDir "results/QC/CollectAlignmentSummaryMetrics/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_AlignmentSummaryMetrics.txt")                     , emit: CollectAlignmentSummary_metrics
    
    script:
    """
    picard CollectAlignmentSummaryMetrics \\
        I=${dedup_bam} \\
        O=${sample_id}_AlignmentSummaryMetrics.txt \\
        R=$fasta \\
        VALIDATION_STRINGENCY=LENIENT \\
        METRIC_ACCUMULATION_LEVEL=ALL_READS
    """
}