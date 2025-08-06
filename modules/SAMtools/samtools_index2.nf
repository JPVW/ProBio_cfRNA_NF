process SAMTOOLS_INDEX2 {
    label 'process_low'
    
    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/RNA_variant_calling/${sample_id}" , mode: 'copy'

    input:
        tuple val(sample_id), path(input_bam)

    output:
        tuple val(sample_id), path(input_bam), path("${input_bam}.bai") , emit: bai

    script:
    """
    samtools index ${input_bam}
    """
}