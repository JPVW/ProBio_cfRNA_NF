process SAMTOOLS_SORT {
    label 'process_low'
    
    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/QC/preseq/${sample_id}" , mode: 'copy'

    input:
        tuple val(sample_id), path(input_bam), path(dedup_bai)

    output:
        tuple val(sample_id), path("${sample_id}_dedup_Aligned.sortedByName.out.bam") , emit: sorted_bam
    script:
    """
    samtools view -h -F 0x904 ${input_bam} | samtools sort -n -o ${sample_id}_dedup_Aligned.sortedByName.out.bam 
    """
}