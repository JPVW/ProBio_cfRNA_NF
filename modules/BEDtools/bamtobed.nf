process BAMTOBED {
    label 'process_low'

    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/STAR_2pass/${sample_id}" , mode: 'copy'
    
    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.bed")                         , emit: bed
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.bed")           , emit: sorted_bed

    script:
    """
    bamToBed -i ${dedup_bam} > ${sample_id}_Aligned.bed
    sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 ${sample_id}_Aligned.bed > ${sample_id}_Aligned.sortedByCoord.bed
    """
}