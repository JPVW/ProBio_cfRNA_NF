process SAMTOOLS_IDXSTAT {
    label 'process_single'
    
    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/STAR_2pass/${sample_id}" , mode: 'copy'

    input:
       tuple val(sample_id), path(input_bam) 
       path(bai)

    output:
        tuple val(sample_id), path("*.idxstats"), emit: idxstats

    script:
    """
    samtools idxstats '$input_bam' > ${sample_id}_Aligned.sortedByCoord.out.bam.idxstats
    """
}