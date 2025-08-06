process BEDTOOLS_INTERSECT {
    label 'process_single'
    
    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/STAR_2pass/${sample_id}" , mode: 'copy'

    input:
       tuple val(sample_id), path(input_bam)
       path gtf_GRCh37

    output:
        tuple val(sample_id), path("*_Aligned_exon_junction.sortedByCoord.out.bam"), emit: junction_bam

    script:
    """
    bedtools intersect -split -abam $input_bam -b $gtf_GRCh37 > ${sample_id}_Aligned_exon_junction.sortedByCoord.out.bam
    """
}