process SAMTOOLS_INDEX {

    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/STAR_2pass/${sample_id}" , mode: 'copy'

    input:
        tuple val(sample_id), path(input_bam)
        tuple val(sample_id), path(input_exon_bam)

    output:
        path("${input_bam}.bai") , emit: bai
        path("${input_exon_bam}.bai") , emit: exon_bai

    script:
    """
    samtools index '$input_bam'
    samtools index '$input_exon_bam'
    """
}