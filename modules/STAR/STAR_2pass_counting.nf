process STAR_2PASS_COUNTING {
    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/STAR_2pass/${sample_id}/dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)
    path(sorted_bam_bai) 
    path(junction_bam)
    path(junction_bam_bai)
    path(gtf)

    output:
    tuple val(sample_id), path("${sample_id}_dedup_genes.csv")                                      , emit: dedup_genes
    tuple val(sample_id), path("${sample_id}_dedup_unique_genes.csv")                               , emit: unique_genes
    tuple val(sample_id), path("${sample_id}_dedup_exon_junction_genes.csv")                        , emit: exon_genes
    tuple val(sample_id), path("${sample_id}_dedup_Aligned.sortedByCoord.out.bam")                  , emit: dedup_bam
    tuple val(sample_id), path("${sample_id}_dedup_exon_junction_Aligned.sortedByCoord.out.bam")    , emit: dedup_exon_bam
    tuple val(sample_id), path("${sample_id}_dedup_unique_mapped.bam")                              , emit: dedup_unique_bam
    tuple val(sample_id), path("${sample_id}_dedup_Aligned.sortedByCoord.out.bam.bai")              , emit: dedup_bai
    tuple val(sample_id), path("${sample_id}_dedup_exon_junction_Aligned.sortedByCoord.out.bam.bai"), emit: dedup_exon_bai
    tuple val(sample_id), path("${sample_id}_dedup_unique_mapped.bam.bai")                          , emit: dedup_unique_bai

    script:
    """
    #!/bin/bash

        umi_tools dedup --method directional -I ${sorted_bam} \\
            --output-stats=${sample_id} \\
            --spliced-is-unique --paired --log=${sample_id}_dedup.log \\
            -S ${sample_id}_dedup_Aligned.sortedByCoord.out.bam

        umi_tools dedup --method directional -I ${junction_bam} \\
            --output-stats=${sample_id} \\
            --spliced-is-unique --paired --log=${sample_id}_split_dedup.log \\
            -S ${sample_id}_dedup_exon_junction_Aligned.sortedByCoord.out.bam

        samtools view -hq 255 ${sample_id}_dedup_Aligned.sortedByCoord.out.bam | samtools sort - -o ${sample_id}_dedup_unique_mapped.bam
        samtools index ${sample_id}_dedup_Aligned.sortedByCoord.out.bam
        samtools index ${sample_id}_dedup_exon_junction_Aligned.sortedByCoord.out.bam
        samtools index ${sample_id}_dedup_unique_mapped.bam

        featureCounts -p -a ${gtf} -s 2 -g gene_id -T 16 -t exon -o ${sample_id}_dedup_gene_featureCounts.txt ${sample_id}_dedup_Aligned.sortedByCoord.out.bam
        cut -f1,6,7- ${sample_id}_dedup_gene_featureCounts.txt | sed 1d > ${sample_id}_dedup_genes.csv

        featureCounts -p -a ${gtf} -s 2 -g gene_id -T 16 -t exon -o ${sample_id}_dedup_gene_unique_featureCounts.txt ${sample_id}_dedup_unique_mapped.bam
        cut -f1,6,7- ${sample_id}_dedup_gene_unique_featureCounts.txt | sed 1d > ${sample_id}_dedup_unique_genes.csv

        featureCounts -p -a ${gtf} -s 2 -g gene_id -T 16 -t exon -o ${sample_id}_dedup_gene_exon_junctions_featureCounts.txt ${sample_id}_dedup_exon_junction_Aligned.sortedByCoord.out.bam
        cut -f1,6,7- ${sample_id}_dedup_gene_exon_junctions_featureCounts.txt | sed 1d > ${sample_id}_dedup_exon_junction_genes.csv
    """
}

    