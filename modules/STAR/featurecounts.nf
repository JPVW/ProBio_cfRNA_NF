process FEATURECOUNTS {
    label 'process_high'
    
    container 'docker://jpvw/probio_cfrna_tools:1.0'
    publishDir "results/STAR_2pass/${sample_id}/dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam) 
    tuple val(sample_id), path(dedup_exon_bam)
    tuple val(sample_id), path(dedup_unique_bam)

    path(gtf)

    output:
    tuple val(sample_id), path("${sample_id}_dedup_genes.csv")                                      , emit: dedup_genes
    tuple val(sample_id), path("${sample_id}_dedup_unique_genes.csv")                               , emit: unique_genes
    tuple val(sample_id), path("${sample_id}_dedup_exon_junction_genes.csv")                        , emit: exon_genes

    script:
    """
    #!/bin/bash
        featureCounts -p -a ${gtf} -s 2 -g gene_id -T 16 -t exon -o ${sample_id}_dedup_gene_featureCounts.txt ${dedup_bam}
        cut -f1,6,7- ${sample_id}_dedup_gene_featureCounts.txt | sed 1d > ${sample_id}_dedup_genes.csv

        featureCounts -p -a ${gtf} -s 2 -g gene_id -T 16 -t exon -o ${sample_id}_dedup_gene_unique_featureCounts.txt ${dedup_unique_bam}
        cut -f1,6,7- ${sample_id}_dedup_gene_unique_featureCounts.txt | sed 1d > ${sample_id}_dedup_unique_genes.csv

        featureCounts -p -a ${gtf} -s 2 -g gene_id -T 16 -t exon -o ${sample_id}_dedup_gene_exon_junctions_featureCounts.txt ${dedup_exon_bam}
        cut -f1,6,7- ${sample_id}_dedup_gene_exon_junctions_featureCounts.txt | sed 1d > ${sample_id}_dedup_exon_junction_genes.csv
    """
}

    