process STAR_AR_COUNTING {
  label 'process_high'  

  container 'docker://jpvw/probio_cfrna_tools:1.0'
  publishDir "results/STAR_AR/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(bam_sorted_aligned)
  path(bai_sorted_aligned) 
  path(AR_gtf)

  output:
  tuple val(sample_id), path("${sample_id}_AR_genes.csv")                                 , optional: true, emit: AR_genes
  tuple val(sample_id), path("${sample_id}_AR_featureCounts.txt")                         , optional: true, emit: AR_gene_featureCounts
  tuple val(sample_id), path("${sample_id}_dedup_AR_genes.csv")                           , optional: true, emit: dedup_AR_genes
  tuple val(sample_id), path("${sample_id}_dedup_AR_Aligned,sortedByCoord.out.bam")       , optional: true, emit: AR_bam_dedup
  tuple val(sample_id), path("${sample_id}_dedup_AR_Aligned,sortedByCoord.out.bam.bai")   , optional: true, emit: AR_bai_dedup
  tuple val(sample_id), path("${sample_id}_dedup_AR_featureCounts.txt")                   , optional: true, emit: dedup_AR_featureCounts
  tuple val(sample_id), path("${sample_id}_AR_dedup.log")                                 , optional: true, emit: AR_log

  script:
  """
  #!/bin/bash

  if [ \$(samtools flagstat "${bam_sorted_aligned}" | head -n 1 | awk '{print \$1}') -gt 0 ]; then
  echo "Sample has reads, performing deduplication and counting"

  umi_tools dedup --method directional \\
    -I ${bam_sorted_aligned} \\
    --output-stats=${sample_id} \\
    --spliced-is-unique --paired \\
    --log=${sample_id}AR_dedup.log \\
    -S "${sample_id}_dedup_AR_Aligned,sortedByCoord.out.bam"

  samtools index ${sample_id}_dedup_AR_Aligned,sortedByCoord.out.bam ${sample_id}_dedup_AR_Aligned,sortedByCoord.out.bam.bai

  featureCounts -p -a $AR_gtf \\
    -s 2 -g gene_id -T 16 -t exon \\
    -o ${sample_id}_dedup_AR_featureCounts.txt ${sample_id}_dedup_AR_Aligned,sortedByCoord.out.bam

  cut -f1,6,7- ${sample_id}_dedup_AR_featureCounts.txt | sed 1d > ${sample_id}_dedup_AR_genes.csv


  featureCounts -p -a $AR_gtf \\
    -s 2 -g gene_id -T 16 -t exon \\
    -o ${sample_id}_AR_featureCounts.txt ${bam_sorted_aligned}

  cut -f1,6,7- ${sample_id}_AR_featureCounts.txt | sed 1d > ${sample_id}_AR_genes.csv

  else
  echo "Sample is empty, skipping deduplication"

  fi
  """
}

    