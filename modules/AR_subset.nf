process AR_SUBSET {
    container 'docker://jpvw/probio_cfrna_tools:1.0'

    publishDir "results/AR_subset/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bam_bai)

    output:
    tuple val(sample_id), path("${sample_id}_unique_mapped_AR.bam")                 , emit: unique
    tuple val(sample_id), path("${sample_id}_multi_mapped.bam")                     , emit: multi
    tuple val(sample_id), path("${sample_id}_unique_mapped_AR.bam")                 , emit: multi_AR
    tuple val(sample_id), path("${sample_id}_merged_AR.bam")                        , emit: merged
    tuple val(sample_id), path("${sample_id}_merged_AR_sorted.bam")                 , emit: AR_sorted
    tuple val(sample_id), path("${sample_id}_merged_AR_sorted_Chr.bam")             , emit: AR_sortedChr
    tuple val(sample_id), path("${sample_id}_merged_AR_sorted_Chr.bam.bai")         , emit: sorted_bai
    tuple val(sample_id), path("${sample_id}_ARsubset_R1.fq.gz")                    , emit: fq1
    tuple val(sample_id), path("${sample_id}_ARsubset_R2.fq.gz")                    , emit: fq2
    tuple val(sample_id), path("${sample_id}_ARsubset_singleton.fq.gz")             , emit: singleton
    tuple val(sample_id), path("${sample_id}_dedup_merged_AR_sorted.bam")           , emit: dedup_bam
    tuple val(sample_id), path("${sample_id}_AR_dedup.log")                         , emit: AR_log

    script:
    """
    samtools view -hbq 255 ${sorted_bam} 'chrX:66753830-67011796'  | samtools sort -n - -o ${sample_id}_unique_mapped_AR.bam
    samtools view -h ${sorted_bam} | grep -vw 'NH:i:1' | samtools view -bS  - -o ${sample_id}_multi_mapped.bam
    samtools index ${sample_id}_multi_mapped.bam
    samtools view -hb ${sample_id}_multi_mapped.bam 'chrX:66753830-67011796'  | samtools sort -n - -o ${sample_id}_multi_mapped_AR.bam
    samtools merge ${sample_id}_merged_AR.bam ${sample_id}_unique_mapped_AR.bam ${sample_id}_unique_mapped_AR.bam
    samtools sort -n ${sample_id}_merged_AR.bam -o ${sample_id}_merged_AR_sorted.bam
    samtools sort ${sample_id}_merged_AR.bam -o ${sample_id}_merged_AR_sorted_Chr.bam
    samtools index ${sample_id}_merged_AR_sorted_Chr.bam ${sample_id}_merged_AR_sorted_Chr.bam.bai
    samtools fastq ${sample_id}_merged_AR_sorted.bam -1 ${sample_id}_ARsubset_R1.fq.gz -2 ${sample_id}_ARsubset_R2.fq.gz -s ${sample_id}_ARsubset_singleton.fq.gz
    umi_tools dedup --method directional -I ${sample_id}_merged_AR_sorted_Chr.bam --output-stats=${sample_id} --spliced-is-unique --paired --log=${sample_id}_AR_dedup.log -S ${sample_id}_dedup_merged_AR_sorted.bam
    """
}