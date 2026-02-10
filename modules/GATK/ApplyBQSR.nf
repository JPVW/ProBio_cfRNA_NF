process APPLYBQSR {
    label 'process_low'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    
    publishDir "results/RNA_variant_calling/applyBQSR/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(split_bam), path(bqsr_table)
    path(fasta)
    path(fai)
    path(dict)
    path(known_vcf)
    path(vcf_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.duplicates_marked.recalibrated.bam"), emit: BQSR_bam
    
    script:
    """
    gatk ApplyBQSR \
        --input ${split_bam} \
        --output ${sample_id}.aligned.duplicates_marked.recalibrated.bam \
        --reference ${fasta} \
        --use-original-qualities \
        --bqsr-recal-file ${bqsr_table} 
    """
}