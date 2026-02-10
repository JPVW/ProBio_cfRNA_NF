process BASERECALIBRATOR {
    label 'process_low'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    
    publishDir "results/RNA_variant_calling/baserecalibrator/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(split_bam)
    path(fasta)
    path(fai)
    path(dict)
    path(known_vcf)
    path(vcf_tbi)

    output:
    tuple val(sample_id), path (split_bam), path("${sample_id}_baserecalibrator.table"), emit: table
    
    script:
    """
    gatk BaseRecalibrator \
        -R ${fasta} \
        -I ${split_bam} \
        --known-sites ${known_vcf} \
        --use-original-qualities \
        -O ${sample_id}_baserecalibrator.table
    """
}