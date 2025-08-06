process HAPLOTYPECALLER {
    label 'process_low'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    
    publishDir "results/RNA_variant_calling/haplotypecaller/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(BQSR_bam) 
    path(fasta)
    path(fai)
    path(dict)
    path(known_vcf)
    path(vcf_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")       , emit: vcf
    tuple val(sample_id), path("${sample_id}.vcf.gz.tbi")   , optional:true, emit: tbi
    tuple val(sample_id), path("${sample_id}.realigned.bam"), optional:true, emit: bam
    
    script:
    """
    gatk HaplotypeCaller \
        --input ${BQSR_bam} \
        --output ${sample_id}.vcf.gz \
        --reference $fasta \
        -dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        --dbsnp $known_vcf
    """
}