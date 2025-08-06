process MERGEVCFS {
    label 'process_low'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    
    publishDir "results/RNA_variant_calling/", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("Combined_samples.vcf.gz")       , emit: vcf
    tuple val(sample_id), path("Combined_samples.vcf.gz.tbi")   , optional: true, emit: tbi

    script:
    def vcf_list = vcfs.collect { f -> "-I ${f}" }.join(' ')

    """
    gatk MergeVcfs \\
        ${vcf_list} \\
        --OUTPUT Combined_samples.vcf.gz 
    """
}