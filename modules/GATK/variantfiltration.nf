process VARIANTFILTRATION {
    label 'process_single'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    
    publishDir "results/RNA_variant_calling/variant_filtered/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}.tbi")   , optional: true, emit: tbi



    script:
    """
    gatk VariantFiltration \\
        --variant $vcf \\
        --output ${sample_id}.vcf.gz \\
        --reference $fasta \\
        --window 35 \\
        --cluster 3 \\
        --filter-name "FS" \\
        --filter "FS > 30.0" \\
        --filter-name "QD" \\
        --filter "QD < 2.0" \\
    """
}