process VEP {
    label 'process_medium'
    
    container "https://depot.galaxyproject.org/singularity/ensembl-vep:114.0--pl5321h2a3209d_0"
    publishDir "results/RNA_variant_calling/variant_filtered/VEP/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)
    val   assembly
    val   species
    val   cache_version
    path  cache

    output:
    tuple val(sample_id), path("*.vep.vcf.gz")       , optional:true, emit: vcf
    tuple val(sample_id), path("*.vep.vcf.gz.tbi")   , optional:true, emit: tbi

    script:
    """
    vep \\
        -i $vcf \\
        -o ${sample_id}.vep.vcf.gz \\
        --assembly $assembly \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache /kyukon/scratch/gent/vo/002/gvo00224/TOBI/Projects/AR-burden/AR-Burden/ProBio_cfRNA_NF/resources/vep_cache/ \\

    tabix -p vcf ${sample_id}.vep.vcf.gz
    """

}