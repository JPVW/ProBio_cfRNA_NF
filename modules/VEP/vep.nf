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
    tuple val(sample_id), path("*.vep.vcf.gz")              ,   optional:true, emit: vcf
    tuple val(sample_id), path("*.vep.vcf.gz.tbi")          ,   optional:true, emit: tbi
    tuple val(sample_id), path("${sample_id}.vep.stats.txt"),   optional:true, emit: stats

    script:
    """
    vep \\
        -i $vcf \\
        -o ${sample_id}.vep.vcf.gz \\
        --assembly $assembly \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache /groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/vep_cache/ \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --everything \\
        --stats_file ${sample_id}.vep.stats.txt \\

    tabix -p vcf ${sample_id}.vep.vcf.gz
    """

}