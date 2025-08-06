process VEP_DOWNLOAD {
    label 'process_medium'
    
    container "https://depot.galaxyproject.org/singularity/ensembl-vep:113.4--pl5321h2a3209d_0"
    publishDir 'resources/', mode: 'copy'

    input:
    val(assembly)
    val(species)
    val(version)

    output:
    path('vep_cache'), emit: cache


    script:
    
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    vep_install \
        --AUTO c \
        --CACHEDIR \${PWD}/vep_cache \
        --SPECIES $species \
        --ASSEMBLY $assembly \
        --CACHE_VERSION $version \
        --NO_UPDATE \
        --NO_HTSLIB \
        --NO_PLUGINS \
        --NO_FA
    """
}