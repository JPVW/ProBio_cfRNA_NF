process SNPEFF_DOWNLOAD {
    label 'process_medium'
    
    container "https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2"
    publishDir 'resources/', mode: 'copy'

    input:
    val(snpeff_db)

    output:
    path('snpeff_cache'), emit: cache

    script:
    """
    snpEff download -v -dataDir \${PWD}/snpeff_cache  $snpeff_db
    """
}