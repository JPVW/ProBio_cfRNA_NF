#!bin/env/usr   nextflow 

process STAR_GENOME {

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.6a--0"
    publishDir 'resources/', mode: 'copy'

    input:
    path fasta_GRCh37
    path gtf_GRCh37
    path sjs

    output:
    path("genomedir") , emit: genomedir

    script:
    """
        STAR \
            --runThreadN 12 \
            --runMode genomeGenerate \
            --genomeDir genomedir/ \
            --genomeFastaFiles $fasta_GRCh37 \
            --sjdbGTFfile $gtf_GRCh37 \
            --limitSjdbInsertNsj 2037800 \
            --sjdbFileChrStartEnd  $sjs \
            --sjdbOverhang 99
    """
}