process SNPEFF_SNPEFF {
    label 'process_medium'
    
    container "https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2"
    publishDir "results/RNA_variant_calling/variant_filtered/snpEff/${sample_id}", mode: 'copy'


    input:
    tuple val(sample_id), path(vcf)
    val   db
    path(cache)

    output:
    tuple val(sample_id), path("${sample_id}.snpeff.ann.vcf"),   emit: vcf
    tuple val(sample_id), path("${sample_id}.snpeff.csv"),       emit: report
    tuple val(sample_id), path("${sample_id}.snpeff.html"),      optional: true, emit: summary_html
    tuple val(sample_id), path("${sample_id}.snpeff.genes.txt"), emit: genes_txt

    script:
    """
    snpEff -Xmx8g \\
        $db \\
        -csvStats ${sample_id}.snpeff.csv \\
        -dataDir \${PWD}/${cache} \\
        $vcf > ${sample_id}.snpeff.ann.vcf
    """

}