process TABIX {
    label 'process_single'
    
    container 'https://depot.galaxyproject.org/singularity/htslib:1.20--h5efdd21_2'
    
    publishDir "results/RNA_variant_calling/haplotypecaller/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(tab)

    output:
    tuple val(sample_id), path("*.tbi"), optional:true, emit: tbi
    tuple val(sample_id), path("*.csi"), optional:true, emit: csi

    script:
    """
    tabix $tab
    """
}