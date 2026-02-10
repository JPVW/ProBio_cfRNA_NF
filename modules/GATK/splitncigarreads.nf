process SPLITNCIGARREADS {
    label 'process_single'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    
    publishDir "results/RNA_variant_calling/splitncigarreads/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(sample_id), path("${sample_id}_SplitNCigar.bam")   , emit: split_bam
    
    script:
    """
    gatk SplitNCigarReads \
        -R ${fasta} \
        -I ${dedup_bam} \
        -O ${sample_id}_SplitNCigar.bam
    """
}