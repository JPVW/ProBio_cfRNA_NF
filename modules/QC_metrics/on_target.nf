process ON_TARGET {
    label 'process_medium'
    
    container 'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0'
    publishDir "results/QC/on_target/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(dedup_bai)
    path(fasta)
    path(fai)
    path(dict)
    path(target)

    output:
    tuple val(sample_id), path("${sample_id}_ontarget_reads_sorted_rmdups.txt")                     , emit: on_target_metrics
    
    script:
    """
    gatk CountReads \
        -R $fasta \
        -I ${dedup_bam} \
        -L $target \
        --read-filter MappedReadFilter \
        --read-filter NotSecondaryAlignmentReadFilter > ${sample_id}_ontarget_reads_sorted_rmdups.txt
    """
}