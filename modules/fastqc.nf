process FASTQC {
    label 'process_medium'
    
    container "staphb/fastqc:latest"
    publishDir "results/fastqc/raw_data/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*.html"), optional:true, emit: html
    tuple val(sample_id), path("*.zip") , optional: true, emit: zip

    script:
    """
    mkdir -p results/fastqc/raw_data/
    fastqc $read1 $read2
    """
}
