#!/usr/bin/env nextflow

process FASTQC {

    container "staphb/fastqc:latest"
    publishDir "results/fastqc/raw_data/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc $read1 $read2
    """
}
