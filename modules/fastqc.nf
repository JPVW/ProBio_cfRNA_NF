#!/usr/bin/env nextflow

process FASTQC {

    container "staphb/fastqc:latest"
    publishDir "results/fastqc/raw_data", mode: 'copy'

    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads} --outdir results/fastqc/raw_data
    """
}
