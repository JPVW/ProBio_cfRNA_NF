#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "clinicalgenomics/trim_galore:0.6.7"
    publishDir "results/trimmed_reads", mode: 'copy'

    input:
    tuple path(read1), path(read2)

    output:
    tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
    path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

    script:
    """
    trim_galore --gzip --illumina  --paired -q 0 -o results/trimmed_reads --fastqc ${read1} ${read2}
    """
}
