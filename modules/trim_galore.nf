process TRIM_GALORE {
    label 'process_high'
    
    container "https://depot.galaxyproject.org/singularity/trim-galore%3A0.6.10--hdfd78af_1"
    publishDir "results/trimmed_reads/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*_trimming_report.txt"), emit: trimming_reports
    tuple val(sample_id), path("*_val_1_fastqc.{zip,html}"), emit: fastqc_reports_1
    tuple val(sample_id), path("*_val_2_fastqc.{zip,html}"), emit: fastqc_reports_2

    script:
    """
    trim_galore --gzip --illumina  --paired -q 0 --fastqc $reads
    """
}
