process EXTRACT_UMI {
    label "process_single"
    label "process_long"

    container "https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0"
    publishDir "results/UMI_extracted/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*.fastq.gz")    , emit: UMI_extracted
    path("*.log")         , emit: log

    script:
    """
    umi_tools extract --extract-method=regex -I ${read2} --bc-pattern "^(?P<umi_1>.{8})(?P<discard_1>.{6}).*" --read2-in=${read1} -S ${sample_id}_val_2_UMIextract.fastq.gz --read2-out=${sample_id}_val_1_UMIextract.fastq.gz -L ${sample_id}_UMIextract.log
    """
}