process STAR_AR {
    label 'process_high'

    container 'docker://jpvw/probio_cfrna_tools:1.0'
    publishDir "results/STAR_AR/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1)
    tuple val(sample_id), path(read2)
    path AR_index
    path AR_gtf

    output:
    tuple val(sample_id), path("${sample_id}_AR_Log.final.out")                            , emit: log_final
    tuple val(sample_id), path("${sample_id}_AR_Log.out")                                  , emit: log_out
    tuple val(sample_id), path("${sample_id}_AR_Log.progress.out")                         , emit: log_progress
    tuple val(sample_id), path("${sample_id}_AR_SJ.out.tab")                               , emit: spl_junc_tab

    tuple val(sample_id), path("${sample_id}_AR_Aligned.sortedByCoord.out.bam")         , emit: bam_sorted_aligned
    tuple val(sample_id), path("${sample_id}_AR_Aligned.sortedByCoord.out.bam.bai")     , optional:true, emit: bai_sorted_aligned

    script:

    """
    #!/bin/bash

    header=\$(zcat ${read1} | head -n 1 || true)
    id=\$(echo "\$header" | cut -f 3-4 -d":" | sed 's/:/./g')
    sample=\$(echo ${read1} | cut -f 1,2,3 -d"_")
    lane=\$(echo "\$header" | cut -f 4 -d":")
    flowcell=\$(echo "\$header" | cut -f 3 -d":" | sed 's/@//' | sed 's/:/_/g')
    barcode=\$(echo "\$header" | cut -f 10 -d":")

    STAR --outSAMattrRGline ID:\$id SM:${sample_id} PL:ILLUMINA PU:\$flowcell.\$lane.\$barcode \\
        --runThreadN 10 \\
        --genomeDir $AR_index \\
        --readFilesIn ${read1} ${read2} \\
        --readFilesCommand zcat \\
        --twopassMode Basic \\
        --outFileNamePrefix ${sample_id}/ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outReadsUnmapped Fastx \\
        --outFilterMultimapNmax 1 \\
        --alignEndsType EndToEnd \\
        --alignSoftClipAtReferenceEnds No \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --quantMode GeneCounts \\
        --peOverlapNbasesMin 10 \\
        --scoreGap -1000000 \\
        --limitBAMsortRAM 1018749472

    # Rename STAR outputs to include sample_id prefix
    for f in ${sample_id}/*; do
        fname=\$(basename \$f)
        mv "\$f" "${sample_id}_AR_\$fname"
    done

    samtools index ${sample_id}_AR_Aligned.sortedByCoord.out.bam ${sample_id}_AR_Aligned.sortedByCoord.out.bam.bai
    """
}