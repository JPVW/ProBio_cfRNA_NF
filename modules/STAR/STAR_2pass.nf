process STAR_2PASS {
    label 'process_high'
    
    container 'docker://jpvw/probio_cfrna_tools:1.0'
    
    publishDir "results/STAR_2pass/${sample_id}/", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path new_index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}_Log.final.out")                , emit: log_final
    tuple val(sample_id), path("${sample_id}_Log.out")                      , emit: log_out
    tuple val(sample_id), path("${sample_id}_Log.progress.out")             , emit: log_progress
    tuple val(sample_id), path("${sample_id}_SJ.out.tab")                   , emit: spl_junc_tab

    tuple val(sample_id), path('*d.out.bam')                                     , optional:true, emit: bam
    tuple val(sample_id), path("${sample_id}_sortedByCoord.out.bam")             , optional:true, emit: bam_sorted
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")     , optional:true, emit: bam_sorted_aligned
    tuple val(sample_id), path('*.tab')                                          , optional:true, emit: tab
    tuple val(sample_id), path('*.ReadsPerGene.out.tab')                         , optional:true, emit: read_per_gene_tab
    tuple val(sample_id), path('*.out.junction')                                 , optional:true, emit: junction


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
        --genomeDir $new_index \\
        --readFilesIn ${read1} ${read2} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}/ \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFilterMultimapNmax 20 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --alignIntronMin 20 \\
        --alignMatesGapMax 1250000 \\
        --alignIntronMax 1250000 \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 12 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimMultimapScoreRange 3 \\
        --chimScoreJunctionNonGTAG -4 \\
        --chimMultimapNmax 20 \\
        --chimNonchimScoreDropMin 10 \\
        --peOverlapNbasesMin 12 \\
        --peOverlapMMp 0.1 \\
        --chimOutJunctionFormat 1 \\
        --quantMode GeneCounts

    # Rename STAR outputs to include sample_id prefix
    for f in ${sample_id}/*; do
        fname=\$(basename \$f)
        mv "\$f" "${sample_id}_\$fname"
    done
    """
}