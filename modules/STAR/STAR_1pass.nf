#!bin/env/usr nextflow

process STAR_1PASS {

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.6a--0"
    publishDir "results/STAR_1pass/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path index_zip

    output:
    tuple val(sample_id), path("${sample_id}_Log.final.out")     , emit: log_final
    tuple val(sample_id), path("${sample_id}_Log.out")           , emit: log_out
    tuple val(sample_id), path("${sample_id}_Log.progress.out")  , emit: log_progress
    tuple val(sample_id), path("${sample_id}_SJ.out.tab")      , emit: spl_junc_tab

    script:
    // --runThreadN 10  needs to be added in final script, can't be used during testing
    """

    STAR    --genomeDir $index_zip \
            --readFilesIn ${read1} ${read2} \
            --readFilesCommand zcat \
            --outFileNamePrefix ${sample_id}/ \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignMatesGapMax 1250000 \
            --alignIntronMax 1250000 \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 \
            --chimNonchimScoreDropMin 10 \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            --chimOutJunctionFormat 1

    # Rename STAR outputs to include sample_id prefix
    for f in ${sample_id}/*; do
        fname=\$(basename \$f)
        mv "\$f" "${sample_id}_\$fname"
    done 
    """
}