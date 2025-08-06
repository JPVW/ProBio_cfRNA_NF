process MERGE_SPLICE_JUNCTIONS {
    label 'process_medium'
    
    publishDir "results", method: 'copy'

    input:
    path sj_files

    output:
    path("SJ.out.pass1_merged.tab"), emit: merged_sjs

    script:
    """
    cat ${sj_files.join(' ')} | awk '\$7 >= 3' | cut -f1-4 | sort -u > SJ.out.pass1_merged.tab
    """
}