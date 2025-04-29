#!/usr/bin/env nextflow

/* Steps overview we need to build:
 * Fastqc + trimming
 * Extract UMI
 * Trim UMI
 * STAR_pass1
 * SJ merge and STAR genome
 * STAR pass 2
 * STAR pass 2 counting
 * AR subset
 * STAR AR 1 pass
 * STAR AR counting
 * HSmetrics
 * ocunt on target reads
 * RNAseqmetrics
 * Collect Alignment summay metrics
 * inner distance
 * preseq
 * Markduplicates etc for variant calling in RNAseq 
 * multiqc
 */ 

log.info """\
    P R O B I O     C F R N A   N F     P I P E L I N E 
    ====================================================
    This pipeline is designed to process targeted RNA-seq data using Nextflow.
    reads: ${params.reads}
    """
    .stripIndent()


// Module INCLUDE statements

include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'

// Primary input 
params.reads = "$projectDir/data/RNA026749_S7_subsamp10000_s100_R{1,2}_001.fastq.gz"

workflow  {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch) 
    
}