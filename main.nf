#!/usr/bin/env nextflow

/* Steps overview we need to build:
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
    """
    .stripIndent()


// Module INCLUDE statements

include { FASTQC } from './modules/fastqc.nf'
include { EXTRACT_UMI } from './modules/extract_UMI.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { STAR_1PASS } from './modules/STAR/STAR_1pass.nf'
include { MERGE_SPLICE_JUNCTIONS } from './modules/merge_sjs.nf'
include { STAR_GENOME } from './modules/STAR/STAR_genome.nf'
include { STAR_2PASS } from './modules/STAR/STAR_2pass.nf'
include { STAR_2PASS_COUNTING } from './modules/STAR/STAR_2pass_counting.nf'
include { MULTIQC } from './modules/multiqc.nf'
include { SAMTOOLS_INDEX } from './modules/SAMtools/samtools_index.nf'
include { BEDTOOLS_INTERSECT } from './modules/BEDtools/bedtools_intersect.nf'


// Primary input 
params.input_csv = "/kyukon/scratch/gent/vo/002/gvo00224/TOBI/Projects/AR-burden/AR-Burden/ProBio_cfRNA_NF/data/paired-end.csv"
params.report_id = "Test_pipline_subsampled"

// Pipeline parameters
params.STAR_1pass_index_zip = "/kyukon/scratch/gent/vo/002/gvo00224/TOBI/Projects/AR-burden/AR-Burden/ARV_fastas/ARV_RNAseq/v0.2/GRCh37/GRCh37.p13_gencodev30_withspikes_STAR_2.7.10b"
params.STAR_2pass_fasta = "/kyukon/scratch/gent/vo/002/gvo00224/TOBI/Projects/AR-burden/AR-Burden/ARV_fastas/ARV_RNAseq/v0.2/GRCh37/GRCh37.p13.genome_with_spikes.fa"
params.STAR_2pass_gtf = "/kyukon/scratch/gent/vo/002/gvo00224/TOBI/Projects/AR-burden/AR-Burden/ARV_fastas/ARV_RNAseq/v0.2/GRCh37/gencode.v30lift37.annotation.with.spikes.gtf"

workflow  {

    // Create input channel from a file path
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [row.sample_id ,file(row.fastq_1), file(row.fastq_2)] }

    // Initial quality control
    FASTQC(read_ch)

    // UMI tag extraction
    EXTRACT_UMI(read_ch)

    // Adapter trimming
    TRIM_GALORE(EXTRACT_UMI.out.UMI_extracted)

    // First round of STAR mapping
    STAR_1PASS(TRIM_GALORE.out.trimmed_reads, file(params.STAR_1pass_index_zip))

    // Combine SJ output from STAR pass 1 into 1 file
    STAR_1PASS.out.spl_junc_tab.map { id, file -> file }
        .collect()
        .set { all_sj_files }

    MERGE_SPLICE_JUNCTIONS(all_sj_files)

    // Generate new index of the genome for STAR to get better alignments around novel splice junctions

    STAR_GENOME(file(params.STAR_2pass_fasta), file(params.STAR_2pass_gtf), MERGE_SPLICE_JUNCTIONS.out.merged_sjs)

    // Run 2nd round of STAR mapping and run SAMtools on bam output and after BEDtools intersect

    STAR_2PASS(TRIM_GALORE.out.trimmed_reads, STAR_GENOME.out.genomedir, file(params.STAR_2pass_gtf))
    BEDTOOLS_INTERSECT(STAR_2PASS.out.bam_sorted_aligned, file(params.STAR_2pass_gtf))
    SAMTOOLS_INDEX(STAR_2PASS.out.bam_sorted_aligned, BEDTOOLS_INTERSECT.out.junction_bam)

    // Rsubread to count genes after 2nd round of STAR
    STAR_2PASS_COUNTING(STAR_2PASS.out.bam_sorted_aligned,SAMTOOLS_INDEX.out.bai, BEDTOOLS_INTERSECT.out.junction_bam.map {id, f -> f}, SAMTOOLS_INDEX.out.exon_bai, file(params.STAR_2pass_gtf))
    

    // Comprehensive QC report
    MULTIQC(
        FASTQC.out.zip.mix(
            FASTQC.out.html,
            EXTRACT_UMI.out.log,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            STAR_1PASS.out.log_final.map { id, f -> f }
        ).collect(), params.report_id
    )
    
}