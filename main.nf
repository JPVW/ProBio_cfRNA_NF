#!/usr/bin/env nextflow


log.info """\
    P R O B I O     C F R N A   N F     P I P E L I N E 
    ====================================================
    This pipeline is designed to process targeted RNA-seq data using Nextflow.
    """
    .stripIndent()


// Module INCLUDE statements

    // Preprocessing
    include { FASTQC } from './modules/fastqc.nf'
    include { EXTRACT_UMI } from './modules/extract_UMI.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'

    // Alignment
    include { STAR_1PASS } from './modules/STAR/STAR_1pass.nf'
    include { MERGE_SPLICE_JUNCTIONS } from './modules/merge_sjs.nf'
    include { STAR_GENOME } from './modules/STAR/STAR_genome.nf'
    include { STAR_2PASS } from './modules/STAR/STAR_2pass.nf'
    include { UMI_TOOLS_DEDUP } from './modules/umi_tools_deduplication.nf'
    include { SAMTOOLS_INDEX as SAMTOOLS_INDEX1} from './modules/SAMtools/samtools_index.nf'
    include { SAMTOOLS_INDEX as SAMTOOLS_INDEX2} from './modules/SAMtools/samtools_index.nf'
    include { FEATURECOUNTS } from './modules/STAR/featurecounts.nf'
    include { SAMTOOLS_SORT } from './modules/SAMtools/samtools_sort.nf'
    include { BEDTOOLS_INTERSECT } from './modules/BEDtools/bedtools_intersect.nf'
    include { BAMTOBED }    from './modules/BEDtools/bamtobed.nf'
    include { AR_SUBSET } from './modules/AR_subset.nf'
    include { STAR_AR } from './modules/STAR/STAR_AR.nf'
    include { STAR_AR_COUNTING } from './modules/STAR/STAR_AR_counting.nf'

    // QC metrics
    include { HSMETRICS } from './modules/QC_metrics/HSmetrics.nf'
    include { ON_TARGET } from './modules/QC_metrics/on_target.nf'
    include { RNASEQMETRICS } from './modules/QC_metrics/RnaSeqMetrics.nf'
    include { COLLECTALIGNMENTMETRICS } from './modules/QC_metrics/CollectAlignmentSummaryMetrics.nf'
    include { INNER_DISTANCE } from './modules/QC_metrics/RSeQC/inner_distance.nf'
    include { JUNCTION_ANNOTATION } from './modules/QC_metrics/RSeQC/junction_annotation.nf'
    include { JUNCTION_SATURATION } from './modules/QC_metrics/RSeQC/junction_saturation.nf'
    include { INFER_EXPERIMENT } from './modules/QC_metrics/RSeQC/infer_experiment.nf'
    include { PRESEQ } from './modules/QC_metrics/preseq.nf'

    // GATK4 RNA variant detection
    include { SPLITNCIGARREADS } from './modules/GATK/splitncigarreads.nf'
    include { BASERECALIBRATOR } from './modules/GATK/BaseRecalibrator.nf'
    include { APPLYBQSR } from './modules/GATK/ApplyBQSR.nf'
    include { HAPLOTYPECALLER } from './modules/GATK/haplotypecaller.nf'
    include { TABIX } from './modules/GATK/tabix.nf'
    include { MERGEVCFS } from './modules/GATK/mergevcfs.nf'
    include { VARIANTFILTRATION } from './modules/GATK/variantfiltration.nf'
    include { SNPEFF_DOWNLOAD } from './modules/SNPEFF/snpeff_download.nf'
    include { SNPEFF_SNPEFF } from './modules/SNPEFF/snpeff.nf'
    include { VEP } from './modules/VEP/vep.nf'
    
    // MultiQC
    include { MULTIQC } from './modules/multiqc.nf'

// Primary input 
params.input_csv = "data/paired-end.csv"
params.report_id = "Test_pipline_subsampled"

// Pipeline parameters
    // provide new indices when running on a new server
params.STAR_1pass_index_zip = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/GRCh37/star2.7.10b"
params.STAR_2pass_fasta = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/GRCh37.p13.genome_with_spikes.fa"
params.fasta_fai = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/GRCh37.p13.genome_with_spikes.fa.fai"
params.fasta_dict = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/GRCh37.p13.genome_with_spikes.dict"
params.STAR_2pass_gtf = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/gencode.v30lift37.annotation.with.spikes.gtf"
params.STAR_AR_index = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/GRCh37/AR_subset_VanAllen"
params.STAR_AR_gtf = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/manual_AR_isoforms_VanAllen.gtf"
params.HSmetrics_target = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/Target_AR_mini_V2.interval_list"
params.HSmetrics_bait = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/Target_AR_mini_V2.interval_list"
params.HSmetrics_refFlat = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/gencode.v30lift37.annotation.with.spikes.refFlat.txt.gz"
params.HSmetrics_bed = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/gencode.v30lift37.annotation.with.spikes.sorted.bed"
params.HSmetrics_bed12 = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/gencode.v30lift37.annotation.with.spikes.bed12"
params.known_vcf = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/renamed.af-only-gnomad.raw.sites.grch37.vcf.gz"
params.known_vcf_tbi = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/renamed.af-only-gnomad.raw.sites.grch37.vcf.gz.tbi"

def db_name = "GRCh37.75"
def species = "homo_sapiens"
def assembly = "GRCh37"
def version = "114"
def cache = "/groups/wyattgrp/users/jvanwelkenhuyzen/pipelines/resources/RNAseq/vep_cache/homo_sapiens" // to change when running on new server

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
    STAR_1PASS(TRIM_GALORE.out.trimmed_reads, 
                file(params.STAR_1pass_index_zip))

    // Combine SJ output from STAR pass 1 into 1 file
    STAR_1PASS.out.spl_junc_tab.map { id, file -> file }
        .collect()
        .set { all_sj_files }

    MERGE_SPLICE_JUNCTIONS(all_sj_files)

    // Generate new index of the genome for STAR to get better alignments around novel splice junctions
    STAR_GENOME(file(params.STAR_2pass_fasta), 
                file(params.STAR_2pass_gtf), 
                MERGE_SPLICE_JUNCTIONS.out.merged_sjs)

    // Run 2nd round of STAR mapping and run SAMtools on bam output and after BEDtools intersect
    STAR_2PASS(TRIM_GALORE.out.trimmed_reads, 
                STAR_GENOME.out.genomedir, 
                file(params.STAR_2pass_gtf))

    BEDTOOLS_INTERSECT(STAR_2PASS.out.bam_sorted_aligned, 
                file(params.STAR_2pass_gtf))
    
    ch_star2_sorted_bai = SAMTOOLS_INDEX1(STAR_2PASS.out.bam_sorted_aligned) 
    ch_star2_exon_bai   = SAMTOOLS_INDEX2(BEDTOOLS_INTERSECT.out.junction_bam)

    
    // Rsubread to count genes after 2nd round of STAR
    UMI_TOOLS_DEDUP(ch_star2_sorted_bai, 
                    ch_star2_exon_bai, 
                    file(params.STAR_2pass_gtf))

    ch_dedup_bam = UMI_TOOLS_DEDUP.out.dedup_bam.join(UMI_TOOLS_DEDUP.out.dedup_bai, by: [0])
    FEATURECOUNTS(UMI_TOOLS_DEDUP.out.dedup_bam, 
                    UMI_TOOLS_DEDUP.out.dedup_exon_bam, 
                    UMI_TOOLS_DEDUP.out.dedup_unique_bam, 
                    file(params.STAR_2pass_gtf))

    
    // Subset reads mapping to the AR locus 'chrX:66753830-67011796' 
    AR_SUBSET(ch_star2_sorted_bai)

    // Remap AR subset reads to specific AR-V fasta file
    STAR_AR(AR_SUBSET.out.fq1, 
            AR_SUBSET.out.fq2, 
            file(params.STAR_AR_index), 
            file(params.STAR_AR_gtf))

    // Count AR-V isoforms
    STAR_AR_COUNTING(STAR_AR.out.bam_sorted_aligned, 
                    STAR_AR.out.bai_sorted_aligned.map { id, f -> f}, 
                    file(params.STAR_AR_gtf))

    // QC METRICS
    HSMETRICS(ch_dedup_bam, 
                file(params.STAR_2pass_fasta), 
                file(params.fasta_fai), 
                file(params.HSmetrics_target))
    //SAMTOOLS_INDEX2(HSMETRICS.out.dupremoved)
    ON_TARGET(ch_dedup_bam, 
                file(params.STAR_2pass_fasta), 
                file(params.fasta_fai), 
                file(params.fasta_dict), 
                file(params.HSmetrics_target))
    RNASEQMETRICS(ch_dedup_bam, 
                    file(params.HSmetrics_refFlat))
    COLLECTALIGNMENTMETRICS(ch_dedup_bam, 
                            file(params.STAR_2pass_fasta))
    INNER_DISTANCE(ch_dedup_bam, 
                    file(params.HSmetrics_bed), 
                    file(params.HSmetrics_bed12))
    INFER_EXPERIMENT(ch_dedup_bam, 
                    file(params.HSmetrics_bed), 
                    file(params.HSmetrics_bed12))
    JUNCTION_ANNOTATION(ch_dedup_bam, 
                        file(params.HSmetrics_bed), 
                        file(params.HSmetrics_bed12))
    JUNCTION_SATURATION(ch_dedup_bam, 
                        file(params.HSmetrics_bed), 
                        file(params.HSmetrics_bed12))
    BAMTOBED(ch_star2_sorted_bai)
    PRESEQ(BAMTOBED.out.sorted_bed)


    // RNAseq Variant calling GATK pipeline
    SPLITNCIGARREADS(ch_dedup_bam, 
                        file(params.STAR_2pass_fasta), 
                        file(params.fasta_fai), 
                        file(params.fasta_dict))

    BASERECALIBRATOR(SPLITNCIGARREADS.out.split_bam, 
                        file(params.STAR_2pass_fasta), 
                        file(params.fasta_fai), 
                        file(params.fasta_dict), 
                        file(params.known_vcf), 
                        file(params.known_vcf_tbi))
    APPLYBQSR(BASERECALIBRATOR.out.table, 
                file(params.STAR_2pass_fasta), 
                file(params.fasta_fai), 
                file(params.fasta_dict), 
                file(params.known_vcf), 
                file(params.known_vcf_tbi))
    HAPLOTYPECALLER(APPLYBQSR.out.BQSR_bam, 
                    file(params.STAR_2pass_fasta), 
                    file(params.fasta_fai), 
                    file(params.fasta_dict), 
                    file(params.known_vcf), 
                    file(params.known_vcf_tbi))
    //TABIX(HAPLOTYPECALLER.out.vcf)
    VARIANTFILTRATION(HAPLOTYPECALLER.out.vcf
                    .join(HAPLOTYPECALLER.out.tbi, failOnMismatch:true, failOnDuplicate:true), 
                    file(params.STAR_2pass_fasta), 
                    file(params.fasta_fai), 
                    file(params.fasta_dict))

    SNPEFF_DOWNLOAD(db_name)
    SNPEFF_SNPEFF(HAPLOTYPECALLER.out.vcf, 
                    db_name, 
                    SNPEFF_DOWNLOAD.out.cache)


    //VEP(HAPLOTYPECALLER.out.vcf, assembly, species, version, cache)


    // Comprehensive QC report
    MULTIQC(
        FASTQC.out.zip
            .map { id, f -> f }
            .mix(
            FASTQC.out.html.map { id, f -> f },
            EXTRACT_UMI.out.log,
            TRIM_GALORE.out.trimming_reports.map { id, f -> f },
            TRIM_GALORE.out.fastqc_reports_1.map { id, f -> f },
            TRIM_GALORE.out.fastqc_reports_2.map { id, f -> f },
            STAR_1PASS.out.log_final.map { id, f -> f },
            STAR_2PASS.out.log_final.map { id, f -> f },
            UMI_TOOLS_DEDUP.out.dedup_log.map { id, f -> f },
            AR_SUBSET.out.AR_log.map { id, f -> f },
            STAR_AR.out.log_final.map { id, f -> f },
            STAR_AR_COUNTING.out.AR_log.map { id, f -> f },
            HSMETRICS.out.metrics.map { id, f -> f },
            ON_TARGET.out.on_target_metrics.map { id, f -> f },
            RNASEQMETRICS.out.RNA_metrics.map { id, f -> f },
            COLLECTALIGNMENTMETRICS.out.CollectAlignmentSummary_metrics.map { id, f -> f },
            INNER_DISTANCE.out.pdf.map { id, f -> f },
            INFER_EXPERIMENT.out.txt.map { id, f -> f },
            JUNCTION_ANNOTATION.out.log.map { id, f -> f },
            JUNCTION_SATURATION.out.distance.map { id, f -> f },
            SPLITNCIGARREADS.out.split_bam.map { id, f -> f },
            BASERECALIBRATOR.out.table.map { id, bam, f -> f },
            APPLYBQSR.out.BQSR_bam.map { id, f -> f },
            HAPLOTYPECALLER.out.vcf.map { id, f -> f },
            VARIANTFILTRATION.out.vcf.map { id, f -> f },
            SNPEFF_SNPEFF.out.report.map { id, f -> f },
        ).collect(), params.report_id
    )
    
}