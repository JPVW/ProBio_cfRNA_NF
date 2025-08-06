# AR Burden nextflow pipeline
![AR_burden_pipeline_overview](https://github.com/user-attachments/assets/1f0cf577-1611-44ef-92fc-12ebb7093ef3)

## INTRODUCTION

<p>The AR Burden pipeline is an bioinformatics pipeline capable of analyzing total RNAseq data, starting from a samplesheet and FATSQ files. The main parts of the pipeline consist of pre-processing, trimming, alignment to the reference genome as well as to known AR splice variants, UMI based deduplication, extensive QC, and detection of genomic variants (based on GATK best practices).<p>

<p>At the moment the pipeline only works with libraries generated with the Takara SMARTer Stranded Total RNA-seq kit. Library preparation can be followed by in-situ hybridization capture for target genes.  After alignment the pipeline will provide count tables per sample following a 2-pass STAR alignment. Additionally AR isoform specific counts will be detected and published in a count table.<p>

1. Read QC (fastqc)
2. UMI extraction (UMI-tools)
3. Adapter and quality trimming (Trim Galore !)
4. Alignment STAR 1st pass (STAR)
5. Merge splice junctions 
6. Regenerate STAR genome (STAR)
7. Alignment STAR 2 pass (STAR)
8. Subsetting to AR genomic locus (SAMtools)
    - Convert bam to fastq 
9. Alignment STAR AR (STAR)
10. UMI deduplication (UMI-tools)
11. Quantification (SubRead)
12. GATK RNAseq variant calling *WIP*
    - SplitNCigarReads
    - Baserecalibrator
    - APPLYBQSR
    - Haplotypecaller
    - SNPEFF
    - Variant filtration
13. Extensive QC metrics
    - preseq
    - picard
    - RSeQC
    - GATK
1. MultiQC

## USAGE

First, prepare a samplesheet (.csv) in the data direcory with the following structure:

"sample_id,fastq_1,fastq_2"

To run the pipeline use: nextflow run main.nf -resume

Before running the pipeline provide the following:
- Indexed fasta file of genome + 92 ERCC spikes
- GTF/dict/bed/bed12 of genome
- STAR indexed genome
- VEP cache

The following files are provided:
- Target interval list of TWIST probes
- AR isoform targets (https://www.nature.com/articles/s41591-021-01244-6)

## CREDIT

This pipeline is written, developed and maintained by Jan Vanwelkenhuyzen (@https://github.com/JPVW).

