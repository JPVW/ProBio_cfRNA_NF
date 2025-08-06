#!/bin/bash
#SBATCH --job-name=Probio_cfRNA_NF
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

module load Nextflow/24.10.2

nextflow run /kyukon/scratch/gent/vo/002/gvo00224/TOBI/Projects/AR-burden/AR-Burden/ProBio_cfRNA_NF/main.nf -profile univ_hpc -resume -with-report
