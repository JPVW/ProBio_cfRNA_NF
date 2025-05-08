# Dockerfile: STAR_UMItools_SAMtools_BEDtools.dockerfile
FROM continuumio/miniconda3

LABEL maintainer="janpietervw@gmail.com"

RUN conda install -y -c bioconda -c conda-forge \
    star=2.7.6a \
    samtools=1.14 \
    bedtools=2.30.0 \
    umi_tools=1.0.1 \
    subread=2.0.0 && \
    conda clean -a