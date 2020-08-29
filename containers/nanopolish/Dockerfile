# Based on https://github.com/jts/nanopolish/blob/master/Dockerfile
# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/variant_callers_pipeline_nanopolish:pandora_paper_tag1

# This container has subsampling software to be used by the pandora_analysis_pipeline
# to build: sudo docker build . -t leandroishilima/variant_callers_pipeline_nanopolish:pandora_paper_tag1
FROM continuumio/miniconda3:4.8.2

# install dependencies
RUN conda install -c conda-forge -c bioconda samtools=1.10 minimap2=2.17 parallel=20200722

# install nanopolish
RUN conda install -c conda-forge -c bioconda nanopolish=0.13.2
