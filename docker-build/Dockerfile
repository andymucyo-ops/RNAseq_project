#base image 
FROM condaforge/mambaforge:24.9.2-0

#container info 
LABEL maintainer="andymucyo@gmail.com"
LABEL description="R environment for DESeq2 analysis - R 4.4 + DESeq2 + clusterProfiler"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive

#set workdir in container 
WORKDIR /project

#import virtual environment config into container
COPY environment.yml .

#create virtual env in container based on environment.yml file
RUN conda update -n base conda && \
    mamba env create -f environment.yml && \
    mamba clean -afy && \
    rm -rf /opt/conda/pkgs/*

#activation of conda env 
ENV PATH="/opt/conda/envs/rnaseq-r/bin:${PATH}"

#tests 
RUN R -e "library(DESeq2); library(clusterProfiler); cat('Environnement pret\n')"

#default command when launching container
CMD ["/bin/bash"]
