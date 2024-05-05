FROM ubuntu:noble-20240114
# FROM continuumio/miniconda3:24.1.2-0


# Installing prerequisites
RUN apt-get update \
  && apt-get install -y --no-install-recommends python3 build-essential=12.10ubuntu1 bedtools=2.31.1+dfsg-2 wget \
  && rm -rf /var/lib/apt/lists/*

# Configure environment
# ENV CONDA_DIR /opt/conda
# ENV PATH $CONDA_DIR/bin:$PATH

#originally installed 
#https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh, which is 
#https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh
#this is python 3.8, but cyvcf2. a necessary package for make_report[3].py, only has builds
#switching to current for python 3.7

# RUN cd /tmp && \
#     mkdir -p $CONDA_DIR && \
#     curl -s https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -o miniconda.sh && \
#     /bin/bash miniconda.sh -f -b -p $CONDA_DIR && \
#     rm miniconda.sh && \
#     $CONDA_DIR/bin/conda config --system --add channels conda-forge && \
#     $CONDA_DIR/bin/conda config --system --set auto_update_conda false && \
#     conda clean -tipsy

# # Switching from hard-links to copies
# RUN conda config --set always_copy true

# RUN conda config --add channels bioconda && \
#     conda install -y -c bioconda mosdepth


# wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh
# bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -f -b -p /opt/conda

WORKDIR /tmp/mosdepth
RUN wget --no-check-certificate https://conda.anaconda.org/bioconda/linux-64/mosdepth-0.3.8-hd299d5a_0.tar.bz2 && \
    tar --bzip2 -xvf mosdepth-0.3.8-hd299d5a_0.tar.bz2 && \
    mv /tmp/mosdepth/bin/mosdepth /usr/bin && \
    rm -r /tmp/mosdepth/

WORKDIR /