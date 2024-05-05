# Using the Ubuntu base image
FROM ubuntu:noble-20240114

# Adding labels for the GitHub Container Registry
LABEL org.opencontainers.image.title="manta"
LABEL org.opencontainers.image.description="Container image for the use of Manta: Structural Variant and indel Caller"
LABEL org.opencontainers.image.version="latest"
LABEL org.opencontainers.image.authors="sanaz.iitd@gmail.com"
LABEL org.opencontainers.image.source=https://github.com/Sanaz01/docker-library


RUN apt-get update \
  && apt-get install -y --no-install-recommends bzip2 build-essential=12.10ubuntu1 wget=1.21.4-1ubuntu4 \
  zlib1g-dev=1:1.3.dfsg-3.1ubuntu2 autoconf=2.71-3 automake=1:1.16.5-1.3ubuntu1 \
  libncurses-dev=6.4+20240113-1ubuntu2 libbz2-dev=1.0.8-5.1 liblzma-dev=5.6.1+really5.4.5-1 \
  libssl-dev=3.0.13-0ubuntu3 libcurl4-gnutls-dev=8.5.0-2ubuntu10.1 bedtools=2.31.1+dfsg-2 \
  && rm -rf /var/lib/apt/lists/*


## Installing HTS_LIB for bgzip and tabix
ENV HTSLIB_INSTALL_DIR=/opt/htslib
WORKDIR /tmp
RUN wget -q --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.9.tar.bz2 && \
    rm /tmp/htslib-1.9.tar.bz2 && \
    cd /tmp/htslib-1.9 && \
    ./configure --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

## Install Manta
ENV manta_version 1.6.0
RUN wget -q --no-check-certificate https://github.com/Illumina/manta/releases/download/v${manta_version}/manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    tar -jxvf manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    rm manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    mv manta-${manta_version}.centos6_x86_64 /usr/local/src/manta

## Install Python2 to run Manta scripts (new versions of Ubuntu do not have Python2)
RUN wget -q --no-check-certificate https://www.python.org/ftp/python/2.7.9/Python-2.7.9.tgz && \
    tar -xzf Python-2.7.9.tgz && \
    rm Python-2.7.9.tgz && \
    cd Python-2.7.9 && \
    ./configure --enable-optimizations && \
    make altinstall

WORKDIR /