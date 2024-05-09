# Using the Ubuntu base image
# FROM ubuntu:noble-20240114
FROM ensemblorg/ensembl-vep

USER root
# # Installing prerequisites
# RUN apt-get update \
#   && apt-get install -y --no-install-recommends python=3.10 build-essential=12.10ubuntu1 wget=1.21.4-1ubuntu3 \
#   zlib1g-dev=1:1.3.dfsg-3.1ubuntu2 autoconf=2.71-3 automake=1:1.16.5-1.3ubuntu1 \
#   libncurses-dev=6.4+20240113-1ubuntu1 libbz2-dev=1.0.8-5ubuntu1 liblzma-dev=5.6.1+really5.4.5-1 \
#   libssl-dev=3.0.13-0ubuntu3 libcurl4-gnutls-dev=8.5.0-2ubuntu9 \
#   && rm -rf /var/lib/apt/lists/*

# # Installing prerequisites
# RUN apt-get update && \
    # apt-get install -y --no-install-recommends python2 perl build-essential=12.10ubuntu1 wget=1.21.4-1ubuntu4 \
    # zlib1g-dev=1:1.3.dfsg-3.1ubuntu2 autoconf=2.71-3 automake=1:1.16.5-1.3ubuntu1 \
    # libncurses-dev=6.4+20240113-1ubuntu2 libbz2-dev=1.0.8-5.1 liblzma-dev=5.6.1+really5.4.5-1 \
    # libssl-dev=3.0.13-0ubuntu3 libcurl4-gnutls-dev=8.5.0-2ubuntu10.1 bedtools=2.31.1+dfsg-2 \
    # && rm -rf /var/lib/apt/lists/*

# # Installing prerequisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends python2 perl build-essential wget \
    zlib1g-dev autoconf automake libncurses-dev libbz2-dev liblzma-dev \
    libssl-dev libcurl4-gnutls-dev bedtools git \
    && rm -rf /var/lib/apt/lists/*


## Installing HTS_LIB for bgzip and tabix
WORKDIR /tmp
ENV HTSLIB_INSTALL_DIR=/opt/htslib
RUN wget -q --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.9.tar.bz2 && \
    rm /tmp/htslib-1.9.tar.bz2 && \
    cd /tmp/htslib-1.9 && \
    ./configure --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix


# Install bcftools
WORKDIR /
RUN wget -q --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xvf bcftools-1.19.tar.bz2 && rm bcftools-1.19.tar.bz2
RUN cd bcftools-1.19 && make
WORKDIR /
ENV PATH="${PATH}:/bcftools-1.19"

# Install Samtools
RUN wget -q --no-check-certificate https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && tar -jxf samtools-1.11.tar.bz2 && rm samtools-1.11.tar.bz2
WORKDIR /samtools-1.11
RUN ./configure && make && make install
WORKDIR /



# Install Python packages. svtools requires Python2
RUN pip2 install svtools==0.5.1

# GIT_SSL_NO_VERIFY=true git clone https://github.com/hall-lab/svtools.git svtools_test


# ## Install VEP 90
# WORKDIR /opt/vep
# RUN cpan install DBI && cpan install Module::Build.pm && cpan install JSON && cpan install Bio::EnsEMBL::Registry
# RUN GIT_SSL_NO_VERIFY=true git clone https://github.com/Ensembl/ensembl-vep.git
# WORKDIR /opt/vep/ensembl-vep
# # RUN git checkout postreleasefix/90
# RUN perl -MCPAN -e 'install DBD::mysql'
# RUN perl INSTALL.pl --NO_UPDATE -a p --PLUGINS
# RUN ln -s /opt/vep/ensembl-vep/vep /usr/bin/variant_effect_predictor.pl
# WORKDIR /
