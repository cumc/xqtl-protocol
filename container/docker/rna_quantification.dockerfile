FROM gaow/base-notebook
MAINTAINER Francois Aguet; Hao Sun

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
	software-properties-common \
        libboost-all-dev \
        libbz2-dev \
        libhdf5-dev \
        libncurses5-dev \
        default-jdk && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
## FIXME: The 'rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*'  command was removed, please add back the specific contents that are to be removed. 


#-----------------------------
# Pipeline components
#-----------------------------

# R packages for TPM QC
RUN R --slave -e "install.packages(c('rlang','RcppEigen','RColorBrewer','ape','reshape2'))"

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
    tar -xf htslib-1.11.tar.bz2 && rm htslib-1.11.tar.bz2 && cd htslib-1.11 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar -xf samtools-1.11.tar.bz2 && rm samtools-1.11.tar.bz2 && cd samtools-1.11 && \
    ./configure --with-htslib=/opt/htslib-1.11 && make && make install && make clean

# STAR v2.7.8a
RUN cd /opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz && \
    tar -xf 2.7.8a.tar.gz && rm 2.7.8a.tar.gz
ENV PATH /opt/STAR-2.7.8a/bin/Linux_x86_64_static:$PATH

# RSEM v1.3.3
RUN cd /opt && \
    wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar -xvf v1.3.3.tar.gz && rm v1.3.3.tar.gz && cd RSEM-1.3.3 && make
ENV PATH /opt/RSEM-1.3.3:$PATH

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz && \
    tar -xf v2.4.1.tar.gz && rm v2.4.1.tar.gz && cd bamtools-2.4.1 && mkdir build && cd build && cmake .. && make && make install && make clean
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

# bamsync
ENV GTEX_PIPELINE 3481dd43b2e8a33cc217155483ce25d5255aafa9
RUN cd /tmp && \
    wget --no-check-certificate https://github.com/broadinstitute/gtex-pipeline/archive/${GTEX_PIPELINE}.zip -O gtex-pipeline.zip && unzip gtex-pipeline.zip && cd gtex-pipeline-${GTEX_PIPELINE} && cd rnaseq/bamsync && make && mv bamsync /usr/local/bin

# Picard tools
RUN mkdir /opt/picard-tools && \
    wget --no-check-certificate -P /opt/picard-tools/ https://github.com/broadinstitute/picard/releases/download/2.25.0/picard.jar

# kallisto
RUN cd /opt && \
    wget https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_linux-v0.46.2.tar.gz && \
    tar -xf kallisto_linux-v0.46.2.tar.gz && rm kallisto_linux-v0.46.2.tar.gz
ENV PATH $PATH:/opt/kallisto_linux-v0.46.2

# bedtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
    tar -xf bedtools-2.30.0.tar.gz && rm bedtools-2.30.0.tar.gz && \
    cd bedtools2 && make && make install && make clean

# UCSC tools
RUN mkdir /opt/ucsc && \
    wget --no-check-certificate -P /opt/ucsc/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph && \
    wget --no-check-certificate -P /opt/ucsc/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod 755 /opt/ucsc/*
ENV PATH /opt/ucsc:$PATH

# python modules
RUN curl https://bootstrap.pypa.io/get-pip.py | python3
RUN pip3 install --upgrade pip setuptools
RUN pip3 install tables numpy "pandas>=1.4.0" scipy pyarrow matplotlib seaborn
# numpy dependencies:
RUN pip3 install pyBigWig

# FastQC
RUN cd /opt && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && mv FastQC FastQC-0.11.9 && cd FastQC-0.11.9 && chmod 775 fastqc
ENV PATH /opt/FastQC-0.11.9:$PATH

# RNA-SeQC
RUN mkdir /opt/rnaseqc && cd /opt/rnaseqc && \
    wget https://github.com/getzlab/rnaseqc/releases/download/v2.4.2/rnaseqc.v2.4.2.linux.gz && \
    gunzip rnaseqc.v2.4.2.linux.gz && mv rnaseqc.v2.4.2.linux rnaseqc && chmod 775 rnaseqc
RUN pip3 install rnaseqc
ENV PATH /opt/rnaseqc:$PATH

# gcloud
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update -y && apt-get install google-cloud-sdk -y

# Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
RUN mv Trimmomatic-0.39.zip /opt/
RUN cd /opt/
RUN unzip Trimmomatic-0.39.zip

# scripts

RUN mv /tmp/gtex-pipeline-${GTEX_PIPELINE}/rnaseq/src /opt/ && rm -rf /tmp/gtex-pipeline*
RUN export PATH=/opt/src:$PATH

ENV PATH /opt/src:$PATH


# Normalization and collapse annotation script
RUN wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/gene_model/collapse_annotation.py && mv collapse_annotation.py /usr/local/bin/ && chmod +x /usr/local/bin/collapse_annotation.py

## The sed is a temporary patch for the issue discussed in https://github.com/pandas-dev/pandas/pull/44632#issuecomment-1029185339
RUN wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/1d6c11c10f9c2e0befeb7076ba9df3a5832d2a0b/qtl/src/eqtl_prepare_expression.py && sed "s/dtype=str/dtype={0:str,1:str}/g"  ./eqtl_prepare_expression.py > /usr/local/bin/eqtl_prepare_expression.py &&  chmod +x /usr/local/bin/eqtl_prepare_expression.py && rm ./eqtl_prepare_expression.py

# gffread
RUN cd /tmp
RUN wget https://github.com/gpertea/gffread/archive/refs/tags/v0.12.7.zip && \
    unzip -o v0.12.7.zip && \
    cd gffread-0.12.7 && make release && \
    mv gffread  /usr/local/bin/ && \
    cd .. && rm -r gffread* v0.12.7*
