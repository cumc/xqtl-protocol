FROM gaow/base-notebook
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"
su -  root # USER root
RUN R --slave -e "install.packages(c('rlang',
                                     'tidyverse',
                                     'BiocManager', 
                                     'RcppEigen',
                                     # For kinship analysis
                                     'igraph'))"
                                     
RUN R --slave -e "BiocManager::install('biomaRt')"
RUN R --slave -e "BiocManager::install('VariantAnnotation')"

# Biopython package was used for summary stats merger script to handle strand flips
# QTL packages was used for Normalization of gene Count Table and TPM in Phenotype Normalization modules
RUN pip install qtl Biopython
                                     
RUN cd /tmp && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
    unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin && rm -rf /tmp/plink*
RUN cd /tmp && wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20211217.zip && \
    unzip plink2_linux_avx2_20211217.zip && mv plink2 /usr/local/bin && rm -rf /tmp/plink2*
RUN cd /tmp && wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip && \
    unzip gcta_1.93.2beta.zip && mv gcta_1.93.2beta/gcta64 /usr/local/bin && rm -rf /tmp/gcta*
#Install bcftools, tabix, and bgzip
RUN cd /tmp && wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.12 && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools
    
RUN cd /tmp && wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -O htslib-1.12.tar.bz2 && \
    tar -xjvf htslib-1.12.tar.bz2 && \
    cd htslib-1.12 && \
    ./configure --prefix=/usr/local/bin && \
    make && \
    make install && \
    cp tabix bgzip htsfile /usr/local/bin && rm -rf /tmp/htslib*
    
#Instal SnpEff that contains SnpSift
RUN cd /tmp && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip &&  \
    mv snpEff /opt && rm -rf /tmp/snpEff*

CMD exec /bin/bash "$@"
