FROM  gaow/base-notebook

# Maintainer and author
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"

# This docker image contains basic bioinfo tools that are used in the data_preprocessing stage of the xqtl pipeline

USER root

## Tidyverse
RUN R --slave -e "install.packages(c('tidyverse','BiocManager', 'RcppEigen'))"

# Genotype

## GWAS_QC

### PLINK

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
    unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin 
    
### PLINK2

RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20211217.zip && \
    unzip plink2_linux_avx2_20211217.zip && mv plink2 /usr/local/bin 
    
### igraph
RUN R --slave -e "install.packages('igraph')"

## GRM
### GCTA

RUN wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip && \
    unzip gcta_1.93.2beta.zip && mv gcta_1.93.2beta/gcta64 /usr/local/bin
## LD
### PLINK (satisfied)

# Phenotype
## Normalization
### QTL packages

RUN pip install qtl

### BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2\
    && tar xvf bcftools-1.2.tar.bz2 \
    && cd bcftools-1.2 \
    && make\
    && make install\
    && cd htslib-1.2.1\
    && make\
    && make install

RUN rm -rf bcftools-1.2 \
    && rm -rf bcftools-1.2.tar.bz2
 
## annotation
RUN R --slave -e "BiocManager::install('biomaRt')"

## Phenotype_formatting & Residual
### bcftools(satisfied)

# Covariate
## See container APEX, PEER, and flashpca

# clear place 
rm -rf /tmp/*