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
RUN cd /tmp \
    && wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2 \
    && tar xvf bcftools-1.2.tar.bz2 \
    && cd bcftools-1.2 \
    && make\
    && make install \
    && cd htslib-1.2.1\
    && make \
    && make install \
    && rm -rf /tmp/bcftools*
    
# Normalization and collapse annotation script
RUN wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/gene_model/collapse_annotation.py && mv collapse_annotation.py /usr/local/bin/

RUN wget wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py && mv eqtl_prepare_expression.py /usr/local/bin/


# gffread
RUN apt-get update && apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN cd /tmp
RUN git clone https://github.com/gpertea/gffread
RUN cd gffread
RUN make release
RUN mv gffread  /usr/local/bin/
RUN cd .. && rm -r gffread

CMD exec /bin/bash "$@"
