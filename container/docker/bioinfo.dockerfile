FROM gaow/base-notebook
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"
su -  root # USER root
RUN R --slave -e "install.packages(c('rlang','tidyverse','BiocManager', 'RcppEigen'))"
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20211217.zip && \
unzip plink2_linux_avx2_20211217.zip && mv plink2 /usr/local/bin
RUN R --slave -e "install.packages('igraph')"
RUN wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip && \
unzip gcta_1.93.2beta.zip && mv gcta_1.93.2beta/gcta64 /usr/local/bin
RUN pip install qtl
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
RUN R --slave -e "BiocManager::install('biomaRt')"
RUN R --slave -e "BiocManager::install('VariantAnnotation')"
RUN pip install Biopython
CMD exec /bin/bash "$@"
