FROM rocker/r-ver:4.2.0
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libxml2 \
libxt6 \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
libpcre3-dev \
libicu-dev \
libjpeg-dev \
libpng-dev \
libxml2-dev \
libglpk-dev \
libssl-dev \
libcurl4-gnutls-dev \
r-cran-rmysql \
wget \
libopenblas-base \
libatlas3-base
RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3
RUN cd /opt && \
wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
tar -xf htslib-1.11.tar.bz2 && rm htslib-1.11.tar.bz2 && cd htslib-1.11 && \
./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
make && make install && make clean
RUN cd /opt && \
wget https://raw.githubusercontent.com/hsun3163/xqtl-pipeline/main/data/cross_reactive_probe_Hop2020.txt
RUN R -e "install.packages('data.table', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('glmnet', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager')"
RUN R -e 'BiocManager::install("zwdzwd/sesame")'
RUN R -e 'sesameData::sesameDataCache()'
RUN R -e 'ExperimentHub::ExperimentHub()[["EH3675"]]' ## To retrieve the missing HM450 data
RUN R -e 'BiocManager::install("preprocessCore", configure.args="--disable-threading")'
RUN R -e 'BiocManager::install("minfi")'
RUN R -e 'BiocManager::install("limma")'
RUN R -e 'BiocManager::install("IlluminaHumanMethylationEPICmanifest")'
RUN R -e 'BiocManager::install("IlluminaHumanMethylation450kmanifest")'
RUN R -e 'BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")'
RUN R -e "BiocManager::install('achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38')"
CMD exec /bin/bash exec sh "$@"
