FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
USER root
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libgsl0-dev

RUN cd /tmp && wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -O htslib-1.12.tar.bz2 && \
tar -xjvf htslib-1.12.tar.bz2 && \
cd htslib-1.12 && \
./configure --prefix=/usr/local/bin && \
make && \
make install && \
cp tabix bgzip htsfile /usr/local/bin && cd ../ && rm -rf /tmp/htslib*

RUN R --slave -e "for (p in c('abind','rlang', 'RcppEigen', 'BiocManager', 'tidyr', 'data.table', 'tibble','modelr','purrr')) if (!(p %in% rownames(installed.packages()))) install.packages(p, repos = 'http://cran.rstudio.com')"

RUN R --slave -e "BiocManager::install('VariantAnnotation')"
# https://github.com/gabraham/plink2R/issues/10
RUN R --slave -e "Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE); remotes::install_github('gabraham/plink2R', subdir='plink2R', ref='d74be015e8f54d662b96c6c2a52a614746f9030d')"
RUN R --slave -e "remotes::install_github('stephenslab/flashr')"
RUN R --slave -e "remotes::install_github('willwerscheid/flashier')"
RUN R --slave -e "remotes::install_github('stephenslab/mashr')"
RUN R --slave -e "remotes::install_github('stephenslab/udr')"
RUN R --slave -e "remotes::install_github('stephenslab/susieR', build_vignettes=FALSE)"
RUN R --slave -e "remotes::install_github('stephenslab/mvsusieR')"
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD exec /bin/bash /bin/bash /entrypoint.sh "$@"
