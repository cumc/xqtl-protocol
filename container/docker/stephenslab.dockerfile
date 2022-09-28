FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
RUN cd /tmp
su -  root # USER root
RUN R --slave -e "install.packages(c('rlang', 'RcppEigen','BiocManager','data.table'))"
RUN R --slave -e "install.packages('remotes')"
RUN R --slave -e "install.packages('tidyr')"
RUN R --slave -e "remotes::install_github('stephenslab/flashr')"
RUN R --slave -e "remotes::install_github('willwerscheid/flashier')"
RUN apt-get --assume-yes update
RUN apt-get --assume-yes install libgsl0-dev
RUN R --slave -e "remotes::install_github('stephenslab/mashr')"
RUN R --slave -e "remotes::install_github('stephenslab/mvsusieR')"
RUN P2R_VERSION=d74be015e8f54d662b96c6c2a52a614746f9030d
RUN cd /opt/
RUN wget https://github.com/gabraham/plink2R/archive/${P2R_VERSION}.zip && \
unzip -o ${P2R_VERSION}.zip && \
R --slave -e "install.packages('plink2R-${P2R_VERSION}/plink2R/',repos=NULL)" && \
R --slave -e "remotes::install_github('stephenslab/susieR')" && \
R --slave -e "remotes::install_github('hadley/devtools', ref='cran')" &&
RUN R --slave -e "for (p in c('abind','data.table', 'tibble','modelr','purrr')) if (!(p %in% rownames(installed.packages()))) install.packages(p, repos = 'http://cran.rstudio.com')"
RUN cd /tmp/
RUN R --slave -e "BiocManager::install('VariantAnnotation')"
RUN R --slave -e "remotes::install_github('stephenslab/udr')"
RUN cd /tmp && wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -O htslib-1.12.tar.bz2 && \
tar -xjvf htslib-1.12.tar.bz2 && \
cd htslib-1.12 && \
./configure --prefix=/usr/local/bin && \
make && \
make install && \
cp tabix bgzip htsfile /usr/local/bin && rm -rf /tmp/htslib*
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD exec /bin/bash /bin/bash /entrypoint.sh "$@"