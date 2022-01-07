FROM gaow/base-notebook AS spython-base
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
RUN cd /tmp
su -  root # USER root
RUN P2R_VERSION=d74be015e8f54d662b96c6c2a52a614746f9030d
RUN wget https://github.com/gabraham/plink2R/archive/${P2R_VERSION}.zip && \
unzip ${P2R_VERSION}.zip && \
R --slave -e "install.packages(c('rlang', 'RcppEigen','devtools','BiocManager','data.table'))" && \
R --slave -e "install.packages('plink2R-${P2R_VERSION}/plink2R/',repos=NULL)" && \
R --slave -e "devtools::install_github('stephenslab/susieR', ref='cran')" && \
R --slave -e "devtools::install_github('hadley/devtools', ref='cran')" &&
RUN R --slave -e "for (p in c('abind','data.table', 'readr', 'dplyr', 'tibble','modelr','purrr')) if (!(p %in% rownames(installed.packages()))) install.packages(p, repos = 'http://cran.rstudio.com')"
RUN R --slave -e "BiocManager::install('VariantAnnotation')"
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh