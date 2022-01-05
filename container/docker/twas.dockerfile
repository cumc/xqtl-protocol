

FROM gaow/base-notebook

MAINTAINER Gao Wang <wang.gao@columbia.edu>

WORKDIR /tmp
USER root

# PLINK

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
    unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin && rm -rf /tmp/*

# GCTA
RUN wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip && \
    unzip gcta_1.93.2beta.zip && mv gcta_1.93.2beta/gcta64 /usr/local/bin && rm -rf /tmp/*

# GEMMA
RUN wget https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz && \
    zcat gemma-0.98.1-linux-static.gz > /usr/local/bin/gemma && \
    chmod +x /usr/local/bin/gemma && rm -rf /tmp/*

# FUSION
ENV P2R_VERSION d74be015e8f54d662b96c6c2a52a614746f9030d
RUN wget https://github.com/gabraham/plink2R/archive/${P2R_VERSION}.zip && \
    unzip ${P2R_VERSION}.zip && \
    R --slave -e "install.packages(c('optparse','RColorBrewer', 'glmnet', 'RcppEigen','devtools','BiocManager','data.table'))" && \
    R --slave -e "install.packages('plink2R-${P2R_VERSION}/plink2R/',repos=NULL)" && \
    R --slave -e "devtools::install_github('stephenslab/susieR', ref='cran')" && \
    R --slave -e "devtools::install_github('hadley/devtools', ref='cran')" && \
    rm -rf /tmp/*

RUN R --slave -e "for (p in c('abind','data.table', 'readr', 'dplyr', 'tibble','modelr','purrr')) if (!(p %in% rownames(installed.packages()))) install.packages(p, repos = 'http://cran.rstudio.com')"

ENV FUSION_VERSION 6fedd22b47f9dab6a790c7779467f4d40ae57704
#RUN cd /usr/local/bin && wget https://raw.githubusercontent.com/gusevlab/fusion_twas/6fedd22b47f9dab6a790c7779467f4d40ae57704/FUSION.compute_weights.R && \
COPY FUSION.compute_weights.R /usr/local/bin/FUSION.compute_weights.R
RUN sed -i '1 s/^/#!\/usr\/bin\/env Rscript\n/' /usr/local/bin/FUSION.compute_weights.R && \
    chmod +x /usr/local/bin/FUSION.compute_weights.R

COPY FUSION.assoc_test.R /usr/local/bin/FUSION.assoc_test.R
RUN sed -i '1 s/^/#!\/usr\/bin\/env Rscript\n/' /usr/local/bin/FUSION.assoc_test.R && \
    chmod +x /usr/local/bin/FUSION.assoc_test.R


# GIT

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get --assume-yes update
RUN apt-get --assume-yes install git
RUN apt-get --assume-yes install libgsl0-dev

#RUN /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
#RUN brew install gsl




# mashr
RUN R --slave -e "install.packages('remote')" && \
    R --slave -e "remotes::install_github('stephenslab/mashr')" && \ 
    R --slave -e "remotes::install_github('stephenslab/mmbr')" 



# flashier
RUN R --slave -e "devtools::install_github('stephenslab/flashr',build_vignettes = TRUE)"
RUN R --slave -e "devtools::install_github('willwerscheid/flashier')"

# UDR
RUN R --slave -e "devtools::install_github('stephenslab/udr',build_vignettes = TRUE)"

 
USER jovyan

# to build: docker build -t gaow/twas -f twas.dockerfile .


