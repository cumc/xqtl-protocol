FROM gaow/base-notebook

MAINTAINER Gao Wang <wang.gao@columbia.edu>

WORKDIR /tmp
USER root

# PLINK

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
    unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin && rm -rf /tmp/*


# R

RUN R --slave -e "for (p in c('rlang','abind','data.table', 'readr', 'dplyr', 'tibble','modelr','purrr')) if (!(p %in% rownames(installed.packages()))) install.packages(p, repos = 'http://cran.rstudio.com')"

# GIT

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get --assume-yes update
RUN apt-get --assume-yes install git
RUN apt-get --assume-yes install libgsl0-dev

#RUN /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
#RUN brew install gsl


# APEX
RUN git clone https://github.com/corbinq/apex.git && \
    cd apex/bin && \
    gunzip apex_Linux_x86_64.gz && \
    mv apex_Linux_x86_64 /usr/local/bin/apex && chmod +x /usr/local/bin/apex

RUN sudo apt-get install tabix


# TensorQTL

RUN pip install tensorqtl
    
USER jovyan

# to build: docker build -t hs3163/apex -f apex.dockerfile .
