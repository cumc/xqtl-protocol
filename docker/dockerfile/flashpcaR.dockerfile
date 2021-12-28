FROM  gaow/base-notebook

# Maintainer and author
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"

# This docker image contains flashpcar and other packages that are used in the data_preprocessing/PCA stage of the xqtl pipeline

USER root

## Tidyverse and other common R packages
RUN R --slave -e "install.packages(c('tidyverse','gridExtra', 'matrixStats','devtools'))"

## flashpcaR
RUN R --slave -e 'devtools::install_github("gabraham/flashpca/flashpcaR")'

### PLINK

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
    unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin 