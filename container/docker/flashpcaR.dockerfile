FROM gaow/base-notebook
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"
su -  root # USER root
RUN R --slave -e "install.packages(c('rlang', 
                                     'tidyverse',
                                     'gridExtra', 
                                     'matrixStats'))"
RUN R --slave -e 'remotes::install_github("gabraham/flashpca/flashpcaR")'
RUN cd /tmp && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
    unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin && rm -rf /tmp/*
CMD exec /bin/bash "$@"