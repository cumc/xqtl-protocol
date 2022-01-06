FROM gaow/base-notebook AS spython-base
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"
su -  root # USER root
RUN R --slave -e "install.packages(c('rlang',,tidyverse','gridExtra', 'matrixStats', 'devtools'))"
RUN R --slave -e 'devtools::install_github("gabraham/flashpca/flashpcaR")'
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin
CMD exec /bin/bash "$@"