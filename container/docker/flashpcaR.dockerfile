FROM gaow/base-notebook
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"
su -  root # USER root
RUN R --slave -e "install.packages(c('rlang', 
                                     'tidyverse',
                                     'gridExtra', 
                                     'matrixStats'))"
RUN R --slave -e 'remotes::install_github("gabraham/flashpca/flashpcaR")'
CMD exec /bin/bash "$@"
