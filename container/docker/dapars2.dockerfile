FROM gaow/base-notebook:latest

LABEL maintainer = "Liucheng Shi<ls3751@cumc.columbia.edu>" 

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

USER root
RUN conda install -y -c bioconda/label/cf201901 bedtools
RUN sudo apt update
RUN sudo apt-get update
# RUN sudo apt upgrade -y --fix-missing
RUN sudo apt-get install -y python2
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python2 get-pip.py
RUN python2.7 -m pip install numpy
RUN python2.7 -m pip install scipy
RUN R -e "install.packages(c('dplyr', 'tidyr', 'doParallel'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN conda install -y -c conda-forge r-vim
RUN conda install -y -c bioconda bioconductor-preprocesscore
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh