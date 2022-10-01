FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun
USER root

# upgrade packages to get R 4.2
RUN apt-get update && apt-get dist-upgrade -y &&  \
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
git-all \
libboost-all-dev \
libgsl-dev \ 
gfortran \
libreadline6-dev \
libx11-dev \
libxt-dev \
libpng-dev \
libjpeg-dev \
libcairo2-dev \
xvfb \
libbz2-dev \
libzstd-dev \
liblzma-dev \
libcurl4-openssl-dev \
texinfo \
texlive \
texlive-fonts-extra \
screen \
libpcre2-dev

RUN R --slave -e "update.packages(ask = FALSE, checkBuilt = TRUE)"
RUN R --slave -e 'install.packages("BiocManager")'
RUN R --slave -e 'BiocManager::install("psichomics")'
RUN wget https://raw.githubusercontent.com/cumc/xqtl-pipeline/main/code/psichomics_annotation_caching_location_update.R && \
Rscript psichomics_annotation_caching_location_update.R
## install sklearn and scipy
RUN pip install sklearn scipy
## SUPPA
RUN pip install SUPPA==2.3
RUN cd /opt && git clone https://github.com/comprna/SUPPA.git
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN echo 'alias SUPPA="python /opt/SUPPA/suppa.py"'  >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh
CMD exec /bin/bash "$@"
