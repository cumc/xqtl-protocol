FROM gaow/base-notebook
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
USER root
RUN apt-get update \
&& apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN R --slave -e "install.packages(c('rlang', 'RcppEigen','devtools'))"
RUN R --slave -e 'devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))'
RUN pip install mofapy2
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh