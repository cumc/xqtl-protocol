FROM gaow/base-notebook AS spython-base
LABEL maintainer="Amanda Tsai<at3535@cumc.columbia.edu>"
RUN cd /tmp
su -  root # USER root
RUN R -e "remotes::install_github('stephenslab/susieR',build_vignettes=FALSE)"
RUN R -e "install.packages('Ckmeans.1d.dp', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('KRIS', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN apt-get update \
&& apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN git clone https://github.com/omerwe/polyfun
RUN mv ./polyfun/* /usr/local/bin/
RUN wget http://www.christianbenner.com/finemap_v1.4_x86_64.tgz
RUN tar xvf finemap_v1.4_x86_64.tgz
RUN mv ./finemap_v1.4_x86_64/finemap_v1.4_x86_64 /usr/local/bin/
RUN wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
RUN tar xzvf ldstore_v2.0_x86_64.tgz
RUN mv ./ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 /usr/local/bin/
RUN python -m pip install -U pip
RUN pip install pandas-plink networkx scipy bitarray rpy2 sklearn
RUN chmod 777 /usr/local/bin/*
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh