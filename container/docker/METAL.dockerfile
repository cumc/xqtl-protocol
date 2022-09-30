FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
RUN cd /tmp
USER root
RUN apt-get update \
&& apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN git clone https://github.com/statgen/METAL
RUN cd METAL
RUN mkdir build && cd build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make
RUN make test
RUN make install
RUN cp ./bin/metal  /usr/local/bin/metal && chmod +x /usr/local/bin/metal
RUN R --slave -e "install.packages(c('BiocManager'))"
RUN R --slave -e "BiocManager::install('VariantAnnotation')"
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh