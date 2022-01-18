FROM gaow/base-notebook
LABEL maintainer="Amanda Tsai<at3535@cumc.columbia.edu>"
RUN cd /tmp
USER root
RUN R --slave -e "install.packages('tidyverse')"
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
git-all \
libboost-all-dev \
libgsl-dev \
zlib1g 
RUN git clone https://github.com/xqwen/torus.git
RUN cd torus/src && make && cp ./torus /usr/local/bin/torus && chmod +x /usr/local/bin/torus
RUN git clone https://github.com/xqwen/dap.git
RUN cd dap/dap_src/ && make && cp ./dap-g /usr/local/bin/dap-g && chmod +x /usr/local/bin/dap-g
RUN git clone https://github.com/xqwen/fastenloc.git
RUN cd fastenloc/src/ && make && cp ./fastenloc /usr/local/bin/fastenloc && chmod +x /usr/local/bin/fastenloc
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh