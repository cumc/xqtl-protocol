FROM gaow/base-notebook
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
RUN su -  root # USER root
RUN apt-get update \
&& apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN git clone https://github.com/corbinq/apex.git \
&& cd apex/bin \
&& gunzip apex_Linux_x86_64.gz \
&& mv apex_Linux_x86_64 /usr/local/bin/apex && chmod +x /usr/local/bin/apex
#Install bcftools, tabix, and bgzip
RUN cd /tmp && wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -O htslib-1.12.tar.bz2 && \
    tar -xjvf htslib-1.12.tar.bz2 && \
    cd htslib-1.12 && \
    ./configure --prefix=/usr/local/bin && \
    make && \
    make install && \
    cp tabix bgzip htsfile /usr/local/bin && rm -rf /tmp/htslib*
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh
