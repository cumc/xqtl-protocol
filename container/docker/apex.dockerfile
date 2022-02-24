FROM gaow/base-notebook
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
su -  root # USER root
RUN apt-get update \
&& apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN git clone https://github.com/corbinq/apex.git \
&& cd apex/bin \
&& gunzip apex_Linux_x86_64.gz \
&& mv apex_Linux_x86_64 /usr/local/bin/apex && chmod +x /usr/local/bin/apex
#Install bcftools, tabix, and bgzip
RUN cd /tmp && wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.12 && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh