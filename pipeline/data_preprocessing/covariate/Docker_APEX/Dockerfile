FROM  gaow/base-notebook

# Maintainer and author
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"

USER root

RUN apt-get update \
    && apt install -y --no-install-recommends  git-all  libboost-all-dev 

RUN git clone https://github.com/corbinq/apex.git \
    && cd apex/bin \
    && gunzip apex_Linux_x86_64.gz \
    && mv apex_Linux_x86_64 apex && chmod +x apex 

WORKDIR ./apex/bin
