# Base image
FROM ubuntu:18.04

# Maintainer and author
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"

#### Build code for installing R ####
# For SOS based xQTL workwlow implementation. With just set an environment for PEER-R
# Modifided from https://github.com/RTIInternational/biocloud_docker_tools/blob/master/peer/v1.3/Dockerfile

RUN useradd docker \
    && mkdir /home/docker \
    && chown docker:docker /home/docker \
    && addgroup docker staff

RUN echo 'deb http://mirror.math.princeton.edu/pub/ubuntu/ bionic main' >> /etc/apt/sources.list \
    && apt-get update \ 
    && apt-get install -y --no-install-recommends \
        dirmngr \
        software-properties-common \
        lsb-release \
        ed \
        less \
        locales \
        wget \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ENV R_BASE_VERSION 3.5.1-2bionic
ENV DEBIAN_FRONTEND noninteractive


RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cpp \
        libgirepository-1.0-1 \
        libglib2.0-0 \
        libelf1 \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libmpx0 \
        curl \
        perl-base \
        gpg-agent \
        libxrender1 \
        yum \
    && apt-get install -y --no-install-recommends libopenblas-base
# Install R 3.5.1

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
    && cp /etc/apt/sources.list /etc/apt/sources.list~ \
    && sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list

# RUN apt-get update 


RUN apt-get update

# RUN apt-get install -y --no-install-recommends r-base r-base-dev r-base-core \
#     && apt-get clean \
#     && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
ENV R_BASE_VERSION 3.5.1-2bionic


# RUN apt-get install -y r-base
RUN wget https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz



RUN apt-get install -y \
        libgcc-6-dev \
        libstdc++-6-dev \
        gcc-5 \
        g++-5 \
        cmake \
        swig \
        gfortran \
        libbz2-dev \
        xfonts-base \
        libxrender1 \ 
        libx11-dev \
        libx11-6 \
        libxt-dev \
        # libXt-devel \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 --slave /usr/bin/g++ g++ /usr/bin/g++-5

RUN wget http://tukaani.org/xz/xz-5.2.2.tar.gz\
    &&tar xzvf xz-5.2.2.tar.gz \
    && cd xz-5.2.2 \
    && ./configure \
    && make \
    && make install

RUN tar -xzf R-3.5.1.tar.gz \
    && rm R-3.5.1.tar.gz \
    && cd R-3.5.1/ \
    && ./configure --prefix=/usr/lib/R/ --with-blas --with-readline=no\
    && make \
    && make install \
    && ln -s /usr/lib/R/bin/R /bin/R \
    && ln -s /usr/lib/R/bin/Rscript /bin/Rscript \
    && cd .. \
    && rm -rf R-3.5.1

WORKDIR /



#### Build code for installing PEER ####    
# Install dependencies and R PEER package

RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz \
   && R CMD INSTALL R_peer_source_1.3.tgz \
   && rm R_peer_source_1.3.tgz 
