FROM ubuntu:18.04
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
ENV R_BASE_VERSION=3.5.1-2bionic
ENV R_BASE_VERSION=3.5.1-2bionic
RUN useradd docker
RUN mkdir /home/docker
RUN chown docker:docker /home/docker
RUN addgroup docker staff
RUN echo 'deb http://mirror.math.princeton.edu/pub/ubuntu/ bionic main' >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
dirmngr \
software-properties-common \
lsb-release \
ed \
less \
locales \
wget \
ca-certificates
RUN rm -rf /var/lib/apt/lists/*
RUN R_BASE_VERSION=3.5.1-2bionic
RUN apt-get update && \
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
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
yum
RUN apt-get install -y --no-install-recommends libopenblas-base
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN cp /etc/apt/sources.list /etc/apt/sources.list~
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list
RUN apt-get update
RUN R_BASE_VERSION=3.5.1-2bionic
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
libxt-dev
RUN rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 --slave /usr/bin/g++ g++ /usr/bin/g++-5
RUN wget http://tukaani.org/xz/xz-5.2.2.tar.gz
RUN tar xzvf xz-5.2.2.tar.gz
RUN cd xz-5.2.2
RUN ./configure
RUN make
RUN make install
RUN cd ..
RUN tar -xzf R-3.5.1.tar.gz
RUN rm R-3.5.1.tar.gz
RUN cd R-3.5.1/
RUN ./configure --prefix=/usr/lib/R/ --with-blas --with-readline=no
RUN make
RUN make install
RUN ln -s /usr/lib/R/bin/R /bin/R
RUN ln -s /usr/lib/R/bin/Rscript /bin/Rscript
RUN cd ..
RUN rm -rf R-3.5.1
RUN cd /
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz
RUN R CMD INSTALL R_peer_source_1.3.tgz
RUN rm R_peer_source_1.3.tgz
RUN echo "cd /" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh
