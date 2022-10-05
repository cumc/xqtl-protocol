FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun
USER root
ENV PATH=/opt/samtools-1.11:$PATH
ENV PATH=/opt/htslib-1.11:$PATH
ENV PATH=/opt/regtools/build:$PATH
ENV PATH=/opt/leafcutter/clustering:$PATH
RUN apt-get update && \
apt-get install -y --no-install-recommends \
git-all \
libboost-all-dev \
libgsl-dev
RUN pip install sklearn scipy
RUN R --slave -e 'install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE, versions = "2.21.1")'
RUN R --slave -e 'install.packages("BiocManager")'
RUN R --slave -e 'BiocManager::install(c("Biobase","DirichletMultinomial"))'
RUN R --slave -e 'remotes::install_github("davidaknowles/leafcutter/leafcutter@63b347a316cc214808b8c734ba181c602e950f06")'
RUN cd /opt && \
git clone https://github.com/griffithlab/regtools && \
cd regtools/ && \
git reset --hard 3d007f5698bb97745aa6e2f96f8aed0e6d3c9d0e && \
mkdir build && \
cd build/ && \
cmake .. && \
make
RUN cd /opt && \
git clone https://github.com/davidaknowles/leafcutter && \
cd leafcutter/  && \
git reset --hard 63b347a316cc214808b8c734ba181c602e950f06
RUN cd /opt && \
wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
tar -xf htslib-1.11.tar.bz2 && rm htslib-1.11.tar.bz2 && cd htslib-1.11 && \
./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
make && make install && make clean
RUN cd /opt && \
wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
tar -xf samtools-1.11.tar.bz2 && rm samtools-1.11.tar.bz2 && cd samtools-1.11 && \
./configure --with-htslib=/opt/htslib-1.11 && make && make install && make clean
CMD exec /bin/bash "$@"
