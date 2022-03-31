FROM ubuntu:18.04
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
RUN useradd docker
RUN mkdir /home/docker
RUN chown docker:docker /home/docker
RUN addgroup docker staff
RUN echo 'deb http://mirror.math.princeton.edu/pub/ubuntu/ bionic main' >> /etc/apt/sources.list
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
	yum \
	dirmngr \
	software-properties-common \
	lsb-release \
	ed \
	less \
	locales \
	wget \
	ca-certificates \
	libopenblas-base \
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
	libxt-dev && rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 --slave /usr/bin/g++ g++ /usr/bin/g++-5
RUN cd /tmp && wget https://cran.r-project.org/src/base/R-3/R-3.5.1.tar.gz && wget http://tukaani.org/xz/xz-5.2.2.tar.gz
RUN tar xzvf xz-5.2.2.tar.gz && cd xz-5.2.2 && ./configure && make && make install && cd .. && rm -rf xz-5.2.2*
RUN tar -xzf R-3.5.1.tar.gz && cd R-3.5.1 && ./configure --prefix=/usr/lib/R/ --with-blas --with-readline=no && \
	make && make install && ln -s /usr/lib/R/bin/R /bin/R && ln -s /usr/lib/R/bin/Rscript /bin/Rscript && \
	cd .. && rm -rf R-3.5.1*
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && R CMD INSTALL R_peer_source_1.3.tgz && rm -f R_peer_source_1.3.tgz
RUN R --slave -e 'install.packages(c("dplyr", "tibble", "readr", "modelr", "purrr"), repos="http://cran.rstudio.com/")'
RUN R --slave -e 'devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))'
RUN pip install mofapy2
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh
