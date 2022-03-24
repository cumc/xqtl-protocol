FROM gaow/base-notebook:latest

LABEL maintainer = "Liucheng Shi<ls3751@cumc.columbia.edu>" 

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

USER root
RUN apt-get update && apt-get install -y python2 && \
    apt-get install -y build-essential wget \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python2 get-pip.py
RUN python2.7 -m pip install numpy scipy
RUN R -e "install.packages(c('doParallel', 'BiocManager'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e 'BiocManager::install("impute")'
RUN conda install -y -c bioconda/label/cf201901 bedtools
WORKDIR /usr/src
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	tar jxf samtools-1.10.tar.bz2 && \
	rm samtools-1.10.tar.bz2 && \
	cd samtools-1.10 && \
	./configure --prefix $(pwd) && \
	make
ENV PATH=${PATH}:/usr/src/samtools-1.10 
# Dapars scripts
RUN for i in DaPars_Extract_Anno.py Dapars2_Multi_Sample.py gtf2bed12.py; do wget https://raw.githubusercontent.com/cumc/xqtl-pipeline/main/code/${i} \
	&& mv ${i} /usr/local/bin/ && chmod +x /usr/local/bin/${i}; done