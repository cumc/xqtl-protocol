FROM gaow/base-notebook:latest

LABEL maintainer = "Liucheng Shi<ls3751@cumc.columbia.edu>" 

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

USER root
RUN apt-get update && apt-get install -y python2 && apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN python2 get-pip.py
RUN python2.7 -m pip install numpy scipy
RUN R -e "install.packages(c('doParallel', 'VIM', 'BiocManager'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN Rscript -e 'BiocManager::install("preprocessCore")'
RUN conda install -y -c bioconda/label/cf201901 bedtools

# Dapars scripts
RUN for i in DaPars_Extract_Anno.py Dapars2_Multi_Sample.py gtf2bed12.py; do wget https://raw.githubusercontent.com/cumc/xqtl-pipeline/main/code/${i} \
	&& mv ${i} /usr/local/bin/ && chmod +x /usr/local/bin/${i}; done
