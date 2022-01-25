FROM gaow/base-notebook:1.0.0

LABEL maintainer = "Liucheng Shi<ls3751@cumc.columbia.edu>" 

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

USER root

RUN conda install -y -c bioconda/label/cf201901 plink2
RUN conda install -y -c bioconda/label/cf201901 bcftools
RUN R -e "install.packages(c('dplyr', 'tidyr', 'doParallel','VIM','preprocessCore'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN mkdir /script
RUN cd /script
RUN wget https://raw.githubusercontent.com/cumc/xqtl-pipeline/main/code/DaPars_Extract_Anno.py &&\
    wget https://raw.githubusercontent.com/cumc/xqtl-pipeline/main/code/Dapars2_Multi_Sample.py &&\
    wget https://github.com/cumc/xqtl-pipeline/blob/main/code/gtf2bed12.py 
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh