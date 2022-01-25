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
RUN wget https://raw.githubusercontent.com/Xu-Dong/Exon_Intron_Extractor/main/scripts/gtf2bed12.py &&\
    wget https://raw.githubusercontent.com/seriousbamboo/DaPars2/master/src/DaPars_Extract_Anno.py &&\
    wget https://raw.githubusercontent.com/seriousbamboo/DaPars2/master/src/DaPars2_Multi_Sample_Multi_Chr.py &&\
    wget https://raw.githubusercontent.com/seriousbamboo/DaPars2/master/src/Dapars2_Multi_Sample.py
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh