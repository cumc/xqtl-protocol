FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
USER root
RUN cd /tmp
## Install Git
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends git-all
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin
RUN pip install tensorqtl fastparquet --no-cache-dir
## Installing the multipy package for the qvalue. Noted: the pip version dont works with python3. Only the github version works.
RUN pip install git+https://github.com/puolival/multipy.git
## Installing R packages for the qvalue
RUN R --slave -e "install.packages(c('BiocManager','tidyr'))"
RUN R --slave -e "BiocManager::install('qvalue')"
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh
