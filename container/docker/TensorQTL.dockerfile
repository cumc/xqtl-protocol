FROM gaow/base-notebook
LABEL MAINTAINER Hao Sun <hs3163@cumc.columbia.edu>
RUN cd /tmp
su -  root # USER root
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200616.zip && \
unzip plink_linux_x86_64_20200616.zip && mv plink /usr/local/bin
RUN pip install tensorqtl
RUN pip install fastparquet
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh