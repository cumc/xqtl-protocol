FROM gaow/base-notebook AS spython-base
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
su -  root # USER root
RUN apt-get update \
&& apt install -y --no-install-recommends  git-all  libboost-all-dev
RUN git clone https://github.com/corbinq/apex.git \
&& cd apex/bin \
&& gunzip apex_Linux_x86_64.gz \
&& mv apex_Linux_x86_64 /usr/local/bin/apex && chmod +x /usr/local/bin/apex
RUN echo "cd /tmp" >> /entrypoint.sh
RUN echo "exec /bin/bash "$@"" >> /entrypoint.sh
RUN chmod u+x /entrypoint.sh
CMD /bin/bash /entrypoint.sh