FROM gaow/base-notebook
LABEL maintainer="Wenhao Gou<wg2364@cumc.columbia.edu>"
su -  root # USER root
RUN pip install qtl
RUN wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2\
&& tar xvf bcftools-1.2.tar.bz2 \
&& cd bcftools-1.2 \
&& make\
&& make install\
&& cd htslib-1.2.1\
&& make\
&& make install
RUN rm -rf bcftools-1.2 \
&& rm -rf bcftools-1.2.tar.bz2
CMD exec /bin/bash "$@"