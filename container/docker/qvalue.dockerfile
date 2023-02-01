FROM gaow/base-otebook
LABEL maintainer="Hao Sun<hs3163@cumc.columbia.edu>"
USER root


# Install R pkg and qvalue
RUN R --slave -e "install.packages(c('tidyverse', 'BiocManager'))"
RUN R --slave -e "BiocManager::install('qvalue')"

CMD exec /bin/bash "$@"
