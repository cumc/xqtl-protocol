FROM gaow/base-notebook

# Maintainer and author
LABEL maintainer="Amanda Tsai<at3535@cumc.columbia.edu>"

WORKDIR /tmp
USER root

# Install R packages
RUN R -e "remotes::install_github('stephenslab/susieR',build_vignettes=FALSE)"
RUN R -e "install.packages('Ckmeans.1d.dp', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('KRIS', dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN curl -so polyfun.yml https://raw.githubusercontent.com/omerwe/polyfun/master/polyfun.yml

# Create the environment:
RUN conda env create -f polyfun.yml
# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "polyfun", "/bin/bash", "-c"]
# The code to run when container is started:
ENTRYPOINT ["conda", "run", "-n", "polyfun", "/bin/bash"]