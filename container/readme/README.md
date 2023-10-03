# Container provisioning pipeline
Containers in the xQTL pipeline are built from conda environments and pushed to container registry using GitHub Actions.  We use the [micromamba container](https://hub.docker.com/r/mambaorg/micromamba) as our base image and install a Conda environment into the container.  Conda packages are only used from [conda-forge](https://anaconda.org/conda-forge), [bioconda](https://anaconda.org/bioconda), and the [personal channel](https://anaconda.org/dnachun) of contributor @danielnachun.

Note: the minimum required version of Singularity for our containers is **3.6.0**.  This release made backwards incompatibile changes to the container signature format that mean any containers built with Singularity 3.6.0 or newer will not work on 3.5.3 or older.  We use the latest version of Apptainer (renamed from Singularity) to build and sign our containers. 

## Pipeline description
The GitHub Actions are located in the [.github/workflows](https://github.com/cumc/xqtl-pipeline/tree/main/.github/workflows) folder, while the Conda environments and the CSV table used to generate those environments are found in this folder.  The environments for a single container are each in their own folder e.g. `bioinfo/bioinfo.yml` – currently each container only has a single environment in it, but our pipeline allows for multiple environments in the same container should the need arise.  We build containers for both Docker and Singularity and push the containers to the registries at [ghcr.io](https://ghcr.io), [quay.io](https://quay.io), and [docker.io](https://docker.io).

### Summary of pipeline steps
1. When `containers.csv` is modified by a pull request, the `update_container_list.yml` [workflow](https://github.com/cumc/xqtl-pipeline/tree/main/.github/workflows/upload_container_list.yml) is triggered which calls a [Python script](https://github.com/cumc/xqtl-pipeline/tree/main/.github/workflows/export_envs.py) that regenerates the environment YML files.
2. Any changes to the YML files are added as a new commit in the pull request.
3. The new commit triggers the `container_build.yml` [workflow](https://github.com/cumc/xqtl-pipeline/tree/main/.github/workflows/container_build.yml), which will build Docker and Singularity containers for any YML files that have been changed.
4. When the pull request is merged, the `container_upload.yml` [workflow](https://github.com/cumc/xqtl-pipeline/tree/main/.github/workflows/container_upload.yml) is triggered.  This pipeline builds the containers again and pushes them to the container registries.

## How to update or build new containers
To update an existing container with new packages or new versions of packages or to build a new container, make a pull request in this repository to change `containers.csv`.  Do **not** directly modify the YML files - they are automatically generated and will be overriden by the next pull request.  To get started:
1. If you have not already done so, create a fork of this repository to your own GitHub account and clone it locally with `git clone git@github.com:YOUR_USERNAME/xqtl-pipeline`.  
2. Checkout a new branch: `git checkout -b my_container_update`.
3. Edit container.csv to make the changes you want to the environments using your preferred spreadsheet editor (Excel, Google Sheets, Numbers, LibreOffice, etc.).  Please refer to the documentation in the next section for details on what to add to the table.
4. Commit the changes you have made with `git commit -m "update my_container"`.
5. Push the branch to the upstream wtih `git push -f origin my_container_update`.
6. Once the checks all pass with no errors, ask @gaow or @danielnachun to merge the pull request.

## Layout of container list
The table requires several pieces of information for each package:
1. The **package name** exactly as it is written for the package in anaconda.org.  Conda package names are not case sensitive but please use all lowercase letters for readability.  Please note the naming conventions for conda packages - for some languages such as R and Perl, most package names begin with `r-`/`perl-`, which Python and C/C++, for example do not.  Also note that Bioconductor R packages begin with `bioconductor-` instead of `r-`.
2. The **URL** of the package home page or source code.  For packages in PyPI (Python), CRAN (R), Bioconductor (R), and CPAN (Perl), please provide use the link to the package repository (i.e. https://pypi.org/project/qtl/) instead of GitHub or a dedicate web page.
3. The main **language** the package is written in.  If you are not sure what to put for some packages, please ask in the pull request.
4. The **version** of the package.  Please use the latest version available in anaconda.org and not the newest version available upstream.  If the latest version of anaconda.org is too old, please open an issue or pull request at https://github.com/danielnachun/recipe_staging to get a newer version on the personal channel of @danielnachun.  If you need to use an older version than the latest version on anaconda.org, please provide an explanation in the pull request.
5. The **git revision** of the package.  This is only needed for packages built from GitHub that are not tagged, or have had many changes since their last tag.  This only applies to some packages hosted on the personal channel of @danielnachun – please ask for help if you not sure what to add here.
6. The **channel** the package is available at on anaconda.org.  The preferred order is conda-forge > bioconda > dnachun.
7. The **environment** the package is needed in.  Multiple environment can be separated by a comma.  Currently we are trying to use the same versions across all environments, but accomodations can be made if this is not possible.
8. The **license** for the package.  Please follow the formating used in the identified column at https://spdx.org/licenses/ – you will be asked to change the format if you do not.  If you are not sure what the license should be for a package, please ask for help with this in the pull request.
9. The URL of the **conda package**.  It should always be of the form https:// anaconda.org / [CHANNEL] / [PACKAGE] e.g. https://anaconda.org/dnachun/hello.

## Troubleshooting
### Packages not found
If you see the error `nothing provides requested PACKAGE_NAME` you have specified a package does not exist.  This could be because the package name is spelled incorrectly, or because there is no available conda package in `conda-forge`, `bioconda` or the personal channel of @danielnachun.  Please open an issue or pull request at https://github.com/danielnachun/recipe_staging to request a new package to be built.
### Dependency conflict
If you see the error `libmamba Could not solve for environment specs` and there are one or more errors that say `The following packages are incompatible`, there is a dependency conflict.  You can test building the environment on an x86-64 Linux system (either a real machine or a virtual machine) with `mamba env create -f ENVIRONMENT.yml`.  In some cases, the dependency conflict cannot be solved by simply changing versions - one or more packages may need to be rebuilt.  Please ask for help if you get stuck!
### Transient network errors
On occasion, you may see errors that a download or upload has failed while building or pushing a container.  There are many different kinds of error messages that can occur with this.  To resolve this issue, click on the "Re-run jobs" button in the upper right corner, and choose "Re-run failed jobs".  Usually the remaining jobs will work after the second try.

## List of packages

