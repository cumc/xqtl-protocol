# ADSP-FG xQTL protocol

Developed for reproducible molecular QTL analysis for the NIH/NIA Alzheimer's Disease Sequencing Project Functional Genomics Consortium.

![QTL Diagram](code/images/complete_workflow.png)


## How to use this resource

### Standardized reference data

Reference data are standardized and curated by the ADSP FGC Standardization Workgroup in coordination with [NIAGCADS](https://www.nia.nih.gov/research/ad-genetics). Please find reference data specifications on [ADSP Dashboard](https://www.niagads.org/adsp/content/adspgcadgenomeresources-v2pdf).

### Software environment

We have prepared containerized software environment through both Docker and Singularity virtualization systems to facilicate software environment setup and to aid in software reproducibility. For those not familiar with this concept please [check out this wiki page of virtualization](https://en.wikipedia.org/wiki/OS-level_virtualization) and [an explanation on Docker website](https://www.docker.com/resources/what-container/).

### Pipeline execution

Pipelines in this repository are written in the [Script of Scripts (SoS) workflow language](https://vatlab.github.io/sos-docs/). Like most other workflow languages, SoS workflows can **distribute and execute computing jobs directly in High Performance Computing cluster**. It can also use **containers (Docker or Singularity)** to help with setting up computational environment and improve reproducibility. Unlike most other workflow languages, SoS workflows are created using SoS Notebooks (based on Ipython Notebook and developed in [Jupyter](https://jupyter.org/)) which allow for both **scientific narrative and pipeline scripts in the same document**. Unlike typical Jupyter Notebooks intended for interactive data analysis, SoS workflows written in Jupyter Notebooks can be executed directly as command line scripts either on a local computer or in a HPC environment. 

We provide this [toy example for running SoS pipeline on a typical HPC cluster environment](https://github.com/cumc/xqtl-pipeline/blob/main/code/misc/Job_Example.ipynb). First time users are encouraged to try it out in order to help setting up the computational environment necessary to run the QTL analysis.

### Source code

- Source code of pipelines and containers implemented in this repository are available at https://github.com/cumc/xqtl-pipeline/tree/main/code. 
- Container configurations (both Docker and Singularity) for required software environments are available at https://github.com/cumc/xqtl-pipeline/tree/main/container.

### Organization of the resource

The website https://cumc.github.io/xqtl-pipeline is generated from files under the `code` folder of the source code repository. The `pipeline` folder contains symbolic links automatically generated for pipeline files under `code.` The logic of the entire xQTL analysis workflow is roughly reflected on the **left sidebar**:

- The **COMMAND GENERATOR** section is reserved for "push button" commands that generate the entire QTL analysis pipeline workflow script from a simple configuration file. Notebooks under these sections are meant to be **executed as command line software** to generate data analysis commands. The generated commands can be executed as is to complete all available analyses or can be used to help customize specific analysis tasks by making modifications to them. The configuration file itself helps centralized control and bookkeeping of workflows executed.
- Other sections in bold contain various types of analysis available, roughly showing in order from upstream to downstream analysis. They are consisted of ***mini-protocols*** as various non-bold, clickable text under each analysis group linking to some notebooks. These notebooks illustrate commands to perform analysis implemented in each mini-protocol. Most of them are "tutorials" in nature and are meant to be **executed interactively in Jupyter or in the command terminal** to run the SoS pipelines line by line. A few are the actual ***pipeline modules*** implementing pipelines in SoS, as will be discussed next.
- *Mini-protocols* can be expanded by clicking on the down arrows to access the SoS workflows implementation of ***pipeline modules***. These are the core pipeline implementations to be **executed as command line software**and are meant to be **self-contained** --- they may be used in other contexts not specific to the xQTL data analysis.

### Getting started

- In order to run the protocol on your computer (or a High Performance Computing cluster), please install Script of Scripts [(see here for a tutorial to set it up)](https://wanggroup.org/orientation/jupyter-setup). 
    - For Linux desktop users you can either install the container [Singularity](https://sylabs.io/singularity/) or [Docker](https://www.docker.com/).
    - For Linux-based HPC users, your system may already have Singularity installed. If not please communicate with the IT support for the HPC. Typically Docker is not allowed on HPC.
    - For Mac desktop users, it is best to install and use [Docker](https://www.docker.com/).
    - For Windows users, you will need to install [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (we have tested it on WSL2 and not on WSL1) and then install singularity within it as instructed in this [post](https://www.blopig.com/blog/2021/09/using-singularity-on-windows-with-wsl2/).
- We have published example data-sets and singularity containers images to [this Synapse folder](https://www.synapse.org/#!Synapse:syn36416559/files/). The instruction for downloading the data programmatically can be found [here](https://help.synapse.org/docs/Upload-and-Download-Data-in-Bulk.2003796248.html). To setup a synapse client, please follow [this post](https://help.synapse.org/docs/Installing-Synapse-API-Clients.1985249668.html).
  - In the `test_data` folder, you can find the data, prefixed with **MWE**,  used to perform unit testing for each module (i.e., whether there is anything wrong within the code).
  - In the `protocol_data` folder, you can find a more sophisticated collection of data, which are used to demonstrate the complete usage of our protocole in [this notebook](code/xqtl_protocol_demo.html) with [source code](https://github.com/cumc/xqtl-pipeline/blob/main/code/xqtl_protocol_demo.ipynb).
  - In the `container` folder above, you can find the Singularity images released for the software environment. If you use Docker (eg on a Linux or Mac Desktop) you **do not** need to download this folder.

### See Also

- Some analysis from ADSP-FG xQTL project using our protocol can be found in the [`cumc/brain-xqtl-analysis` github repo](https://github.com/cumc/brain-xqtl-analysis)

## Our team

This repository is developed by the ADSP FG Brain xQTL consortium.

### Developers

Lead developers

- Hao Sun, Department of Neurology, Columbia University
- Gao Wang, Department of Neurology, Columbia University

Contributors

- Xuanhe Chen, Department of Biostatistics, Columbia University
- Wenhao Gou, Department of Biostatistics, Columbia University
- Yuqi Miao, Department of Biostatistics, Columbia University
- Liucheng Shi, Department of Biostatistics, Columbia University
- Amanda Tsai, Department of Biostatistics, Columbia University  

### Leadership

Brain xQTL project leadership

- Philip De Jager, Department of Neurology, Columbia University
- Carlos Crunchaga, Department of Psychiatry, Neurology and Genetics, Washington University in St. Louis

Brain xQTL methods and data integration work group

- Gao Wang (work group leader), Department of Neurology, Columbia University
- Xiaoling Zhang, Departments of Medicine and Biostatistics, Boston University
- Edoardo Marcora, Departments of Neuroscience, Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai
- Fanny Leung (also leads data standardization WG), Department of the Pathology and Laboratory Medicine, University of Pennsylvania
