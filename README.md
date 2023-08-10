# ADSP FunGen-xQTL computational protocol

Developed for reproducible molecular QTL analysis for the NIH/NIA Alzheimer's Disease Sequencing Project Functional Genomics Consortium.

![QTL Diagram](code/images/complete_workflow.png)


## How to use this resource

### Standardized reference data

Reference data are standardized and curated by the ADSP FGC Standardization Workgroup in coordination with [NIAGCADS](https://www.nia.nih.gov/research/ad-genetics). Please find reference data specifications on [ADSP Dashboard](https://www.niagads.org/adsp/content/adspgcadgenomeresources-v2pdf).

### Software environment

We have prepared containerized software environment through both Docker and Singularity virtualization systems to facilicate software environment setup and to aid in software reproducibility. For those not familiar with this concept please [check out this wiki page of virtualization](https://en.wikipedia.org/wiki/OS-level_virtualization) and [an explanation on Docker website](https://www.docker.com/resources/what-container/).

### Pipeline execution

Pipelines in this repository are written in the [Script of Scripts (SoS) workflow language](https://vatlab.github.io/sos-docs/). Like most other workflow languages, SoS workflows can **distribute and execute computing jobs directly in High Performance Computing cluster**. It can also use **containers (Docker or Singularity)** to help with setting up computational environment and improve reproducibility. Unlike most other workflow languages, SoS workflows are created using SoS Notebooks (based on Ipython Notebook and developed in [Jupyter](https://jupyter.org/)) which allow for both **scientific narrative and pipeline scripts in the same document**. Unlike typical Jupyter Notebooks intended for interactive data analysis, SoS workflows written in Jupyter Notebooks can be executed directly as command line scripts either on a local computer or in a HPC environment. 

We provide this [toy example for running SoS pipeline on a typical HPC cluster environment](https://github.com/cumc/xqtl-pipeline/blob/main/code/misc/Job_Example.ipynb). First time users are encouraged to try it out in order to help setting up the computational environment necessary to run the analysis in this protocol.

### Source code

- Source code of pipelines and containers implemented in this repository are available at https://github.com/cumc/xqtl-pipeline/tree/main/code. 
- Container specifications for required software environments are available at https://github.com/cumc/xqtl-pipeline/tree/main/container.

### Organization of the resource

The website https://cumc.github.io/xqtl-pipeline is generated from files under the `code` folder of the source code repository. The `pipeline` folder contains symbolic links automatically generated for pipeline files under `code.` Users are encouraged to clone this repository, and execute from the root of the repository folders by typing 

```
sos run pipelines/<pipeline_file>.ipynb
```

that is, executing the symbolic links directly to perform the analysis. 

The logic of the entire xQTL analysis workflow is roughly reflected on the **left sidebar**:

- The **GETTING STARTED**  section serves as the main landing page or index of the xQTL protocol, guiding users through the various pipelines implemented in this repository. It's structured to mirror the logic of the xQTL analysis we've crafted. Because this page provides pointers to other sections, users can primarily focus here without having to sift through the rest of the pages pages in this repository.
- The **COMMAND GENERATOR** section is designed as a one-stop hub for "push button" commands, enabling users to generate the full QTL analysis pipeline workflow scripts from a straightforward configuration file. Notebooks within this section are intended to be **executed as command line software** for data analysis command generation. Users can then execute these generated commands directly to conduct all preset analyses. Alternatively, they can tweak the commands to cater to particular analysis requirements. The configuration file serves a dual purpose: streamlining control and maintaining a record of executed workflows.
- Other sections in bold fonts provide an array of available analyses, presented roughly from upstream to downstream processes. Most of these sections feature ***mini-protocols***, represented as clickable, non-bold text under each analysis category, leading to specific notebooks. These notebooks detail the commands necessary for the analyses defined in the respective mini-protocols. Predominantly tutorial-based, they are designed to be **executed interactively in Jupyter or via the command terminal**, allowing users to navigate through the SoS pipelines step by step. A few of these sections serve as actual ***pipeline modules*** which we'll discuss next (see below).
- *Mini-protocols*, as mentioned earlier, can be expanded by clicking the downward arrows, revealing the SoS implementations of ***pipeline modules***. These represent the crux of the pipeline implementations and are intended to be **executed as command line software**. They're also **self-contained**, allowing for reusability beyond the specific context of xQTL data analysis.

### Getting started

- In order to run the xQTL protocol on your computer (or a High Performance Computing cluster), please install Script of Scripts [(see here for a tutorial to set it up with `micromamba`)](https://wanggroup.org/orientation/jupyter-setup). 
    - For Linux and Mac desktop users you can either install the container [`Singularity`](https://sylabs.io/singularity/) (Singularity **version >=3.9.4**) or [`Docker`](https://www.docker.com/). In the xQTL project we primarily use `Singularity`.
    - For Linux-based HPC users, your system may already have `Singularity` installed. If not please communicate with the IT support for the HPC. Typically Docker is not allowed on HPC.
    - For Windows users, you will need to install [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (we have tested it on WSL2 and not on WSL1) and then install `Singularity` within WSL as instructed in this [post](https://www.blopig.com/blog/2021/09/using-singularity-on-windows-with-wsl2/).
- We have published example data-sets and `Singularity` containers images to [this Synapse folder](https://www.synapse.org/#!Synapse:syn36416559/files/). The instruction for downloading the data programmatically can be found [here](https://help.synapse.org/docs/Upload-and-Download-Data-in-Bulk.2003796248.html). To setup a synapse client, please follow [this post](https://help.synapse.org/docs/Installing-Synapse-API-Clients.1985249668.html).
  - In the `test_data` folder, you can find the data, prefixed with **MWE** (Minimal Working Example), used to perform unit testing for each module (i.e., whether there is anything wrong within the code).
  - In the `protocol_data` folder, you can find a more sophisticated collection of data, which are used to demonstrate the complete usage of our protocole in [this notebook](https://cumc.github.io/xqtl-pipeline/code/xqtl_protocol_demo.html) with [source code](https://github.com/cumc/xqtl-pipeline/blob/main/code/xqtl_protocol_demo.ipynb).
  - In the `container` folder above, you can find the Singularity images released for the software environment. If you use Docker (eg on a Linux or Mac Desktop) you **do not** need to download this folder.
- Please clone this repository https://github.com/cumc/xqtl-pipeline to your computer. This is the source code for this resource. All pipelines are symbolic links under `pipeline` folder to various notebooks under `code` folder. You can follow our mini-protocols to run the pipelines under `pipeline` folder.

### See Also

- Some analysis from FunGen-xQTL project using our protocol can be found in the [`cumc/brain-xqtl-analysis` github repo](https://github.com/cumc/brain-xqtl-analysis)

## Our team

This repository is developed by the NIA FunGen-xQTL consortium.

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

FunGen leadership

- Philip De Jager, Department of Neurology, Columbia University
- Carlos Crunchaga, Department of Psychiatry, Neurology and Genetics, Washington University in St. Louis

FunGen-xQTL methods and data integration working group

- Gao Wang (work group leader), Department of Neurology, Columbia University
- Xiaoling Zhang, Departments of Medicine and Biostatistics, Boston University
- Edoardo Marcora, Departments of Neuroscience, Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai
- Fanny Leung (also leads data standardization WG), Department of the Pathology and Laboratory Medicine, University of Pennsylvania
