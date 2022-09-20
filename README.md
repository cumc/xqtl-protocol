# xQTL analysis pipeline

Developed for reproducible QTL analysis for the NIH/NIA Alzheimer's Disease Sequencing Project Functional Genomics Consortium.

## How to use this resource

### Standardized reference data

Reference data are standardized and curated by the ADSP FGC Standardization Workgroup in coordination with [NIAGCADS](https://www.nia.nih.gov/research/ad-genetics). Please find reference data specifications on [ADSP Dashboard](https://www.niagads.org/adsp/content/adspgcadgenomeresources-v2pdf).

### Pipeline execution

Pipelines in this repository are written in the [Script of Scripts (SoS) workflow language](https://vatlab.github.io/sos-docs/). Like most other workflow languages, SoS workflows can **distribute and execute computing jobs directly in High Performance Computing cluster**. It can also use **containers (Docker or Singularity)** to help with setting up computational environment and improve reproducibility. Unlike most other workflow languages, SoS workflows are created using SoS Notebooks (based on Ipython Notebook and developed in [Jupyter](https://jupyter.org/)) which allow for both **scientific narrative and pipeline scripts in the same document**. Unlike typical Jupyter Notebooks intended for interactive data analysis, SoS workflows written in Jupyter Notebooks can be executed directly as command line scripts either on a local computer or in a HPC environment. 

We provide this [toy example for running SoS pipeline on a typical HPC cluster environment](https://github.com/cumc/xqtl-pipeline/blob/main/code/misc/Job_Example.ipynb). First time users are encouraged to try it out in order to help setting up the computational environment necessary to run the QTL analysis.

### Source code

- Source code of pipelines and containers implemented in this repository are available at https://github.com/cumc/xqtl-pipeline/tree/main/code. 
- Container configurations are available at https://github.com/cumc/xqtl-pipeline/tree/main/container.

### Getting started

- Working examples and containers are available through a request to access [this Synapse folder](https://www.synapse.org/#!Synapse:syn36416559/files/).
  - In the `test_data` folder, you can find the data, prefixed with **MWE**,  used to perform unit testing for each module (i.e., whether there is anything wrong within the code).
  - In the `protocol_data` folder, you can find a more sophisticated collection of data, which are used to demonstrate the usage of xqtl-protocol in this [WIP notebook that outline the comprehensive xqtl-analysis procedure](https://github.com/cumc/xqtl-pipeline/blob/main/code/xqtl_protocol_demo.ipynb).
- Under the `container` folder above, you can find the `singularity` image release for the software environment. You can also build the singularity image from configuration files at: https://github.com/cumc/xqtl-pipeline/tree/main/container/singularity.


### Organization of the resource

The website https://cumc.github.io/xqtl-pipeline is generated from files under the `code` folder of the source code repository. The `pipeline` folder contains symbolic links automatically generated for pipeline files under `code.` The logic of the entire xQTL analysis workflow is roughly reflected on the **left sidebar**:

- The **COMMAND GENERATOR** section is reserved for "push button" commands that generate the entire QTL analysis pipeline workflow script from a simple configuration file. Notebooks under these sections are meant to be **executed as command line software** to generate data analysis commands. The generated commands can be executed as is to complete all available analyses or can be used to help customize specific analysis tasks by making modifications to them. The configuration file itself helps centralized control and bookkeeping of workflows executed.
- Other sections in bold contain various types of analysis available, roughly showing in order from upstream to downstream analysis. We will refer to them as ***analysis groups***, which are further divided into ***protocols*** by various non-bold, clickable text under each analysis group linking to some notebooks. These notebooks illustrate commands to perform analysis implemented in the protocol. Most of them are "tutorials" in nature and are meant to be **executed interactively in Jupyter or in the command terminal** to run the SoS pipelines line by line. A few are the actual ***pipeline modules*** implementing pipelines in SoS, as will be discussed next.
- *Protocols* can be expanded by clicking on the down arrows to access the SoS workflows implementation of ***pipeline modules***. These are the core pipeline implementations to be **executed as command line software**and are meant to be **self-contained** --- they may be used in other contexts not specific to the xQTL data analysis. Each of these pipeline

## xQTL workflow schema (WIP)

To perform a complete analysis from molecular phenotype calling up to the xqtl discovery as demonstrated in the xqtl-protocol paper(in preparation), please read the ***mini protocols*** in the following orders:
### Molecular phenotype calling
1. [Reference data procession](https://cumc.github.io/xqtl-pipeline/code/data_preprocessing/reference_data.html)
2. [Gene expression calling for eQTL](https://cumc.github.io/xqtl-pipeline/code/molecular_phenotypes/bulk_expression.html)
3. [Alternative splicing calling for sQTL](https://cumc.github.io/xqtl-pipeline/code/molecular_phenotypes/splicing.html)
### xQTL discovery
1. [Genotype Processing](https://cumc.github.io/xqtl-pipeline/code/data_preprocessing/genotype_preprocessing.html)
2. [Covariates Processing](https://cumc.github.io/xqtl-pipeline/code/data_preprocessing/covariate_preprocessing.html)
3. [Phenotype Processing](https://cumc.github.io/xqtl-pipeline/code/data_preprocessing/phenotype_preprocessing.html)
4. [CisQTL Association Scan](https://cumc.github.io/xqtl-pipeline/code/association_scan/cisQTL_scan.html)
5. [TransQTL Association Scan](https://cumc.github.io/xqtl-pipeline/code/association_scan/transQTL_scan.html)

***Noted that mini protocols 1,2,3 of xQTL discovery are interwined with each other. Please read the mini protocol for the transition point from one to another.***

![QTL Diagram](code/images/complete_workflow.png)


## See Also
- Some example analysis using our pipeline can be found in the [brain-xqtl-analysis github repo](https://github.com/cumc/brain-xqtl-analysis)

## Contributors

This repository is developed by the ADSP FG Brain xQTL consortium.

### Developers

Lead developers

- Hao Sun, Department of Neurology, Columbia University
- Gao Wang, Department of Neurology, Columbia University

Contributors

- Wenhao Gou, Department of Biostatistics, Columbia University
- Liucheng Shi, Department of Biostatistics, Columbia University
- Xuanhe Chen, Department of Biostatistics, Columbia University
- Amanda Tsai, Department of Biostatistics, Columbia University  

Brain xQTL project leadership

- Philip De Jager, Department of Neurology, Columbia University
- Carlos Crunchaga, Department of Psychiatry, Neurology and Genetics, Washington University in St. Louis

Brain xQTL methods and data integration work group

- Gao Wang (work group leader), Department of Neurology, Columbia University
- Xiaoling Zhang, Departments of Medicine and Biostatistics, Boston University
- Edoardo Marcora, Departments of Neuroscience, Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai

