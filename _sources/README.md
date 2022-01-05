# xQTL analysis pipeline

Developed for reproducible QTL analysis for the NIH/NIA Alzheimer's Disease Sequencing Project Functional Genomics Consortium

https://cumc.github.io/xqtl-pipeline

## How to use this resource

The left sidebar lists the analysis to be performed for QTL analysis, roughly in sequential orders.

- The **COMPLETE PIPELINES** section is reserved for publishing the code used for various QTL analysis. Contents in this section can be generated automatically from workflow notebooks in the other sections, as will be discussed next.
- The rest sections in bold are various types of analysis available, from generating the molecular phenotypes to performing some of the selected multi-omics data integration analysis.
- Text under each bold section title shows the complete workflow commands to perform analysis implemented using the Script of Scripts (SoS) language, as will be discussed next.
- The workflows can be expanded by clicking on the down arrows to check out the SoS workflows implementing each task as an analysis module. These are the core pipeline implementations. Each of these pipeline modules are documented with some background information, required input, expected output, and a minimal working example, followed by the actual code implementation.

## Overall xQTL workflow schema

![QTL Diagram](images/complete_workflow.png)

## Genotype, phenotype and covariates data preprocessing workflow schema

![Data preprocessing diagram](images/data_preprocessing.png)

## Contributors

This repository is developed by the ADSP FG Brain xQTL consortium.

### Lead developers

Project leader

- Philip De Jager, Department of Neurology, Columbia University

Brain xQTL methods and analysis subgroup

- Gao Wang, Department of Neurology, Columbia University
- Hao Sun, Department of Neurology, Columbia University
- Wenhao Gou, Department of Biostatistics, Columbia University
- Amanda Tsai, Department of Biostatistics, Columbia University  
- Xiaoling Zhang, Departments of Medicine and Biostatistics, Boston University

### Consortium collaborators

Single-cell nucleotide RNA-seq

- 

- 
