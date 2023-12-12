## Title



## Authors




## Abstract




## Introduction




### Development of the protocol




### Applications of the method




### Comparison with other methods




### Experimental Design




#### Reference data (Step 1)
##### A.  Reference Data

Please refer to the protocol website for more information on this miniprotocol.
#### Molecular Phenotypes (Step 2)
##### A.  RNA-seq expression


Our pipeline follows the GTEx pipeline from the [GTEx](https://gtexportal.org/home/aboutGTEx#staticTextAnalysisMethods)/[TOPMed](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) project for the quantification of expression from RNA-seq data. Either paired end or single end fastq.gz files may be used as the initial input. Different read [strandedness options](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) are possible, including rf, fr, unstranded or strand_missing, based on library types from Signal et al [[cf. Signal et al (2022)](https://doi.org/10.1186/s12859-022-04572-7)]. Strand detection steps are included in this pipeline to automaticaly detect the strand of the input fastq file by levradging the gene count table ouptut from [STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf). Read length data is required and is set to a default value of 100 for reads of zero length. 



We recommend using fastqc to quality control reads, followed by the removal of adaptors with fastp if necessary. After quality control has been conducted, STAR may be used to align reads to the reference genome. We combine this mapping step with an additionaly quality control step using Picard to mark duplicate reads and collect other metrics. 



Gene-level RNA expression is called with RNA-SeQC v2.4.2 using a reference gtf that has been collapsed to contain genes instead of gene transcripts [[cf. DeLuca et al., Bioinformatics, 2012](https://doi.org/10.1093/bioinformatics/bts196)]. Reads are only included if they are uniquely mapped (mapping quality of 255 for STAR BAMs), if they are aligned in proper pairs, if the read alignment distance was less than or equal to six (i.e. alignments must not contain more than six non-reference bases), and if they are fully contained within exon boundaries. Reads overlapping introns are not included. Exon level read counts are also called by RNA-SeQC. If a read overlaps multiple exons, then a fractional value equal to the portion of the read contained within the exon is allotted. We call transcript-level RNA expression using RSEM v1.3.0 with a transcript reference gtf. 

Quality control of the TPM matrices follows methods outlined by [GTEx V8](https://gtexportal.org/home/aboutGTEx#staticTextAnalysisMethods). First, genes are removed from the matrices if over 20% of samples have have a TPM expression level of 10% or less. Sample level filtering includes three checks to detect sample outliers. Samples are only removed if they are marked as outliers in all three checks. 



The first looks at Relative Log Expression (RLE). It is assumed that most gene expression values in a sample should be near a mean value and only a few genes are differentially expressed. To calculate RLE, for each gene *i*, calculate its median in *N* samples as *Medi*. Then for each sample *j* and expression value *eij*, count the difference between *eij* and *Medi* (*deij = eij-Medi*). Then create boxplots for each sample based on *deij* and sort by interquartile range. Samples with larger interquartile ranges are more likely to be outliers. 



The second quality control check looks at heirarchical clustering of samples. Samples are expected to have short distances between others and therefore should cluster homogeneously, so distant samples are expected to be outliers. Distance is calculated as 1-spearman correlation for heirarchical clustering. The top 100 genes sorted by variance are used to calculate Mahalanobis distance. Chi2 p-values are then calculated based on Mahalanobis distance. Clusters with 60% or more samples with Bonferroni corrected p-values less than 0.05 are marked as outliers. 



The third and last quality control check looks at D-statistics which represent the average correlation between a sample's expression and other sampels. Samples with low D-statistics are likely to be outliers. 









The normalization step follows steps used by the GTeX pipeline. Genes are first filtered to keep genes where TPM is greater than 10% in at least 20% of the samples. They are also kept if read counts is greater than 6 in at least 20% of the samples. The filtered data is then normalized using the Trimmed Mean of M-value (TMM) method. 



##### B.  Alternative splicing from RNA-seq data


Our pipeline calls alternative splicing events from RNA-seq data using leafcutter and psichomics to call the RNA-seq data from `fastq.gz` data which has been mapped to a reference genome using STAR with the wasp option. It implements the GTEx pipeline for GTEx/TOPMed project. Please refer to [this page](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) for detail. The choice of pipeline modules in this project is supported by internal (unpublished) benchmarks from GTEx group.



We use two different tools to quantify the many types of splicing events which are outlined in [[cf. Wang et al (2008)](https://doi.org/10.1038/nature07509)] and [[cf. Park et al (2018)](https://doi.org/10.1016/j.ajhg.2017.11.002)]. The first, leafcutter, quantifies the usage of alternatively excised introns. This collectively captures skipped exons, 5’ and 3’ alternative splice site usage and other complex events [[cf. Li et al 2018](https://doi.org/10.1038/s41588-017-0004-9)]. This method was previously applied to ROSMAP data as part of the Brain xQTL version 2.0.  The second, psichomics, quantifies specific splicing events [[cf. Agostinho et al. 2019](https://doi.org/10.1093/nar/gky888)].


Quality control and normalization are performed on output from the leafcutter and psichomics tools. The raw output data is first converted to bed format. Quality control involves the removal of features with high missingness across samples (default 40%) and the replacement of NA values in the remaining samples with mean existing values. Then introns with less than a minimal variation (default of 0.005) are removed from the data. Quantile-Quantile normalization is performed on the quality controlled data. 
#### Data Pre-processing (Step 3)
##### A.  Genotype data preprocessing


The goal of this module is to perform QC on VCF files, including 



1. Handling the formatting of multi-allelic sites. 

2. Genotype and variant level filtering based on genotype calling qualities. 

3. Known/novel variants annotation.

4. Summary statistics before and after QC, in particular the ts/tv ratio, to assess the effectiveness of QC.



1 and 2 will change the genotype data. 3 and 4 above are for explorative analysis on the overall quality assessment of genotype data in the VCF files. We annotate known and novel variants because ts/tv are expected to be different between known and novel variants, and is important QC metric to assess the effectiveness of our QC.



### Multi-allelic sites



Mult-allelic sites can be problematic in many ways for downstreams analysis, even of they are handled in terms of formatting after QC. We provide an optional workflow module to keep only bi-allelic sites from data, although by default we will include these sites in the VCF file we generate.

This notebook includes workflow for



- Compute kinship matrix in sample and estimate related individuals

- Genotype and sample QC: by MAF, missing data and HWE

- LD pruning for follow up PCA analysis on genotype, as needed



A potential limitation is that the workflow requires all samples and chromosomes to be merged as one single file, in order to perform both sample and variant level QC. However, in our experience using this pipeline with 200K exomes with 15 million variants, this pipeline works on the single merged PLINK file.

Steps to generate a PCA include 



- removing related individuals

- pruning variants in linkage disequilibrium (LD)

- perform PCA analysis on genotype of unrelated individuals

- excluding outlier samples in the PCA space for individuals of homogeneous self-reported ancestry. These outliers may suggest poor genotyping quality or distant relatedness.



### Limitations



1. Some of the PCs may capture LD structure rather than population structure (decrease in power to detect associations in these regions of high LD)

2. When projecting a new study dataset to the PCA space computed from a reference dataset: projected PCs are shrunk toward 0 in the new dataset

3. PC scores may capture outliers that are due to family structure, population structure or other reasons; it might be beneficial to detect and remove these individuals to maximize the population structure captured by PCA (in the case of removing a few outliers) or to restrict analyses to genetically homogeneous samples

For each chromosome, we compute a GRM using data excluding this chromosome. Computation is implemented using `GCTA` software package.

The module streamlines conversion between PLINK and VCF formats, specifically:



1. Conversion between VCF and PLINK formats

2. Split data (by specified input, by chromosomes, by genes)

3. Merge data (by specified input, by chromosomes)
##### B.  Phenotype data preprocessing


This pipeline is based on [`pyqtl`, as demonstrated here](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/eqtl_prepare_expression.py).



### Alternative implementation



Previously we use `biomaRt` package in R instead of code from `pyqtl`. The core function calls are:



```r

    ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "$[ensembl_version]")

    ensembl_df <- getBM(attributes=c("ensembl_gene_id","chromosome_name", "start_position", "end_position"),mart=ensembl)

```



We require ENSEMBL version to be specified explicitly in this pipeline. As of 2021 for the Brain xQTL project, we use ENSEMBL version 103.
##### C.  Covariate Data Preprocessing


This workflow implements 3 procedures for hidden factor analysis from omcis data:



1. The [Probabilistic Estimation of Expression Residuals (PEER) method](https://github.com/PMBio/peer/wiki/Tutorial), a method also used for GTEx eQTL data analysis. 

2. Factor analysis using Bi-Cross validation, Owen, Art & Wang, Jingshu. (2015). Bi-Cross-Validation for Factor Analysis. Statistical Science. 31. 10.1214/15-STS539. with software package `APEX` (Corbin Quick, Li Guan, Zilin Li, Xihao Li, Rounak Dey, Yaowu Liu, Laura Scott, Xihong Lin, bioRxiv 2020.12.18.423490; doi: https://doi.org/10.1101/2020.12.18.423490)

3. PCA with automatic determination of the number of factors to use. This is mainly inspired by a [recent benchmark from Jessica Li's group](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02761-4).





Overall, we will pick PCA based approach for the xQTL project, although additional considerations should be taken for single-cell eQTL analysis as investigated in [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02873-5).

### Expertise needed to implement the protocol




### Limitations




## Materials



### Software



### Hardware



## Procedure




> CRITICAL  
  To improve readability, the code outlined here highlights the notebook to run for each step and not the necessary parameters in some cases. Please refer to the protocol website for information on what parameters to include.
### 1. Reference data
#### A.  Reference Data

Please refer to the protocol website for more information on this miniprotocol.
### 2. Molecular Phenotypes
#### A.  RNA-seq expression

##### i. Perform data quality summary via `fastqc`
Timing <4 min

```

!sos run RNA_calling.ipynb fastqc \
    --cwd ../../output_test \
    --samples ../../PCC_sample_list_subset \
    --data-dir /restricted/projectnb/amp-ad/ROSMAP_PCC_AC/PCC/ \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    -c ../csg.yml  -q neurology
```


##### ii. Cut adaptor (Optional)
Timing ~10 min

```
!sos run RNA_calling.ipynb fastp_trim_adaptor \
    --cwd ../../output_test \
    --samples ../../PCC_sample_list_subset \
    --data-dir /restricted/projectnb/amp-ad/ROSMAP_PCC_AC/PCC/ \
    --STAR-index ../../reference_data/STAR_Index/ \
    --gtf ../../reference_data/reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.ERCC.gtf \
    --reference-fasta ../../reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta \
    --ref-flat ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformated.ERCC.gtf.ref.flat \
    --container ../../rna_quantification.sif \
    -c ../csg.yml  -q neurology
```


##### iii. Read alignment via STAR and QC via Picard
Timing <2 hours

```
!sos run RNA_calling.ipynb STAR_align \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \

```



```
!sos run RNA_calling.ipynb STAR_align \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \

```


##### iv. Call gene-level RNA expression via rnaseqc
Timing <30 min

```
!sos run RNA_calling.ipynb rnaseqc_call \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \

```


##### v. Call transcript level RNA expression via RSEM
Timing <X hours


```
!sos run RNA_calling.ipynb rsem_call \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \

```


##### vi. Multi-sample RNA-seq QC
Timing <15min

```
!sos run bulk_expression_QC.ipynb qc \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \

```


##### vii. Multi-sample read count normalization
Timing <10min

```
!sos run bulk_expression_normalization.ipynb normalize \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \

```


#### B.  Alternative splicing from RNA-seq data

##### i. Splicing Quantification

Timing <30min

```
!sos run splicing_calling.ipynb leafcutter \
    --container oras://ghcr.io/cumc/leafcutter_apptainer:latest \

```


Timing ~30min

```
!sos run splicing_calling.ipynb psichomics \
    --container oras://ghcr.io/cumc/psichomics_apptainer:latest \

```


##### ii. Splicing QC and Normalization
Timing  ~30min

```
!sos run splicing_normalization.ipynb leafcutter_norm \
    --container oras://ghcr.io/cumc/leafcutter_apptainer:latest \

```


Timing ~20min

```
!sos run splicing_normalization.ipynb psichomics_norm \
    --container oras://ghcr.io/cumc/psichomics_apptainer:latest \

```


### 3. Data Pre-processing
#### A.  Genotype data preprocessing


```
sos run VCF_QC.ipynb rename_chrs \
    --cwd reference_data --container bioinfo.sif
```



```
sos run VCF_QC.ipynb dbsnp_annotate \
    --cwd reference_data --container bioinfo.sif
```



```
sos run VCF_QC.ipynb qc    \
    --cwd MWE/output/genotype_1 --container bioinfo.sif -J 1 -c csg.yml -q csg
```



```
sos run VCF_QC.ipynb qc    \
    --cwd MWE/output/genotype_4 --container bioinfo.sif --add-chr
```



```

```



```

```



```

```


##### Perform QC on both rare and common variants

```
sos run xqtl-pipeline/pipeline/GWAS_QC.ipynb qc_no_prune \
   --container /mnt/vast/hpc/csg/containers/bioinfo.sif \

```


##### Sample match with genotype
Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb genotype_phenotype_sample_overlap \
        --container containers/bioinfo.sif \

```


##### Kinship QC

Timing <2 min

```
sos run pipeline/GWAS_QC.ipynb king \
    --container containers/bioinfo.sif \

```


##### Prepare unrelated individuals data for PCA

Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb qc \
   --container containers/bioinfo.sif \

```


Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb qc \
   --container containers/bioinfo.sif \

```



```
sos run GWAS_QC.ipynb qc_no_prune \
    --container container/bioinfo.sif
```


##### Estimate kinship in the sample


```
sos run GWAS_QC.ipynb king \
    --container container/bioinfo.sif
```


##### Sample selection and QC the genotype data for PCA


```
sos run GWAS_QC.ipynb qc \
    --container container/bioinfo.sif
```


##### PCA analysis for unrelated samples
Timing <2 min

```
sos run pipeline/PCA.ipynb flashpca \
   --container containers/flashpcaR.sif \

```



```

```


##### Projection of related individuals


```
sos run PCA.ipynb project_samples \
  --container container/flashpcaR.sif

```



```

```


##### Finalize genotype QC by PCA for homogenous population


```
sos run GWAS_QC.ipynb qc_no_prune \
    --container container/bioinfo.sif
```



```
sos run GWAS_QC.ipynb qc_no_prune \
    --container container/bioinfo.sif
```



```
sos run genotype_formatting.ipynb merge_plink \
    --container container/bioinfo.sif
```


##### Split data by population


```

```


##### For each population, do variant level and sample level QC on unrelated individuals, in preparation for PCA analysis


```
    sos run GWAS_QC.ipynb qc \
        --container container/bioinfo.sif
```


##### For each population, extract previously selected variants from related individuals in preparation for PCA, only applying missingness filter at sample level, if applicable


```
    sos run GWAS_QC.ipynb qc_no_prune \
        --container container/bioinfo.sif -s force
```


##### For each population, run PCA analysis for unrelated samples


```
    sos run PCA.ipynb flashpca \
        --container container/flashpcaR.sif
```


##### For each population, run projection of related individuals if applicable


```
    sos run PCA.ipynb project_samples \
      --container container/flashpcaR.sif
```



```

```


##### For each population do their own QC to finalize 


```
sos run GRM.ipynb grm \
    --container container/bioinfo.sif
```


##### Merge separated bed files into one


```
sos run pipeline/genotype_formatting.ipynb vcf_to_plink
    --container /mnt/vast/hpc/csg/containers/bioinfo.sif \

```



```
sos run xqtl-pipeline/pipeline/genotype_formatting.ipynb merge_plink \
    --container /mnt/vast/hpc/csg/containers/bioinfo.sif \

```


##### Genotype data partition by chromosome

Timing <1 min

```
sos run pipeline/genotype_formatting.ipynb genotype_by_chrom \
    --container containers/bioinfo.sif 
```


#### B.  Phenotype data preprocessing


```
sos run gene_annotation.ipynb annotate_coord_gene \
    --container container/rna_quantification.sif --phenotype-id-type gene_name
```


Timing <1 min

```
sos run pipeline/gene_annotation.ipynb annotate_coord_protein \
    --container containers/rna_quantification.sif --sep "," 
```


##### Partition by chromosome


```
sos run pipeline/phenotype_formatting.ipynb phenotype_by_chrom \
    --container containers/bioinfo.sif
```



```
sos run pipeline/phenotype_formatting.ipynb partition_by_chrom \
    --container containers/rna_quantification.sif
```



```
sos run pipeline/phenotype_formatting.ipynb partition_by_chrom \
    --container containers/rna_quantification.sif
```



```
sos run phenotype_imputation.ipynb flash \

```


#### C.  Covariate Data Preprocessing


```
sos run pipeline/PEER_factor.ipynb PEER \
   --container containers/PEER.sif  \

```



```

```



```

```



```

```


##### Compute residule on merged covariates and perform hidden factor analysis

Timing <1 min

```
sos run pipeline/covariate_hidden_factor.ipynb Marchenko_PC \
   --container containers/PCAtools.sif
```


##### APEX


```
sos run pipeline/BiCV_factor.ipynb BiCV \
   --container containers/apex.sif  \

```


##### Merge Covariates and Genotype PCA

Timing <1 min

```
sos run pipeline/covariate_formatting.ipynb merge_genotype_pc \
    --container containers/bioinfo.sif
```



## Timing




| Step(Major Section) | Substep(Miniprotocol) | Time|
|------|-----|----|
|Reference data| Reference Data| ~4 hours|
|Molecular Phenotypes| RNA-seq expression| <3.5 hours|
| | Alternative splicing from RNA-seq data| <2 hours|
|Data Pre-processing| Genotype data preprocessing| < X minutes|
| | Phenotype data preprocessing| < X minutes|
| | Covariate Data Preprocessing| < X minutes|

## Troubleshooting




## Anticipated Results




#### A.  Reference Data

Our pipeline uses the following reference data for RNA-seq expression quantification:

#### A.  RNA-seq expression

The final output contained QCed and normalized expression data in a bed.gz file. This file is ready for use in TensorQTL.
#### B.  Alternative splicing from RNA-seq data

The final output contains the QCed and normalized splicing data from leafcutter and psichonics.
#### A.  Genotype data preprocessing

#### B.  Phenotype data preprocessing

# Phenotype data preprocessing

#### C.  Covariate Data Preprocessing

# Covariate Data Preprocessing


## Figures




## Tables




## Supplementary Information




## Author Contributions Statements



## Acknowledgements



## Competing Interests



## References



1. Signal et al. 2022. https://doi.org/10.1186/s12859-022-04572-7 
2. DeLuca et al. 2012. https://doi.org/10.1093/bioinformatics/bts196 
3. Wang et al. 2008. https://doi.org/10.1038/nature07509 
4. Park et al. 2018. https://doi.org/10.1016/j.ajhg.2017.11.002 
5. Li et al. 2018. https://doi.org/10.1038/s41588-017-0004-9 
6. Agostinho et al. 2019. https://doi.org/10.1093/nar/gky888 

## Keywords




