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


Quality control and normalization are performed on output from the leafcutter and psichomics tools. The raw output data is first converted to bed format. Quality control involves the removal of features with high missingness across samples (default 40%) and any remaining missing values are retained. Then introns with less than a minimal variation (default of 0.005) are removed from the data. Quantile-Quantile normalization is performed on the quality controlled data. Imputation of missing values is done separately afterwards. 
#### Data Pre-processing (Step 3)
##### A.  Genotype data preprocessing


A major challenge in biomedical research is the quality control (QC) of sequencing data. False positive variant calls can hinder the ability to detect disease associated variants or introduce spurious associations, therefore the need for a rigorous QC. Our pipeline focuses on QC after the variant calling stage and requires project Variant Calling Format (pVCF) as input files. We have defined default thresholds for genotype and variant-level hard filtering based on recommendations from the UK Biobank team and a thorough review of the literature [[cf. Carson et al. BMC Bioinformatics (2014)](https://doi.org/10.1186/1471-2105-15-125),[cf. Lek et al. Nature (2016)](https://doi.org/10.1038/nature19057),[cf. Szustakowski et al. Nature Genetics (2021)](https://doi.org/10.1038/s41588-021-00885-0)]. Bcftools is used in our QC steps. We first handle multi-allelic sites by splitting them into bi-allelic records. We include an optional workflow to keep only bi-allelic sites in the data. Variants are then annotated based on dbSNP data. Genotypes are kept if they have a Genotype Depth (DP) >= 10 and a Genotype Quality (GQ) >= 20. Variants are included if at least one sample has an allelic balance (AB) >= 0.15 for Single Nucleotide Variants (SNVs) and AB>=0.2 for indels, variant missigness is below 20% and Hardy-Weinberg Equilibrium p-value is > 1e-08. Allele balance is calculated for heterozygotes as the number of bases supporting the least-represented allele over the total number of base observations. Output summary statistics, such as transistion/transversion ratios (TS/TV ratio) are calculated to determine the effectiveness of QC. 

We include steps for the formatting of genotype files. This includes the conversion between VCF and PLINK formats, the splitting of data (by specified input, chromosomes or genes) and the merging of data (by specified input, or by chromosomes).
##### B.  Phenotype data preprocessing


We use a gene coordinate annotation pipeline based on [`pyqtl`, as demonstrated here](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/eqtl_prepare_expression.py). This adds genomic coordinate annotations to gene-level molecular phenotype files generated in `gct` format and converts them to `bed` format for downstreams analysis.


Empirical Bayes Matrix Factorization (EBMF) is the primary imputation method we use to impute molecular phenotype data [[cf. Qi et al., medRxiv,  2023](https://doi.org/10.1101/2023.11.29.23299181)]. We also provide a collection of other methods to impute missing omics data values including missing forest, XGboost, k-nearest neighbors, soft impute, mean imputation, and limit of detection.

We include a collection of workflows to format molecular phenotype data. These include workflows to separate phenotypes by chromosome, by user-provided regions, a workflow to subset bam files and a workflow to extract samples from phenotype files.

##### C.  Covariate Data Preprocessing


Our covariate preprocessing steps merge genotypic principal components and fixed covariate files into one file for downstream QTL analysis. 

We provide two different procedures for hidden factor analysis from omics data in our pipeline. The first is the [Probabilistic Estimation of Expression Residuals (PEER) method](https://github.com/PMBio/peer/wiki/Tutorial), a method also used for GTEx eQTL data analysis. The second, and the one use for our main analyses, is a PCA based approach with automatic determination of the number of factors to use. This is mainly inspired by a recent benchmark from Jessica Li's group [[cf. Zhou et al., Genome Biology, 2022](https://doi.org/10.1186/s13059-022-02761-4)]. Please note that additional considerations should be taken for single-cell eQTL analysis as investigated by [[cf. Xue et al., Genome Biology, 2023](https://doi.org/10.1186/s13059-023-02873-5)].
#### QTL Association Testing (Step 4)
##### A.  QTL Association Analysis


We perform QTL association testing using TensorQTL [[cf. Taylor-Weiner et al (2019)](https://doi.org/10.1186/s13059-019-1836-7)]. An additional protocol was added to test for quantile QTL associations.
#### Advanced cis-QTL Analysis (Step 5)

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
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
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

##### i. Quality Control
Timing X min

```
sos run VCF_QC.ipynb qc    \
    --container oras://ghcr.io/cumc/bioinfo_apptainer:latest
```


##### ii. VCF to PLINK
Timing X min

```
sos run pipeline/genotype_formatting.ipynb vcf_to_plink
    --container /mnt/vast/hpc/csg/containers/bioinfo.sif \

```


##### iii. Perform QC on both rare and common variants

```
sos run xqtl-protocol/pipeline/GWAS_QC.ipynb qc_no_prune \
   --container /mnt/vast/hpc/csg/containers/bioinfo.sif \

```


##### iv. Partition Genotype Data by Chromosome
Timing <1 min

```
sos run pipeline/genotype_formatting.ipynb genotype_by_chrom \
    --container containers/bioinfo.sif 
```


##### v. Sample match with genotype
Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb genotype_phenotype_sample_overlap \
        --container containers/bioinfo.sif \

```


##### vi. Kinship QC

Timing <2 min

```
sos run pipeline/GWAS_QC.ipynb king \
    --container containers/bioinfo.sif \

```


##### vii. Prepare unrelated individuals data for PCA

Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb qc \
   --container containers/bioinfo.sif \

```


##### viii. PCA analysis for unrelated samples
Timing <2 min

```
sos run pipeline/PCA.ipynb flashpca \
   --container containers/flashpcaR.sif \

```


#### B.  Phenotype data preprocessing

##### i. Cooridnate Annotation
Timing <1 min

```
!sos run gene_annotation.ipynb annotate_coord_gene \
    --container  oras://ghcr.io/cumc/rna_quantification_apptainer:latest --phenotype-id-type gene_name \

```


Timing <1 min

```
!sos run gene_annotation.ipynb annotate_coord_protein \
    --container  oras://ghcr.io/cumc/rna_quantification_apptainer:latest --sep ","  \

```


##### iii. Partition by chromosome
Timing < 1 min

```
!sos run phenotype_formatting.ipynb phenotype_by_chrom \
    --container oras://ghcr.io/cumc/bioinfo_apptainer:latest \

```


#### C.  Covariate Data Preprocessing

##### i. Merge Covariates and Genotype PCs

Timing <1 min

```
sos run pipeline/covariate_formatting.ipynb merge_genotype_pc \
    --container containers/bioinfo.sif
```


##### ii. Generate Hidden Factors
Timing X min

```
!sos run covariate_hidden_factor.ipynb PEER \
   --container oras://ghcr.io/cumc/factor_analysis_apptainer:latest \

```


Timing <1 min

```
!sos run covariate_hidden_factor.ipynb Marchenko_PC \
    --container oras://ghcr.io/cumc/pcatools_apptainer:latest \

```


### 4. QTL Association Testing
#### A.  QTL Association Analysis

##### i. Cis TensorQTL Command 

```
sos run pipeline/TensorQTL.ipynb cis \
    --container containers/TensorQTL.sif 
```


##### ii. Trans TensorQTL Command 


```
sos run xqtl-protocol/pipeline/TensorQTL.ipynb trans \
    --container containers/TensorQTL.sif 
```


### 5. Advanced cis-QTL Analysis

## Timing




| Step(Major Section) | Substep(Miniprotocol) | Time|
|------|-----|----|
|Reference data| Reference Data| ~4 hours|
|Molecular Phenotypes| RNA-seq expression| <3.5 hours|
| | Alternative splicing from RNA-seq data| <2 hours|
|Data Pre-processing| Genotype data preprocessing| < X minutes|
| | Phenotype data preprocessing| < 12 minutes|
| | Covariate Data Preprocessing| < 3 minutes|
|QTL Association Testing| QTL Association Analysis| < X minutes|

## Troubleshooting




## Anticipated Results




####  Reference Data

Please refer to the protocol website for more information on these results.
####  RNA-seq expression

The final output contained QCed and normalized expression data in a bed.gz file. This file is ready for use in TensorQTL.
####  Alternative splicing from RNA-seq data

The final output contains the QCed and normalized splicing data from leafcutter and psichomics.
####  Genotype data preprocessing

####  Phenotype data preprocessing

# Phenotype data preprocessing

Phenotype preprocessing should result in a phenotype file formatted and ready for use in TensorQTL.
####  Covariate Data Preprocessing

Processed covariate data includes a file with covariates and hidden factors for use in TensorQTL.
####  QTL Association Analysis

TensorQTL will produce empirical and standardized cis/trans results.

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
7. Carson et al. 2014. https://doi.org/10.1186/1471-2105-15-125 
8. Lek et al. 2016. https://doi.org/10.1038/nature19057 
9. Szustakowski et al. 2021. https://doi.org/10.1038/s41588-021-00885-0 
10. Qi et al. 2023. https://doi.org/10.1101/2023.11.29.23299181 
11. Zhou et al. 2022. https://doi.org/10.1186/s13059-022-02761-4 
12. Xue et al. 2023. https://doi.org/10.1186/s13059-023-02873-5 
13. Taylor-Weiner et al. 2019. https://doi.org/10.1186/s13059-019-1836-7 

## Keywords




