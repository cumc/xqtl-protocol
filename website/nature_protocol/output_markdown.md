## Title



## Authors




## Abstract




## Introduction




### Development of the protocol




### Applications of the method




### Comparison with other methods




### Experimental Design




#### Molecular Phenotypes (Step 1)
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
#### Data Pre-processing (Step 2)
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




### 1. Molecular Phenotypes
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
#no STAR wasp filter:
!sos run RNA_calling.ipynb STAR_align \
    --cwd ../../output_test/star_output \
    --samples ../../PCC_sample_list_subset \
    --data-dir /restricted/projectnb/amp-ad/ROSMAP_PCC_AC/PCC/ \
    --STAR-index ../../reference_data/STAR_Index/ \
    --gtf ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.ERCC.gtf \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    --reference-fasta ../../reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta \
    --ref-flat ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformated.ERCC.gtf.ref.flat \
    -s build -c ../csg.yml  -q neurology 
```



```
#include STAR wasp filter:
!sos run RNA_calling.ipynb STAR_align \
    --cwd ../../output_test/star_output_wasp \
    --samples ../../PCC_sample_list_subset \
    --data-dir /restricted/projectnb/amp-ad/ROSMAP_PCC_AC/PCC/ \
    --STAR-index ../../reference_data/STAR_Index/ \
    --gtf ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.ERCC.gtf \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    --reference-fasta ../../reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta \
    --ref-flat ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformated.ERCC.gtf.ref.flat \
    --varVCFfile ../../reference_data/ZOD14598_AD_GRM_WGS_2021-04-29_all.recalibrated_variants.leftnorm.filtered.AF.WASP.vcf \
    -s build  -c ../csg.yml  -q neurology 
```


##### iv. Call gene-level RNA expression via rnaseqc
Timing <30 min

```
!sos run RNA_calling.ipynb rnaseqc_call \
    --cwd ../../output_test/star_output \
    --samples ../../PCC_sample_list_subset \
    --data-dir /restricted/projectnb/amp-ad/ROSMAP_PCC_AC/PCC/ \
    --STAR-index ../../reference_data/STAR_Index/ \
    --gtf ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.ERCC.gtf \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    --reference-fasta ../../reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta \
    --ref-flat ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformated.ERCC.gtf.ref.flat \
    --bam_list ../../output_test/star_output/PCC_sample_list_subset_bam_list \
    -s build  -c ../csg.yml  -q neurology 
```


##### v. Call transcript level RNA expression via RSEM
Timing <X hours


```
!sos run RNA_calling.ipynb rsem_call \
    --cwd ../../output_test/star_output \
    --samples ../../PCC_sample_list_subset \
    --data-dir /restricted/projectnb/amp-ad/ROSMAP_PCC_AC/PCC/ \
    --STAR-index ../../reference_data/STAR_Index/ \
    --gtf ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.ERCC.gtf \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    --reference-fasta ../../reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta \
    --ref-flat ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformated.ERCC.gtf.ref.flat \
    --bam_list ../../output_test/star_output/PCC_sample_list_subset_bam_list \
    --RSEM-index ../../reference_data/RSEM_Index \
    -s build  -c ../csg.yml  -q neurology 
```


##### vi. Multi-sample RNA-seq QC
Timing <15min

```
!sos run bulk_expression_QC.ipynb qc \
    --cwd ../../output_test/rnaseqc_qc \
    --tpm-gct ../../output_test/star_output/PCC_sample_list_subset.rnaseqc.gene_tpm.gct.gz \
    --counts-gct ../../output_test/star_output/PCC_sample_list_subset.rnaseqc.gene_readsCount.gct.gz \
    --DSFilterPercent 0.1 \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    -s force -c ../csg.yml  -q neurology 
```


##### vii. Multi-sample read count normalization
Timing <10min

```
!sos run bulk_expression_normalization.ipynb normalize \
    --cwd ../../output_test/normalize \
    --tpm-gct ../../output_test/rnaseqc_qc/PCC_sample_list_subset.rnaseqc.low_expression_filtered.outlier_removed.tpm.gct.gz \
    --counts-gct ../../output_test/rnaseqc_qc/PCC_sample_list_subset.rnaseqc.low_expression_filtered.outlier_removed.geneCount.gct.gz \
    --annotation-gtf ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.ERCC.gtf \
    --sample-participant-lookup ../../PCC_sample_subset_map_test \
    --count-threshold 1 \
    --tpm_threshold 0.1 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm  \
    --container oras://ghcr.io/cumc/rna_quantification_apptainer:latest \
    -s force -c ../csg.yml -q neurology
```


#### B.  Alternative splicing from RNA-seq data

##### i. Splicing Quantification

Timing <30min

```
!sos run splicing_calling.ipynb leafcutter \
    --cwd ../../output_test/leafcutter \
    --samples ../../PCC_sample_list_subset_leafcutter \
    --data-dir ../../output_test/star_output_wasp \
    --container oras://ghcr.io/cumc/leafcutter_apptainer:latest \
    -c ../csg.yml -q neurology
```


Timing ~30min

```
!sos run splicing_calling.ipynb psichomics \
    --cwd ../../output_test/psichomics/ \
    --samples ../../PCC_sample_list_subset_leafcutter \
    --data-dir ../../output_test/star_output_wasp \
    --splicing_annotation ../../reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.ERCC.SUPPA_annotation.rds \
    --container oras://ghcr.io/cumc/psichomics_apptainer:latest \
    -c ../csg.yml -q neurology
```


##### ii. Splicing QC and Normalization
Timing  ~30min

```
!sos run splicing_normalization.ipynb leafcutter_norm \
    --cwd ../../output_test/leafcutter/normalize \
    --ratios ../../output_test/leafcutter/PCC_sample_list_subset_leafcutter_intron_usage_perind.counts.gz \
    --container oras://ghcr.io/cumc/leafcutter_apptainer:latest \
    -s force -c ../csg.yml -q neurology
```


Timing ~20min

```
!sos run splicing_normalization.ipynb psichomics_norm \
    --cwd ../../output_test/psichomics/normalize \
    --ratios ../../output_test/psichomics/psi_raw_data.tsv \
    --container oras://ghcr.io/cumc/psichomics_apptainer:latest \
    -c ../csg.yml -q neurology
```


### 2. Data Pre-processing
#### A.  Genotype data preprocessing


```
sos run VCF_QC.ipynb rename_chrs \
    --genoFile reference_data/00-All.vcf.gz \
    --cwd reference_data --container bioinfo.sif
```



```
sos run VCF_QC.ipynb dbsnp_annotate \
    --genoFile reference_data/00-All.add_chr.vcf.gz \
    --cwd reference_data --container bioinfo.sif
```



```
sos run VCF_QC.ipynb qc    \
    --genoFile data/MWE/MWE_genotype.vcf     \
    --dbsnp-variants data/reference_data/00-All.add_chr.variants.gz  \
    --reference-genome data/reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta   \
    --cwd MWE/output/genotype_1 --container bioinfo.sif -J 1 -c csg.yml -q csg
```



```
sos run VCF_QC.ipynb qc    \
    --genoFile data/mwe/mwe_genotype_list    \
    --dbsnp-variants data/reference_data/00-All.add_chr.variants.gz  \
    --reference-genome data/reference_data/GRCh38_full_analysis_set_plus_decoy_hla.noALT_noHLA_noDecoy_ERCC.fasta   \
    --cwd MWE/output/genotype_4 --container bioinfo.sif --add-chr
```



```
grep Ts/Tv MWE_genotype.leftnorm.known_variant.snipsift_tstv | rev | cut -d',' -f1 | rev
```



```
grep Ts/Tv MWE_genotype.leftnorm.filtered.*_variant.snipsift_tstv | rev | cut -d',' -f1 | rev
```



```
grep Ts/Tv MWE_genotype.leftnorm.novel_variant.snipsift_tstv | rev | cut -d',' -f1 | rev
grep Ts/Tv MWE_genotype.leftnorm.filtered.novel_variant.snipsift_tstv | rev | cut -d',' -f1 | rev
```


##### Perform QC on both rare and common variants

```
sos run xqtl-pipeline/pipeline/GWAS_QC.ipynb qc_no_prune \
   --cwd Genotype \
   --genoFile Genotype/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.bed \
   --geno-filter 0.1 \
   --mind-filter 0.1 \
   --hwe-filter 1e-08   \
   --mac-filter 0 \
   --container /mnt/vast/hpc/csg/containers/bioinfo.sif \
   -J 1 -q csg -c csg.yml --mem 150G
```


##### Sample match with genotype
Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb genotype_phenotype_sample_overlap \
        --cwd output/sample_meta \
        --genoFile input/protocol_example.genotype.chr21_22.fam  \
        --phenoFile input/protocol_example.protein.csv \
        --container containers/bioinfo.sif \
        --mem 5G
```


##### Kinship QC

Timing <2 min

```
sos run pipeline/GWAS_QC.ipynb king \
    --cwd output/kinship \
    --genoFile input/protocol_example.genotype.chr21_22.bed \
    --name pQTL \
    --keep-samples output/sample_meta/protocol_example.protein.sample_genotypes.txt \
    --container containers/bioinfo.sif \
    --no-maximize-unrelated \
    --mem 40G
```


##### Prepare unrelated individuals data for PCA

Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb qc \
   --cwd output/cache \
   --genoFile output/kinship/protocol_example.genotype.chr21_22.pQTL.unrelated.bed \
   --mac-filter 5 \
   --container containers/bioinfo.sif \
   --mem 16G
```


Timing <1 min

```
sos run pipeline/GWAS_QC.ipynb qc \
   --cwd output/cache \
   --genoFile input/protocol_example.genotype.chr21_22.bed \
   --keep-samples output/sample_meta/protocol_example.protein.sample_genotypes.txt \
   --name pQTL \
   --mac-filter 5 \
   --container containers/bioinfo.sif \
   --mem 40G
```



```
sos run GWAS_QC.ipynb qc_no_prune \
    --cwd output/genotype \
    --genoFile output/genotype/chr1_chr6.20220110.related.bed \
    --keep-variants output/genotype/chr1_chr6.20220110.unrelated.for_pca.filtered.prune.in \
    --maf-filter 0 --geno-filter 0 --mind-filter 0.1 \
    --name for_pca \
    --container container/bioinfo.sif
```


##### Estimate kinship in the sample


```
sos run GWAS_QC.ipynb king \
    --cwd output \
    --genoFile data/rename_chr22.bed \
    --kinship 0.13 \
    --name 20220110 \
    --container container/bioinfo.sif
```


##### Sample selection and QC the genotype data for PCA


```
sos run GWAS_QC.ipynb qc \
    --cwd output \
    --genoFile output/rename_chr22.20220110.unrelated.bed \
    --maf-filter 0.01 \
    --name for_pca \
    --container container/bioinfo.sif
```


##### PCA analysis for unrelated samples
Timing <2 min

```
sos run pipeline/PCA.ipynb flashpca \
   --cwd output/genotype_pca \
   --genoFile output/cache/protocol_example.genotype.chr21_22.pQTL.plink_qc.prune.bed \
   --container containers/flashpcaR.sif \
   --mem 16G
```



```
%preview ~/tmp/25-Jan-2022/output/pca/MWE_pheno.pca.scree.png
```


##### Projection of related individuals


```
```
sos run PCA.ipynb project_samples \
  --cwd output/pca \
  --genoFile output/rename_chr22.20220110.related.for_pca.filtered.extracted.bed \
  --phenoFile data/MWE_pheno.txt \
  --pca-model output/pca/MWE_pheno.pca.rds \
  --label-col RACE \
  --pop-col RACE \
  --maha-k 2 \
  --container container/flashpcaR.sif
```
```



```
%preview ~/tmp/25-Jan-2022/output/pca/MWE_pheno.pca.projected.pc.png
```


##### Finalize genotype QC by PCA for homogenous population


```
sos run GWAS_QC.ipynb qc_no_prune \
    --cwd output \
    --genoFile output/rename_chr22.20220110.unrelated.bed \
    --remove-samples output/pca/MWE_pheno.pca.projected.outliers \
    --name no_outlier \
    --container container/bioinfo.sif
```



```
sos run GWAS_QC.ipynb qc_no_prune \
    --cwd output \
    --genoFile output/rename_chr22.20220110.related.bed \
    --remove-samples output/pca/MWE_pheno.pca.projected.outliers \
    --keep-variants output/rename_chr22.20220110.unrelated.no_outlier.filtered.bim \
    --maf-filter 0 --geno-filter 0 --mind-filter 0.1 --hwe-filter 0 \
    --name no_outlier \
    --container container/bioinfo.sif
```



```
sos run genotype_formatting.ipynb merge_plink \
    --genoFile output/rename_chr22.20220110.unrelated.no_outlier.filtered.bed \
               output/rename_chr22.20220110.related.no_outlier.filtered.extracted.bed \
    --cwd output/genotype_final \
    --name chr22_20220110_qced \
    --container container/bioinfo.sif
```


##### Split data by population


```
pheno = read.table("data/MWE_pheno.txt", header = TRUE, stringsAsFactors=F)
for (i in 1:3){
 race = subset(pheno, RACE == i)
 race_id = cbind(race[,2],race[, 2])
 write.table(race_id, paste0("output/ID.", "race", i), quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)
}
```


##### For each population, do variant level and sample level QC on unrelated individuals, in preparation for PCA analysis


```
for i in race1 race3; do
    sos run GWAS_QC.ipynb qc \
        --cwd output \
        --genoFile output/rename_chr22.20220110.unrelated.bed \
        --keep-samples output/ID.$i \
        --name for_pca_$i \
        --container container/bioinfo.sif
```


##### For each population, extract previously selected variants from related individuals in preparation for PCA, only applying missingness filter at sample level, if applicable


```
for i in race1 race3; do
    sos run GWAS_QC.ipynb qc_no_prune \
        --cwd output \
        --genoFile output/rename_chr22.20220110.related.bed \
        --keep-variants output/rename_chr22.20220110.unrelated.for_pca_$i.filtered.prune.in \
        --keep-samples output/ID.$i \
        --maf-filter 0 --geno-filter 0 --mind-filter 0.1 --hwe-filter 0 \
        --name for_pca_$i \
        --container container/bioinfo.sif -s force
```


##### For each population, run PCA analysis for unrelated samples


```
for i in race1 race3; do
    sos run PCA.ipynb flashpca \
        --name $i \
        --cwd output/pca \
        --genoFile output/rename_chr22.20220110.unrelated.for_pca_$i.filtered.prune.bed \
        --phenoFile data/MWE_pheno.txt \
        --label-col RACE \
        --pop-col RACE \
        --maha-k 2 \
        --k 5 \
        --container container/flashpcaR.sif
```


##### For each population, run projection of related individuals if applicable


```
for i in race3; do
    sos run PCA.ipynb project_samples \
      --name $i \
      --cwd output/pca \
      --genoFile output/rename_chr22.20220110.related.for_pca_$i.filtered.extracted.bed \
      --phenoFile data/MWE_pheno.txt \
      --pca-model output/pca/MWE_pheno.$i.pca.rds \
      --label-col RACE \
      --pop-col RACE \
      --maha-k 2 \
      --k 5 \
      --container container/flashpcaR.sif
```



```
%preview ~/tmp/19-Jan-2022/output/pca/MWE_pheno.race3.pca.projected.pc.png 
```


##### For each population do their own QC to finalize 


```
sos run GRM.ipynb grm \
    --cwd output \
    --genotype-list data/genotype/mwe_genotype.list \
    --container container/bioinfo.sif
```


##### Merge separated bed files into one


```
sos run pipeline/genotype_formatting.ipynb vcf_to_plink
    --genoFile `ls vcf_qc/*.leftnorm.bcftools_qc.vcf.gz` \
    --cwd Genotype/ \
    --keep_samples ./ROSMAP_sample_list.txt
    --container /mnt/vast/hpc/csg/containers/bioinfo.sif \
    -J 22 -q csg -c csg.yml --mem 120G
```



```
sos run xqtl-pipeline/pipeline/genotype_formatting.ipynb merge_plink \
    --genoFile `ls *.leftnorm.bcftools_qc.bed` \
    --name ROSMAP_NIA_WGS.leftnorm.bcftools_qc  \
    --cwd Genotype/ \
    --container /mnt/vast/hpc/csg/containers/bioinfo.sif \
    -J 5 -q csg -c csg.yml --mem 300G
```


##### Genotype data partition by chromosome

Timing <1 min

```
sos run pipeline/genotype_formatting.ipynb genotype_by_chrom \
    --genoFile input/protocol_example.genotype.chr21_22.bed \
    --cwd output \
    --chrom `cut -f 1 input/protocol_example.genotype.chr21_22.bim | uniq | sed "s/chr//g"` \
    --container containers/bioinfo.sif 
```


#### B.  Phenotype data preprocessing


```
sos run gene_annotation.ipynb annotate_coord_gene \
    --cwd output \
    --phenoFile data/MWE.pheno_log2cpm.tsv.gz \
    --annotation-gtf reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.gene.ERCC.gtf \
    --sample-participant-lookup data/sampleSheetAfterQC.txt \
    --container container/rna_quantification.sif --phenotype-id-type gene_name
```


Timing <1 min

```
sos run pipeline/gene_annotation.ipynb annotate_coord_protein \
    --cwd output/phenotype \
    --phenoFile input/protocol_example.protein.csv \
    --annotation-gtf reference_data/Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.ERCC.gtf \
    --phenotype-id-type gene_name \
    --sample-participant-lookup output/sample_meta/protocol_example.protein.sample_overlap.txt \
    --container containers/rna_quantification.sif --sep "," 
```


##### Partition by chromosome


```
sos run pipeline/phenotype_formatting.ipynb phenotype_by_chrom \
    --cwd output/phenotype_by_chrom \
    --phenoFile output/phenotype/protocol_example.protein.bed.gz \
    --chrom `for i in {21..22}; do echo chr$i; done` \
    --container containers/bioinfo.sif
```



```
sos run pipeline/phenotype_formatting.ipynb partition_by_chrom \
    --cwd output  \
    --phenoFile MWE.log2cpm.mol_phe.bed.gz \
    --region-list ROSMAP_PCC.methylation.M.renamed.region_list \
    --container containers/rna_quantification.sif
```



```
sos run pipeline/phenotype_formatting.ipynb partition_by_chrom \
    --cwd mQTL_perchrom  \
    --phenoFile ROSMAP_arrayMethylation_covariates.sesame.methyl.beta.sample_matched.bed_BMIQ.bed.filter_na.bed.softImputed.bed.gz \
    --region-list ROSMAP_PCC.methylation.M.renamed.region_list \
    --container containers/rna_quantification.sif
```



```
sos run phenotype_imputation.ipynb flash \
    --phenoFile ./proteomics/rosmap/test1.bed.gz \
    --cwd ./proteomics/rosmap/ \
    --mem 40G \
    --walltime 100h
```


#### C.  Covariate Data Preprocessing


```
sos run pipeline/PEER_factor.ipynb PEER \
   --cwd output \
   --phenoFile ALL.log2cpm.bed.chr12.mol_phe.bed.gz  \
   --container containers/PEER.sif  \
   --N 3
```



```
tree ./output
```



```
%preview ./output/MWE.Cov_PEER.PEER_diagnosis.pdf -s png
```



```
cat ./output/MWE.Cov_PEER.PEER.cov.stdout
```


##### Compute residule on merged covariates and perform hidden factor analysis

Timing <1 min

```
sos run pipeline/covariate_hidden_factor.ipynb Marchenko_PC \
   --cwd output/covariate \
   --phenoFile output/phenotype/protocol_example.protein.bed.gz  \
   --covFile output/covariate/protocol_example.samples.protocol_example.genotype.chr21_22.pQTL.plink_qc.prune.pca.gz \
   --mean-impute-missing \
   --container containers/PCAtools.sif
```


##### APEX


```
sos run pipeline/BiCV_factor.ipynb BiCV \
   --cwd output \
   --phenoFile ALL.log2cpm.bed.chr12.mol_phe.bed.gz  \
   --container containers/apex.sif  \
   --N 3
```


##### Merge Covariates and Genotype PCA

Timing <1 min

```
sos run pipeline/covariate_formatting.ipynb merge_genotype_pc \
    --cwd output/covariate \
    --pcaFile output/genotype_pca/protocol_example.genotype.chr21_22.pQTL.plink_qc.prune.pca.rds \
    --covFile  input/protocol_example.samples.tsv \
    --tol_cov 0.4  \
    --k `awk '$3 < 0.8' output/genotype_pca/protocol_example.genotype.chr21_22.pQTL.plink_qc.prune.pca.scree.txt | tail -1 | cut -f 1 ` \
    --container containers/bioinfo.sif
```



## Timing




| Step | Substep | Time|
|------|-----|----|
|Molecular Phenotypes| RNA-seq expression| <3.5 hours|
| | Alternative splicing from RNA-seq data| <2 hours|
|Data Pre-processing| Genotype data preprocessing| < X minutes|
| | Phenotype data preprocessing| < X minutes|
| | Covariate Data Preprocessing| < X minutes|

## Troubleshooting




## Anticipated Results




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




