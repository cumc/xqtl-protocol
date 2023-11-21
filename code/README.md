# xQTL protocol core modules

## Basic QTL analysis

1. `molecular_phenotypes`: quantification and QC of molecular phenotypes
2. `data_preprocessing`: prepare (including QC) $X$, $Y$, $Z$ for association analysis
3. `association_scan`: xQTL calling

## Advanced QTL analysis and integration

1. `cis_analysis`: SuSiE, fSuSiE, mvSuSiE, ColocBoost (for mvSuSiE and ColocBoost input can be QTL individual level data, xQTL summary statistics data and GWAS summary statistics), univariate TWAS weights, mr.mash, mmQTL (individual level data)
2. `trans_analysis`: SuSiE_RSS (genome-wide, including GWAS), polyfun, MRAID
3. `multivariate_genome`: MASH, METAL (genome-wide multivariate analysis)
4. `enrichment`: GREGOR, sLDSC, and a customized implementation using VEP
5. `pecotmr_analysis`: pair-wise enrichment, colocalization, TWAS and MR using the pecotmr framework
6. `rare_xqtl`: watershed pipeline

## Utilities

1. `command_generator`: an attempt to make a "one-click analysis" possible
2. `post_processing`: formatting of analysis results
3. `misc`: various scripts for maintenance and management

## Others

1. `quantile_qtl`: an experiment of applying quantile regression to xQTL studies
2. `integrative_prototype`: approaches we have tried but decided to not pursue further 
