# Specific information on sections that need work in xQTL nature protocol paper

## Need Gao to complete or requires his input:
 * "comparisons with other methods" part of the introduction. See https://www.nature.com/nprot/for-authors/preparing-your-submission.
 * do we want to remove parts that are described as obselete (biomaRt) in gene_annotation.ipynb?
 * do we still need a notebook for the production of LD blocks and reference panel in the reference data section. Or is this what ld_prune_reference.ipynb is for?


## Parts that need work (marked with "FIXME" in the notebook):
 * reference_data.ipynb - code or reference to notebook for production of LD blocks and reference panel
 * RNA_calling.ipynb fastqc error - issue with final output file name. 
 * RNA_calling.ipynb rsem_call error
 * covariate_preprocessing.ipynb - timing
 * covariate_hidden_factor.ipynb - Marchenko PC citation and may need updated links to examples
 * phenotype_preprocessing.ipynb - timing and parts relating to phenotype_imputation.ipynb
 * phenotype_formatting.ipynb - updates needed to pipeline. also needs timing and information on the other commands
 * phenotype_imputation.ipynb - update notebook once container is ready
 * gene_annotation.ipynb - timing and list out other commands in the steps