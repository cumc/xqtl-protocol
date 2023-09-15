# This needs pgenlibr package
# devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")

# read PLINK files
read_pvar <- function(pgen){
  pvarf <- paste0(tools::file_path_sans_ext(pgen), ".pvar")
  pvardt <- data.table::fread(pvarf, skip = "#CHROM")
  pvardt <- dplyr::rename(pvardt, "chrom" = "#CHROM", "pos" = "POS",
                "alt" = "ALT", "ref" = "REF", "id" = "ID")
  pvardt <- pvardt[, c("chrom", "id", "pos", "alt", "ref")]
  return(pvardt)
}

read_bim <- function(bed) {
  bimf <- paste0(tools::file_path_sans_ext(bed), ".bim")
  bim <- data.table::fread(bimf)
  colnames(bim) <- c("chrom", "id", "gpos", "pos", "a1", "a0")
  return(bim)
}

read_psam <- function(pgen) {
  psamf <- paste0(tools::file_path_sans_ext(pgen), ".psam")
  psam = data.table::fread(psamf, header=T)
  colnames(psam)[1:2] = c("FID", "IID")
  return(psam)
}

read_fam <- function(bed) {
    famf <- paste0(tools::file_path_sans_ext(bed), ".fam")
    return(data.table::fread(famf, header = F))
}

# open pgen/pvar PLINK 2 data format
open_pgen <- function(pgenf){
    return(pgenlibr::NewPgen(pgenf))
} 

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
open_bed <- function(bed){
    raw_s_ct <- nrow(read_fam(bed))
    return(pgenlibr::NewPgen(bed, raw_sample_ct = raw_s_ct))
}

read_pgen <- function(pgen, variantidx = NULL, meanimpute = F ) {
  if (is.null(variantidx)){
    variantidx <- 1: pgenlibr::GetVariantCt(pgen)}

  pgenlibr::ReadList(pgen,
                     variant_subset = variantidx,
                     meanimpute = meanimpute)
}

## genotype matrix preprocesing
compute_maf <- function(geno){
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x))==1) return(T)
  else return(F)
}

filter_X <- function(X, missing_rate_thresh, maf_thresh) {
    rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, compute_maf) < maf_thresh)
    if (length(rm_col)) X <- X[, -rm_col]
    rm_col <- which(apply(X, 2, is_zero_variance))
    if (length(rm_col)) X <- X[, -rm_col]
    return(mean_impute(X))
}


load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names 
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           imiss_cutoff = 0,
                                           y_as_matrix = FALSE) {
    library("dplyr")
    library("readr")
    library("stringr")
    library("purrr")
 
    ### Load genotype
    geno = plink2R::read_plink(genotype)
    rownames(geno$bed) = read.table(text = rownames(geno$bed), sep= ":" )$V2
    ### FIXME: complete this -- let's just summarize genotype data is enough
    print("Dimension of input genotype data is")

    ### Load phenotype and covariates
    covar_list = covariate%>%mutate(covar = map(path, ~read_delim(.x,"\t")%>%select(-`#id`)%>%na.omit%>%t()))

    # FIXME: 
    ### First, create the phenotype Y matrix
    ### with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column
    ### NOTE: natually, missing data in Y should be allowed because as a matrix, some columns of Y may not have data for some samples and that is totally okay
    traits = phenotype_list$tissue
    
    ### Secondly, find matching samples between phenotype Y matrix, list of covariate matrices Z, and genotype X matrix, 
    ### and actually match the samples by name which should be the rownames
    ### NOTE: it has to be a 3 way matching for all of the X, Y and Z.
    ### Each column of Y should match the corresponding Z matrix -- if the sample in the column does not exist in the corresponding Z matrix we have to set it NA; then we only keep the rows of Y that has at least non-NA values
    ### Hao, it is not clear to me if it is the case below because am not sure if the 3rd `mutate` assigns the variables sequentially
    ### So I'm not sure if the logic is 100% correct
    phenotype_list = inner_join(phenotype,covar_list, by = "tissue")
    phenotype_list = phenotype_list %>% mutate(flag = map(path.x, read_gene_pheno)) %>% dplyr::filter(lengths(flag) >1 ) %>% mutate(Y = map2(path.x,covar, ~read_gene_pheno(.x)%>%select(-c(`#chr`,start,end,${id_name}),      rownames(.y))%>%t%>%as.matrix))%>%
                                    mutate(covar = map2(covar, Y , ~.x[intersect(.x%>%rownames,rownames(.y)),] ), # add a flag so that those returning NA can be filtered out
                                            X_data = map(covar,~X[intersect(rownames(.x),rownames(X)),]),
                                            covar = map2(covar, X_data , ~.x[intersect(.x%>%rownames,rownames(.y)),] ),
                                            Y = map2(covar, Y ,~ .y[intersect(.x%>%rownames,rownames(.y)),] ),
                                            #### Get residue for each of tissue
                                            Y_resid = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale%>%t%>%as_tibble))
    ### Next, regress covariate from each column of Y separately, take residual and center/scale it
    Y_resid = phenotype_list%>%select(Y_resid)%>%tidyr::unnest(Y_resid)%>%t%>%as.matrix # in this part unnest will broad those matrix at smalled dimension (fewer rows/columns), so we should leave out those rows/ columns
    colnames(Y_resid) = phenotype_list$tissue
    ### Next, filter X matrix based on the remaining sample genotypes
    N = nrow(X)
    maf_cutoff = max(mac_cutoff/(2.0*N), maf_cutoff)
    X = filter_X(X, imiss_cutoff, maf_cutoff)
    ### FIXME: complete this -- let's just summarize genotype data is enough
    print("Dimension of filtered genotype data is")
    ### Finally, regress covariate from X to create a list of X residual matrices
    ### Note: here each 
    X_list = phenotype_list%>%mutate( X_resid = map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale))%>%pull(X_resid)
    ## FIXME: return quantities should be these:
    return (list(
            residual_X_scaled =, # is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names
            residual_Y_scaled =, # if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X_scaled in terms of which samples are missing.
            X = # is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y_scaled matrix form (not list); but the matrices inside residual_X_scaled would be subsets of sample name of residual_Y_scaled matrix form (not list).
            ))
}