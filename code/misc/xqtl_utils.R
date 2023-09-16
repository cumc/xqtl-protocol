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

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle  <- "--file="
  match   <- grep(needle,cmdArgs)
  if (length(match) > 0) {
    ## Rscript
    path <- cmdArgs[match]
    path <- gsub("\\~\\+\\~", " ", path)
    return(normalizePath(sub(needle, "", path)))
  } else {
    ## 'source'd via R console
    return(sys.frames()[[1]]$ofile)
  }
}

load_script <- function() {
  fileName <- thisFile()
  return(ifelse(!is.null(fileName) && file.exists(fileName),
                readChar(fileName,file.info(fileName)$size),""))
}
             
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names 
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           imiss_cutoff = 0,
                                           y_as_matrix = FALSE
                                          ) {
       library("plink2R") # This is to add Rcpp and other dependency
    library("dplyr")
    library("readr")
    library("stringr")
    library("purrr")

    ### Load genotype
    geno = read_plink(genotype)
    rownames(geno$bed) = read.table(text = rownames(geno$bed), sep= ":" )$V2
    
    
     ### Load phenotype and covariates
    phenotype_list = tibble(covariate_path = covariate ,phenotype_path =phenotype ) %>%
        mutate(covar = map(covariate_path, ~read_delim(.x,"\t")%>%select(-1)%>%na.omit%>%t()),
              Y = map2(phenotype_path,covar, ~read_delim(.x,"\t")%>%select(-4)%>%select(rownames(.y))%>%t()%>%as.matrix), # Y and covar 1st match
              Y = map(Y, ~.x%>%na.omit),    # remove na where Y raw data has na which block regression
              dropped_sample = map2(covar, Y , ~rownames(.x)[!rownames(.x) %in% rownames(.y)])  ,
              covar = map2(covar, Y , ~.x[intersect(.x%>%rownames,rownames(.y)),]), # remove the dropped samples from Y
              X_data = map(covar,~ filter_X( geno$bed[intersect(rownames(.x),rownames(geno$bed)),], imiss_cutoff, max(maf_cutoff, mac_cutoff/(2*length(intersect(rownames(.x),rownames(geno$bed))) ) ))   ))   # Y ( cov )  and specific X and covar match, filter X variants based on the overlapped samples.
              #### Get residue Y for each of condition

   phenotype_list = phenotype_list%>%mutate(Y_resid = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%
                                                           scale%>%t%>%as_tibble)) ## T so that it can be unnest
    if(y_as_matrix){
    Y_resid = phenotype_list%>%select(Y_resid)%>%tidyr::unnest(Y_resid)%>%t%>%as.matrix
    }else{Y_resid = map(phenotype_list$Y_resid,~.x%>%t) # Transpose back 
    }
    
    maf_cutoff = max(maf_cutoff,mac_cutoff/(2*nrow(geno$fam)))
    X = filter_X(geno$bed, imiss_cutoff, maf_cutoff) ## Filter X for mvSuSiE
    print(paste0("Dimension of input genotype data is row:", nrow(X), " column: ",ncol(X) ))
    
    X_list = phenotype_list%>%mutate( X_resid = map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale))%>%pull(X_resid)
    return (list(
            residual_Y_scaled = Y_resid, # if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X_scaled in terms of which samples are missing.
            residual_X_scaled = X_list , # is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names
            X = X ,  # is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y_scaled matrix form (not list); but the matrices inside residual_X_scaled would be subsets of sample name of residual_Y_scaled matrix form (not list).
            dropped_sample = phenotype_list$dropped_sample , # A documentation of dropped samples
            traits = phenotype_list$phenotype_path
            ))
    }
               
    # FIXME: 
    ### First, create the phenotype Y matrix
    ### with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column
    ### NOTE: natually, missing data in Y should be allowed because as a matrix, some columns of Y may not have data for some samples and that is totally okay
    
    
    ### Secondly, find matching samples between phenotype Y matrix, list of covariate matrices Z, and genotype X matrix, 
    ### and actually match the samples by name which should be the rownames
    ### NOTE: it has to be a 3 way matching for all of the X, Y and Z.
    ### Each column of Y should match the corresponding Z matrix -- if the sample in the column does not exist in the corresponding Z matrix we have to set it NA; then we only keep the rows of Y that has at least non-NA values
    ### Hao, it is not clear to me if it is the case below because am not sure if the 3rd `mutate` assigns the variables sequentially
    ### So I'm not sure if the logic is 100% correct
   
    ### Next, regress covariate from each column of Y separately, take residual and center/scale it
       ### Next, filter X matrix based on the remaining sample genotypes

    ### FIXME: complete this -- let's just summarize genotype data is enough
    ### Finally, regress covariate from X to create a list of X residual matrices
    ## FIXME: return quantities should be these: