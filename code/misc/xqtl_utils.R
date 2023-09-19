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
             
compute_cov_flash <- function(Y, error_cache = NULL){
    covar <- diag(ncol(Y))
    tryCatch({
    fl <- flashier::flash(Y, var.type = 2, prior.family = c(flashier::prior.normal(), flashier::prior.normal.scale.mix()), backfit = TRUE, verbose.lvl=0)
    if(fl$n.factors==0){
      covar <- diag(fl$residuals.sd^2)
    } else {
      fsd <- sapply(fl$fitted.g[[1]], '[[', "sd")
      covar <- diag(fl$residuals.sd^2) + crossprod(t(fl$flash.fit$EF[[2]]) * fsd)
    }
    if (nrow(covar) == 0) {
      covar <- diag(ncol(Y))
      stop("Computed covariance matrix has zero rows")
    }
    }, error = function(e) {
      if (!is.null(error_cache)) {
        saveRDS(list(data=Y, message=warning(e)), error_cache)
        warning("FLASH failed. Using Identity matrix instead.")
        warning(e)
      } else {
        stop(e)
      }
    })
    s <- apply(Y, 2, sd, na.rm=T)
    if (length(s)>1) s = diag(s)
    else s = matrix(s,1,1)
    covar <- s%*%cov2cor(covar)%*%s
    return(covar)
}

compute_cov_diag <- function(Y){
    covar <- diag(apply(Y, 2, var, na.rm=T))
    return(covar)
}
             

read_pheno <- function(file,region){
data.table::fread(cmd = paste0("tabix -h ",file ," ",region))%>%as_tibble() 
}





             
load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names 
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           imiss_cutoff = 0,
                                           y_as_matrix = FALSE) {
    library("plink2R")
    library("dplyr")
    library("readr")
    library("stringr")
    library("purrr")

    ## Load genotype
    geno = read_plink(genotype)
    rownames(geno$bed) = read.table(text = rownames(geno$bed), sep= ":" )$V2
    
    ## Load phenotype and covariates and perform some pre-processing
    ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
    phenotype_list = tibble(covariate_path = covariate, phenotype_path =phenotype) %>%
        mutate(covar = map(covariate_path, ~read_delim(.x,"\t")%>%select(-1)%>%na.omit%>%t()),
              Y = map2(phenotype_path,covar, ~read_pheno(.x,region )%>%select(-4)%>%select(rownames(.y))%>%t()%>%as.matrix), # Y and covar 1st match
              Y = map(Y, ~.x%>%na.omit),    # remove na where Y raw data has na which block regression
              dropped_sample = map2(covar, Y , ~rownames(.x)[!rownames(.x) %in% rownames(.y)])  ,
              covar = map2(covar, Y , ~.x[intersect(.x%>%rownames,rownames(.y)),]), # remove the dropped samples from Y
              X_data = map(covar,~ filter_X( geno$bed[intersect(rownames(.x),rownames(geno$bed)),], imiss_cutoff, max(maf_cutoff, mac_cutoff/(2*length(intersect(rownames(.x),rownames(geno$bed))) ) ))   ))   
              
    ## Get residue Y for each of condition
    phenotype_list = phenotype_list%>%mutate(Y_resid = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%
                                                           scale%>%t%>%as_tibble)) ## T so that it can be unnest
    if(y_as_matrix){
        Y_resid = phenotype_list%>%select(Y_resid)%>%tidyr::unnest(Y_resid)%>%t%>%as.matrix
        colnames(Y_resid) = conditions
        print(paste("Dimension of Y matrix:", nrow(Y_resid), ncol(Y_resid)))
    } else {
        Y_resid = map(phenotype_list$Y_resid,~.x%>%t) # Transpose back 
        names(Y_resid) = conditions
    }
    X = filter_X(geno$bed, imiss_cutoff, maf_cutoff) ## Filter X for mvSuSiE
    ## Get residue X for each of condition
    print(paste0("Dimension of input genotype data is row:", nrow(X), " column: ", ncol(X) ))
    X_list = phenotype_list%>%mutate( X_resid = map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale))%>%pull(X_resid)
    ## residual_Y_scaled: if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X_scaled in terms of which samples are missing.
    ## residual_X_scaled: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
    ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y_scaled matrix form (not list); but the matrices inside residual_X_scaled would be subsets of sample name of residual_Y_scaled matrix form (not list).
    return (list(
            residual_Y_scaled = Y_resid,
            residual_X_scaled = X_list,
            X = X,
            dropped_sample = phenotype_list$dropped_sample  # keep track of dropped samples due to a lack of overlap btw X,Y,Z
            ))
}