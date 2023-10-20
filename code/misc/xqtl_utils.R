# This needs pgenlibr package
# devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
# and PLINK2R
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE); remotes::install_github('gabraham/plink2R', subdir='plink2R', ref='d74be015e8f54d662b96c6c2a52a614746f9030d')

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
    rm_col <- which(apply(X, 2, compute_maf) <= maf_thresh)
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
             
tabix_region <- function(file, region){
    data.table::fread(cmd = paste0("tabix -h ", file, " ", region))%>%as_tibble() 
}

load_regional_association_data <- function(genotype, # PLINK file
                                           phenotype, # a vector of phenotype file names 
                                           covariate, # a vector of covariate file names corresponding to the phenotype file vector
                                           region, # a string of chr:start-end
                                           conditions, # a vector of strings
                                           maf_cutoff = 0,
                                           mac_cutoff = 0,
                                           imiss_cutoff = 0,
                                           y_as_matrix = FALSE,
                                           indel = TRUE) {
    library(plink2R)
    library(dplyr)
    library(readr)
    library(stringr)
    library(purrr)

    ## Load genotype
    geno = read_plink(genotype)
    rownames(geno$bed) = read.table(text = rownames(geno$bed), sep= ":")$V2
    ## if indel is true, remove the indel in the genotype
    if (indel==TRUE){
    geno_bim = geno$bim%>%rename("chrom" = "V1","variant_id" = "V2","alt" = "V5","ref"="V6")%>%mutate(indel = ifelse(grepl("[^ATCG]",alt)=="TRUE"|grepl("[^ATCG]",ref)=="TRUE"|nchar(alt)!=nchar(ref),1, 0))
    geno_bed = geno$bed[,geno_bim$indel==0]}
    else {
    geno_bed = geno$bed
    }
    ## Load phenotype and covariates and perform some pre-processing
    ### including Y ( cov ) and specific X and covar match, filter X variants based on the overlapped samples.
    data_list = tibble(covariate_path = covariate, phenotype_path =phenotype) %>%
        mutate(covar = map(covariate_path, ~read_delim(.x,"\t")%>%select(-1)%>%na.omit%>%t()),
        Y = map2(phenotype_path,covar, ~{
          y_data <- tabix_region(.x, region)%>%select(-4)%>%select(rownames(.y))%>%t()%>%as.matrix
          return(y_data)
          }),
        Y = map(Y, ~.x%>%na.omit),    # remove na where Y raw data has na which block regression
        dropped_sample = map2(covar, Y , ~rownames(.x)[!rownames(.x) %in% rownames(.y)]),
        covar = map2(covar, Y , ~.x[intersect(.x%>%rownames,rownames(.y)),]), # remove the dropped samples from Y
        X_data = map(covar,~ filter_X( geno_bed[intersect(rownames(.x),rownames(geno_bed)),], imiss_cutoff, max(maf_cutoff, mac_cutoff/(2*length(intersect(rownames(.x),rownames(geno_bed))) ) ))   ))
              
    ## Get residue Y for each of condition and its mean and sd
    data_list = data_list%>%mutate(Y_resid_mean = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%mean),
                               Y_resid_sd = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%sd),
                               Y_resid = map2(Y,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale%>%t%>%as_tibble)) ## T so that it can be unnest
    if(y_as_matrix) {
        Y_resid = data_list%>%select(Y_resid)%>%tidyr::unnest(Y_resid)%>%t%>%as.matrix
        colnames(Y_resid) = conditions
        print(paste("Dimension of Y matrix:", nrow(Y_resid), ncol(Y_resid)))
    } else {
        Y_resid = map(data_list$Y_resid,~.x%>%t) # Transpose back 
        names(Y_resid) = conditions
    }
    # Get X matrix for union of samples
    all_samples = map(data_list$covar, ~rownames(.x))%>%unlist%>%unique()
    maf_cutoff = max(maf_cutoff,mac_cutoff/(2*length(all_samples)))
    X = filter_X(geno_bed[all_samples,], imiss_cutoff, maf_cutoff) ## Filter X for mvSuSiE
    #
    maf_list = lapply(data_list$X_data, function(x) apply(x, 2, compute_maf))
    ## Get residue X for each of condition and its mean and sd
    print(paste0("Dimension of input genotype data is row:", nrow(X), " column: ", ncol(X) ))
    X_list = data_list%>%mutate(X_resid_mean= map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%data.frame()%>%apply(.,2,mean)),
                               X_resid_sd= map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%data.frame()%>%apply(.,2,sd)),
                               X_resid = map2(X_data,covar,~.lm.fit(x = cbind(1,.y), y = .x)$residuals%>%scale))%>%select(X_resid_mean,X_resid_sd,X_resid)
    ## residual_Y_scaled: if y_as_matrix is true, then return a matrix of R conditions, with column names being the names of the conditions (phenotypes) and row names being sample names. Even for one condition it has to be a matrix with just one column. if y_as_matrix is false, then return a list of y either vector or matrix (CpG for example), and they need to match with residual_X_scaled in terms of which samples are missing.
    ## residual_X_scaled: is a list of R conditions each is a matrix, with list names being the names of conditions, column names being SNP names and row names being sample names.
    ## X: is the somewhat original genotype matrix output from `filter_X`, with column names being SNP names and row names being sample names. Sample names of X should match example sample names of residual_Y_scaled matrix form (not list); but the matrices inside residual_X_scaled would be subsets of sample name of residual_Y_scaled matrix form (not list).
    return (list(
            residual_Y_scaled = Y_resid,
            residual_X_scaled = X_list$X_resid,
            residual_Y_sd = data_list$Y_resid_sd,
            residual_X_sd = X_list$X_resid_sd,
            dropped_sample = data_list$dropped_sample,
            covar = data_list$covar,
            Y = data_list$Y,
            X = X,
            X_data = data_list$X_data,
            maf = maf_list
            ))
}

load_regional_finemapping_data <- function(...) {
  dat <- load_regional_association_data(...)
  return (list(
          residual_Y_scaled = dat$residual_Y_scaled,
          residual_X_scaled = dat$residual_X_scaled,
          residual_Y_sd = dat$residual_Y_sd,
          residual_X_sd = dat$residual_X_sd,
          X = dat$X,
          dropped_sample = dat$dropped_sample,
          maf = dat$maf
          ))
}

load_regional_quantile_data <- function(...) {
  dat <- load_regional_association_data(...)
  return (list(
          Y = dat$Y,
          X_data = dat$X_data,
          covar = dat$covar,
          dropped_sample = dat$dropped_sample,
          maf = dat$maf
          ))
}

post_process_susie <- function(fobj, fdat, r, signal_cutoff = 0.7) {
    library(susieR)
    library(dplyr)
    get_cs_index <- function(snps_idx, fitted_data) {
        idx <- tryCatch(
            which(
                pmap(list(a = fitted_data$sets$cs), function(a) snps_idx %in% a) %>% unlist()
            ),
            error = function(e) NA_integer_
        )
        if(length(idx) == 0) return(NA_integer_)
        return(idx)
    }
    # collect results
    eff_idx = which(fobj$V>0)
    if (length(eff_idx)>0) {
        fobj$analysis_script = load_script()
        fobj$cs_corr = get_cs_correlation(fobj, X=fdat$residual_X_scaled[[r]])
        fobj$cs_snps = gsub("_",":",names(fobj$pip[unlist(fobj$sets$cs)]))
        fobj$phenotype_name = colnames(fdat$residual_Y_scaled[[r]])
        fobj$dropped_samples = fdat$dropped_sample[[r]]
        fobj$sample_names = rownames(fdat$residual_Y_scaled[[r]])
        fobj$variant_names = gsub("_",":",names(fobj$pip))
        variants_index = c(which(fobj$pip >= signal_cutoff), unlist(fobj$sets$cs)) %>% unique %>% sort
        if (length(variants_index)==0) {
            variants_index = which.max(fobj$pip)
        }
        maf = fdat$maf[[r]][variants_index]
        variants = gsub("_",":",names(fobj$pip)[variants_index])
        pip = fobj$pip[variants_index]
        cs_info = map_int(variants_index, ~get_cs_index(.x, fobj))
        cs_index = ifelse(is.na(cs_info), 0, str_replace(names(fobj$sets$cs)[cs_info], "L", "") %>% as.numeric)
        Y_resid_sd = fdat$residual_Y_sd[[r]]
        X_resid_sd = fdat$residual_X_sd[[r]]
        univariate_res = univariate_regression(fdat$residual_X_scaled[[r]][, variants_index, drop=F], fdat$residual_Y_scaled[[r]])
        fobj$top_loci = cbind(variants, maf, univariate_res$betahat*Y_resid_sd/X_resid_sd, univariate_res$sebetahat*Y_resid_sd/X_resid_sd, pip, cs_index)
        colnames(fobj$top_loci) = c("variant_id", "maf", "bhat", "sbhat", "pip", "cs_index")
        rownames(fobj$top_loci) = NULL
        # trim effects
        fobj$alpha = fobj$alpha[eff_idx,,drop=F]
        fobj$mu = fobj$mu[eff_idx,,drop=F]
        fobj$mu2 = fobj$mu2[eff_idx,,drop=F]
        fobj$V = fobj$V[eff_idx]
        # trim results
        fobj$Xr = NULL
        fobj$fitted = NULL
        colnames(fobj$lbf_variable) = NULL
        colnames(fobj$alpha) = NULL
        colnames(fobj$mu) = NULL
        colnames(fobj$mu2) = NULL
        names(fobj$X_column_scale_factors) = NULL
        names(fobj$pip) = NULL
        class(fobj) = "list"
    } else {
        fobj = load_script()
    }
    return(fobj)
}