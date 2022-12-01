#' Function to perform differential analysis using linear mixture models
#'
#' @param metadata Data frame containing the required metadata
#' @param provided.matrix Matrix or sparse matrix class to use as input. Overrides useMatrix
#' @param subset Character vector containing names of features to subset for differential analysis
#' @param sample.column Character vector of length one indicating the column containing the sample identity
#' @param cluster.column Character vector of length one indicating the column containing the cluster identity
#' @param selected.clusters Character vector indicating the clusters to be included in the analysis
#' @param treatment.column Character vector of length one indicating the column containing the treatments
#' @param treatment.levels Character vector of length two indicating which treatments to compare
#' @param ncores Integer indicating how many cores to use for parallel computing
#' @return Data frame containing the results of the differential analysis
#'

DiffLMM = function(metadata,
                   provided.matrix,
                   subset = NULL,
                   sample.column,
                   cluster.column,
                   selected.clusters,
                   treatment.column,
                   treatment.levels,
                   ncores = 1){

  if(class(metadata) != "data.frame"){
    stop("Specified metadata is not data.frame class")
  }

  if(!sample.column %in% colnames(metadata)){
    stop("Specified sample column not found in @cellColData slot")
  }

  if(!cluster.column %in% colnames(metadata)){
    stop("Specified pseudotime column not found in @cellColData slot")
  }

  if(!(treatment.column %in% colnames(metadata))){
    stop("Specified treatment column not found in @cellColData slot.")
  }

  if(!sum(treatment.levels %in% unique(metadata[,treatment.column])) == 2){
    stop("Specified treatment levels not found in specified treatment column in @cellColData slot")
  }

  if(sum(selected.clusters %in% unique(metadata[,cluster.column])) != length(selected.clusters)){
    message(paste0(setdiff(selected.clusters, unique(metadata[,cluster.column])), " not found in specified cluster column. Proceeding with: ",intersect(selected.clusters, unique(metadata[,cluster.column]))))
    if(!sum(selected.clusters %in% unique(metadata[,cluster.column])) == 0)
      stop("None of the selected clusters found in specified cluster.column")
  }


  m <- provided.matrix


  if(!is.null(subset)){
    if(nrow(m[rownames(m) %in% subset,]) == 0){
      stop("None of the specified features in subset is found in matrix")
    }else{
      m = m[rownames(m) %in% subset,]
    }
  }

  message("------- Running differential analysis based on linear mixture models (LMM) -------")
  out = mclapply(rownames(m), function(x){
    dev = as.vector(t(m[x,]))
    names(dev) = colnames(m)
    temp = metadata
    temp = temp[,c(sample.column, cluster.column, treatment.column)]
    temp[,x] = dev[rownames(temp)]
    temp = temp[temp[,treatment.column] %in% treatment.levels,]
    temp = temp[temp[,cluster.column] %in% selected.clusters,]
    temp = temp[complete.cases(temp),]

    temp = data.frame(motif = temp[,x], treatment = temp[,treatment.column], sample = temp[,sample.column])

    lmm1 = suppressMessages(lmer(data = temp, formula = motif ~ treatment + (1|sample), REML = F))
    lmm2 = suppressMessages(lmer(data = temp, formula = motif ~ (1|sample), REML = F))
    out = anova(lmm1,lmm2)

    pval = out$`Pr(>Chisq)`[2]
    delta = mean(temp$motif[temp$treatment == treatment.levels[2]], na.rm = T) - mean(temp$motif[temp$treatment == treatment.levels[1]], na.rm = T)
    fc = mean(temp$motif[temp$treatment == treatment.levels[2]], na.rm = T) / mean(temp$motif[temp$treatment == treatment.levels[1]], na.rm = T)

    final = c(pval = pval, delta = delta, fc = fc)
    return(final)
  }, mc.cores = ncores)

  names(out) = rownames(m)
  out = as.data.frame(do.call(rbind,out))
  out$fdr = p.adjust(out$pval, method = "fdr")
  out$feature = rownames(out)
  message("------- Done! -------")

  return(out)
}
