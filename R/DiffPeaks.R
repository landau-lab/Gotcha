#' Function to perform differential analysis on peaks using mixed-effects logistic regression
#'
#' @param metadata Data frame containing the required metadata
#' @param provided.matrix Matrix or sparse matrix class to use as input. Overrides useMatrix
#' @param log.OR.threshold Limit testing to peaks which show at least an odds-ratio of X (log-scale) between the two groups of cells.
#' @param min.pct Only test peaks that are detected in a minimum fraction of min.pct cells in either of the two populations
#' @param subset Character vector containing names of features to subset for differential analysis
#' @param sample.column Character vector of length one indicating the column containing the sample identity
#' @param cluster.column Character vector of length one indicating the column containing the cluster identity
#' @param selected.clusters Character vector indicating the clusters to be included in the analysis
#' @param treatment.column Character vector of length one indicating the column containing the treatments
#' @param treatment.levels Character vector of length two indicating which treatments to compare
#' @param fixed.effect Character vector indicating the predictor variables in the model that are assumed to have the same effect on the response variable (treatment) for all the units (peaks) under study. Each variable should be a column name in the metadata dataframe.
#' @param random.effect Character vector indicating the predictor variables in the model that are assumed to have different effects on the response variable (treatment) for different units (peaks) under study. Each variable should be a column name in the metadata dataframe.
#' @param ncores Integer indicating how many cores to use for parallel computing
#' @return Data frame containing the results of the differential analysis
#'
#'
DiffPeaks = function(metadata,
                     provided.matrix,
                     log.OR.threshold = 0.1,
                     min.pct = 0.01,
                     subset = NULL,
                     sample.column,
                     cluster.column,
                     selected.clusters,
                     treatment.column,
                     treatment.levels,
                     fixed.effect,
                     random.effect,
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

  cells.to.use <- rownames(metadata)[metadata[,cluster.column]%in%selected.clusters]
  metadata <- metadata[cells.to.use,]
  m <- Signac::BinarizeCounts(provided.matrix[,cells.to.use])

  if(!is.null(subset)){
    if(nrow(m[rownames(m) %in% subset,]) == 0){
      stop("None of the specified features in subset is found in matrix")
    }else{
      m = m[rownames(m) %in% subset,]
    }
  }

  #Get cells from each group / treatment.levels
  group1=stats::na.omit(rownames(metadata)[metadata[,treatment.column]==treatment.levels[[1]]])
  group2=stats::na.omit(rownames(metadata)[metadata[,treatment.column]==treatment.levels[[2]]])


  #a = Number of exposed cases (MUT and peak accessible)
  #b = Number of exposed non-cases (MUT and peak non-accessible)
  #c = Number of unexposed cases (WT and peak accessible)
  #d = Number of unexposed non-cases (WT and peak non-accessible)
  ## Odds-Ratio = a/b/c/d = ad/bc

  #group2 = MUT, rowSums = nb of cells expressing the peak
  a = rowSums(m[,group2])
  b = length(group2)-a
  c = rowSums(m[,group1])
  d = length(group1)-c

  OR = ((a*d)+1) / ((b*c)+1)
  log.OR = log(OR)
  upper95CI = log.OR + 1.96*(sqrt(1/a + 1/b+ 1/c + 1/d))
  lower95CI = log.OR - 1.96*(sqrt(1/a + 1/b+ 1/c + 1/d))
  CI_include1 = ifelse(lower95CI<=1 & upper95CI>=1, TRUE, FALSE)

  #Compute number of cells expressing the Peak out of total number of cells in group

  pct.1 = rowSums(m[,group1])/length(group1)
  pct.2 = rowSums(m[,group2])/length(group2)

  #Compute per group average expression of Peaks
  avgExp.gr1 = rowMeans(m[,group1])
  avgExp.gr2 = rowMeans(m[,group2])
  avg_log2FC = log((avgExp.gr2 + 1/10000), base=2) - log((avgExp.gr1 + 1/10000), base = 2)

  features.meta = data.frame(pct.1=pct.1, pct.2=pct.2, avgExp.gr1=avgExp.gr1, avgExp.gr2=avgExp.gr2, avg_log2FC=avg_log2FC, log.OR=log.OR, lower95CI = lower95CI, upper95CI = upper95CI, CI_include1 = CI_include1, accessible_grp1 = c, accessible_grp2 = a)

  #Apply subsetting based on min.pct and log.OR.threshold - return filtered list of peaks
  subset.features = rownames(features.meta[features.meta$pct.1>min.pct & features.meta$pct.2>min.pct & abs(features.meta$log.OR) > log.OR.threshold & CI_include1==FALSE,])

  message(paste0("------- Running differential analysis based on general mixture models (GMM) on ", length(subset.features),  " Peaks -------"))

  print(mean(abs(features.meta$log.OR)))

  m <- m[subset.features,]

  out = mclapply(rownames(m), function(x){

    dev = as.vector(t(m[x,]))
    names(dev) = colnames(m)
    temp = metadata

    temp = temp[,c(c(sample.column, cluster.column, treatment.column),fixed.effect, random.effect)]
    temp[,x] = dev[rownames(temp)]
    temp = temp[temp[,treatment.column] %in% treatment.levels,]
    temp = temp[temp[,cluster.column] %in% selected.clusters,]
    temp = temp[complete.cases(temp),]

    temp = as.data.frame(temp)

    names(temp)[1] = 'Sample'
    names(temp)[3] = 'treatment'
    names(temp)[length(names(temp))] = 'peaks'

    temp$peaks = as.factor(temp$peaks)
    temp$treatment = as.factor(temp$treatment)
    temp$Sample = as.factor(temp$Sample)

    #scale fixed.effect
    temp[,fixed.effect]=scale(temp[,fixed.effect])


    fmla1 <- stats::as.formula(object = paste("treatment ~ peaks +",
                                       paste(fixed.effect, collapse = "+"),
                                       '+',
                                       paste('(1|',random.effect,')', collapse = "+")
    ))

    fmla2 <- stats::as.formula(object = paste("treatment ~ ",
                                       paste(fixed.effect, collapse = "+"),
                                       '+',
                                       paste('(1|',random.effect,')', collapse = "+")
    ))

    glmm1 = lme4::glmer(data = temp, formula = fmla1, family = binomial)

    glmm2 = lme4::glmer(data = temp, formula = fmla2, family = binomial)

    out = lmtest::lrtest(glmm1,glmm2)

    pval = out$`Pr(>Chisq)`[2]

    chisq = out$Chisq[2]

    final = c(pval = pval, chisq = chisq)

    return(final)
  }, mc.cores = ncores)

  names(out) = rownames(m)
  out = as.data.frame(do.call(rbind,out))

  out$fdr = p.adjust(out$pval, method = "fdr")
  out$feature = rownames(out)

  out = merge(out, features.meta, by='row.names')
  message("------- Done! -------")

  return(out)
}


