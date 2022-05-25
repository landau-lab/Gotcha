DiffLMM = function(archrProject,
                   useMatrix = "MotifMatrix",
                   provided.matrix = NULL,
                   subset = NULL,
                   sample.column,
                   cluster.column,
                   selected.clusters,
                   treatment.column,
                   treatment.levels,
                   ncores = 1){

  if(class(archrProject) != "ArchRProject"){
    stop("Specified object is not an ArchrProject")
  }

  if(is.null(provided.matrix)){
    if(!useMatrix %in% getAvailableMatrices(archrProject)){
      stop("Specified matrix not found in ArchRProject")
    }
  }

  if(!sample.column %in% colnames(archrProject@cellColData)){
    stop("Specified sample column not found in @cellColData slot")
  }

  if(!cluster.column %in% colnames(archrProject@cellColData)){
    stop("Specified pseudotime column not found in @cellColData slot")
  }

  if(!(treatment.column %in% colnames(archrProject@cellColData))){
    stop("Specified treatment column not found in @cellColData slot.")
  }

  if(!sum(treatment.levels %in% unique(archrProject@cellColData[,treatment.column])) == 2){
    stop("Specified treatment levels not found in specified treatment column in @cellColData slot")
  }

  if(sum(selected.clusters %in% unique(archrProject@cellColData[,cluster.column])) != length(selected.clusters)){
    message(paste0(setdiff(selected.clusters, unique(archrProject@cellColData[,cluster.column])), " not found in specified cluster column. Proceeding with: ",intersect(selected.clusters, unique(archrProject@cellColData[,cluster.column]))))
    if(!sum(selected.clusters %in% unique(archrProject@cellColData[,cluster.column])) == 0)
      stop("None of the selected clusters found in specified cluster.column")
  }

  if(is.null(provided.matrix)){
    mm <- getMatrixFromProject(archrProject, useMatrix = useMatrix)
    if(useMatrix == "MotifMatrix"){
      m <- mm@assays@data$z # Use z score
    }
    if(useMatrix == "GeneScoreMatrix"){
      m <- mm@assays@data$GeneScoreMatrix # Use gene scores
      rownames(m) = mm@elementMetadata$name # Add gene names to matrix
    }
  }else{
    m <- provided.matrix
  }

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
    temp = as.data.frame(archrProject@cellColData)
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
