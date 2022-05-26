#' Function to measure co-accessibility in subgroups of cells of interest
#'
#' @param archrProject ArchR project containing the output of the AddGenotypingArchr function
#' @param clusters Column name where the cluster labels are stored in the ArchR object metadata
#' @param select.clusters Clusters to be included
#' @param genotype Column name where the genotype labels are stored in the ArchR object metadata
#' @param genotype.levels Genotypes to be included in the analysis
#' @param peakSubset Set of peaks to include in the analysis
#' @param maxDist Maximum distance in base pairs to consider a pair of peaks
#' @param sample.column Column name indicating the sample identity
#' @param reducedDims Character vector of length one indicating the reduced dimensions to use
#' @param corCutOff Integer indicating the minimum correlation value to retain in the output
#' @param ncores Integer indicating how many cores to use for parallel computing
#' @return List class object containing the co-accessibility values and loops for each of the cell groups
#' @examples

DiffCoAccess = function(archrProject,
                        clusters,
                        select.clusters = NULL,
                        genotype,
                        genotype.levels,
                        peakSubset = NULL,
                        maxDist,
                        sample.column = NULL,
                        reducedDims = "IterativeLSI",
                        corCutOff = 0,
                        ncores = 1){

  if(class(archrProject) != "ArchRProject"){
    stop("Provided object is not ArchRProject class")
  }

  if(!(clusters %in% colnames(archrProject@cellColData))){
    stop("Provided cluster column name is not present in @cellColData")
  }

  if(!(genotype %in% colnames(archrProject@cellColData))){
    stop("Provided treatment column name is not present in @cellColData")
  }

  if(sum(genotype.levels %in% archrProject@cellColData[,genotype]) != 2){
    stop("Not all treatment levels are present in @cellColData treatment column")
  }

  message("------- Retrieving and binarizing cell by peak matrix... -------")
  pm = getMatrixFromProject(ArchRProj = archrProject, useMatrix = "PeakMatrix", binarize = T)

  message("------- Subsetting cell by peak matrix by cells with available genotype information... -------")

  if(!is.null(select.clusters)){
    pm = subset(pm,select = pm@colData[,genotype] %in% genotype.levels & pm@colData[,clusters] %in% select.clusters)
  }else{
    pm = subset(pm,select = pm@colData[,genotype] %in% genotype.levels)
  }

  peakSet <- getPeakSet(archrProject)

  if(!is.null(peakSubset)){

    if(class(peakSubset) == "SummarizedExperiment"){
      peakSubset = peakSubset[!is.na(peakSubset@assays@data$FDR$x),]

      peakSet <- peakSet[peakSet %over% makeGRangesFromDataFrame(peakSubset[peakSubset@assays@data$FDR$x < 0.05,]@elementMetadata),]
    }

    if(class(peakSubset) == "GRanges"){
      peakSet <- subsetByOverlaps(peakSet, peakSubset)
    }
  }
  message("------- ",length(peakSet@seqnames), " peaks and ",ncol(pm@assays@data$PeakMatrix)," cells retained -------")

  chr <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))

  peakSummits <- resize(peakSet, 1, "center")
  peakWindows <- resize(peakSummits, maxDist, "center")

  o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
  o <- o[o[,1] != o[,2], ]
  o$seqnames <- as.character(seqnames(peakSet))[o[,1]]
  o$idx1 <- peakSet$idx[o[,1]]
  o$idx2 <- peakSet$idx[o[,2]]

  o$Peak.1.class = NA
  o$Peak.2.class = NA
  o$Peak.1.distToTSS = NA
  o$Peak.2.distToTSS = NA
  o$group1.correlation = NA
  o$group2.correlation = NA
  o$Peak.1.nearestGene = NA
  o$Peak.2.nearestGene = NA
  o$correlation.delta = NA
  o$Odds.ratio = NA
  o$pvalue = NA
  o$detected.counts = NA

  # Calculate observed correlation within genotype
  message("----- Calculating correlation differences between genotypes for each cluster... -----")

  out = list()
  for(i in unique(pm@colData[,clusters])){

    message(paste0("------- Processing cluster ",i," -------"))

    message("------- Generating cell agregates and getting co-accessibility using ArchR... -------")

    archrProject = archrProject[!is.na(archrProject@cellColData[,genotype]) & archrProject@cellColData[,clusters] == i,]

    proj.1 = archrProject[as.vector(archrProject@cellColData[,clusters] == i & archrProject@cellColData[,genotype] %in% genotype.levels[1]),]
    proj.2 = archrProject[as.vector(archrProject@cellColData[,clusters] == i & archrProject@cellColData[,genotype] %in% genotype.levels[2]),]

    cells.1 = nrow(proj.1@cellColData)
    cells.2 = nrow(proj.2@cellColData)

    message("------- Getting co-accessibility for ",i," ", genotype.levels[1]," using ArchR -------")
    proj.1 = addCoAccessibility(
      ArchRProj = proj.1,
      maxDist = maxDist,
      k = round(sqrt(cells.1)),
      reducedDims = reducedDims)

    message("------- Getting co-accessibility for ",i," ", genotype.levels[2], " using ArchR -------")
    proj.2 = addCoAccessibility(
      ArchRProj = proj.2,
      maxDist = maxDist,
      k = round(sqrt(cells.2)),
      reducedDims = reducedDims)


    coA.1 <- getCoAccessibility(ArchRProj = proj.1,
                                corCutOff = 0,
                                resolution = 1,
                                returnLoops = F)


    coA.2 <- getCoAccessibility(ArchRProj = proj.2,
                                corCutOff = 0,
                                resolution = 1,
                                returnLoops = F)

    coA.1.loops <- getCoAccessibility(ArchRProj = proj.1,
                                      corCutOff = corCutOff,
                                      resolution = 1,
                                      returnLoops = T)

    coA.2.loops <- getCoAccessibility(ArchRProj = proj.2,
                                      corCutOff = corCutOff,
                                      resolution = 1,
                                      returnLoops = T)

    rownames(coA.1) = paste0(as.character(coA.1$seqnames),"_",coA.1$queryHits,"_",coA.1$subjectHits)
    rownames(coA.2) = paste0(as.character(coA.2$seqnames),"_",coA.2$queryHits,"_",coA.2$subjectHits)
    common = intersect(rownames(coA.1),rownames(coA.2))
    coA.1 = coA.1[common,]
    coA.2 = coA.2[common,]

    delta = coA.2$correlation - coA.1$correlation
    names(delta) = common

    message("------- Delta correlations calcualated! -------")

    message("------- Performing statistical test for co-accessibility differences -------")

    pm.temp = pm

    for(x in chr){

      message("Subsetting by ",x)

      tempWT = subset(pm.temp,
                      subset = as.character(pm@rowRanges@seqnames) == x,
                      select = pm.temp@colData[,clusters] == i & pm.temp@colData[,genotype] == genotype.levels[1])

      tempMUT = subset(pm.temp,
                       subset = as.character(pm@rowRanges@seqnames) == x,
                       select = pm.temp@colData[,clusters] == i & pm.temp@colData[,genotype] == genotype.levels[2])

      idx <- BiocGenerics::which(o$seqnames == x)

      rowCorCpp <- function(idxX, idxY, X, Y) {
        .Call('_ArchR_rowCorCpp', PACKAGE = 'ArchR', idxX, idxY, X, Y)}

      corVals_WT <- rowCorCpp(idxX = o[idx, ]$idx1, idxY = o[idx, ]$idx2, X = as.matrix(tempWT@assays@data$PeakMatrix), Y = as.matrix(tempWT@assays@data$PeakMatrix))
      corVals_MUT <- rowCorCpp(idxX = o[idx, ]$idx1, idxY = o[idx, ]$idx2, X = as.matrix(tempMUT@assays@data$PeakMatrix), Y = as.matrix(tempMUT@assays@data$PeakMatrix))

      obs.delta = corVals_MUT - corVals_WT

      o[idx, ]$group1.correlation <- as.numeric(corVals_WT)
      o[idx, ]$group2.correlation <- as.numeric(corVals_MUT)
      o[idx, ]$correlation.delta <- as.numeric(obs.delta)

      message("------- Performing Fisher test for connections for each peak in ",x," -------")

      peaks = unique(c(coA.1$queryHits, coA.2$queryHits))
      peaks = peaks[peaks %in% findOverlaps(getPeakSet(archrProject),peakSubset)@from]

      fisher.out = pbmclapply(peaks, function(z){

        group1 = table(factor(coA.1[coA.1$seqnames == x & coA.1$queryHits == z,]$FDR < 0.25 & coA.1[coA.1$seqnames == x & coA.1$queryHits == z,]$correlation > corCutOff, levels = c("TRUE","FALSE")))
        group2 = table(factor(coA.2[coA.2$seqnames == x & coA.2$queryHits == z,]$FDR < 0.25 & coA.2[coA.2$seqnames == x & coA.2$queryHits == z,]$correlation > corCutOff, levels = c("TRUE","FALSE")))

        df = cbind(group2, group1)

        ft = fisher.test(cbind(group2, group1))

        return(c(pval = ft$p.value, odds.ratio = ft$estimate, detected.counts = sum(df)))
      }, mc.cores = ncores)

      if(length(fisher.out)>0){

        fisher.out = as.data.frame(do.call(rbind,fisher.out))
        colnames(fisher.out) = c("pvalue","odds.ratio","detected.counts")
        fisher.out$chr = x
        fisher.out$peakIndex = peaks
        fisher.out$FDR = p.adjust(fisher.out$pvalue, method = "fdr")

      }
    }
    rownames(o) = paste0(o$seqnames,"_",o$queryHits,"_",o$subjectHits)

    cell.cov.1 = rowSums(pm.temp@assays@data$PeakMatrix[,pm.temp@colData[,genotype] == genotype.levels[1]])
    cell.cov.2 = rowSums(pm.temp@assays@data$PeakMatrix[,pm.temp@colData[,genotype] == genotype.levels[2]])
    names(cell.cov.1) = names(cell.cov.2) = peakSet$idx

    o$cell.coverage.1 = cell.cov.1[as.character(o$idx1)] + cell.cov.1[as.character(o$idx2)]
    o$cell.coverage.2 = cell.cov.2[as.character(o$idx1)] + cell.cov.2[as.character(o$idx2)]
    o$cell.coverage.total = o$cell.coverage.1+o$cell.coverage.2

    peakTypes = peakSet$peakType
    names(peakTypes) = peakSet$idx

    o$Peak.1.class = peakTypes[as.character(o$idx1)]
    o$Peak.2.class = peakTypes[as.character(o$idx2)]

    distances = peakSet$distToTSS
    names(distances) = peakSet$idx

    o$Peak.1.distToTSS = distances[as.character(o$idx1)]
    o$Peak.2.distToTSS = distances[as.character(o$idx2)]

    nearest = peakSet$nearestGene
    names(nearest) = peakSet$idx

    o$Peak.1.nearestGene = nearest[as.character(o$idx1)]
    o$Peak.2.nearestGene = nearest[as.character(o$idx2)]

    o$FDR = p.adjust(o$pvalue, method = "fdr")

    o$correlation.delta.agg = delta[rownames(o)]

    out[[i]] = list(diffCoAccess=o,diffConnections = fisher.out, coAccessGroup1 = coA.1.loops, coAccessGroup2 = coA.2.loops, coAccess1.no.loop = coA.1, coAccess2.no.loop = coA.2)
  }

  message("------- Done! -------")

  return(out)

}
