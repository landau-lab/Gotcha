#' Function to plot co-accessibility of subgroups of cells of interest
#'
#' @param archrProject ArchR project containing the output of the AddGenotypingArchr function
#' @param DCA Output obtained from the DiffCoAccess function
#' @param coords GRanges indicating the coordinates of the region to be included in the plot
#' @param features GRanges list class object containing features to be included in the track plot (i.e., motif sites, gene annotations, locus control regions)
#' @param condition.column Name of the column in the ArchR metadata containing a condition to subset the cells included
#' @param condition Levels of the condition column to be included in the analysis 
#' @param cluster.genotype.column Name of the column in the ArchR metadata containing a cluster and treatment labels
#' @param cluster.genotype.levels  Levels of the column in the ArchR metadata containing a cluster and treatment labels to be compared
#' @param min.cor Integer indicating the minimum correlation value to retain a connection in the plot
#' @param path Path to the folder where to store the plots as .pdf files
#' @param plot Logical, whether to visualize the plots in the R studio viewer
#' @return Prints .pdf files of the plots in the indicated path. If plot = T, it will also plot it in the R studio viewer
#' @examples

PlotDCA = function(archrProject,
                   DCA,
                   coords,
                   features,
                   condition.column = NULL,
                   condition = NULL,
                   cluster.genotype.column,
                   cluster.genotype.levels,
                   min.cor = 0.1,
                   path = NULL,
                   plot = TRUE){

  if(class(archrProject) != "ArchRProject"){
    stop("Provided object is not ArchRProject class")
  }

  if(class(DCA) != "list"){
    stop("DCA object must be a list (output of DiffCoAccessibility function)")
  }

  if(class(coords)[1] != "GRanges"){
    stop("coords object must be of GRanges class")
  }

  if(!(cluster.genotype.column %in% colnames(archrProject@cellColData))){
    stop("Provided cluster plus genotype column name is not present in @cellColData")
  }

  if(!(condition.column %in% colnames(archrProject@cellColData))){
    stop("Provided condition column name is not present in @cellColData")
  }

  if(sum(cluster.genotype.levels %in% archrProject@cellColData[,cluster.genotype.column]) != 2){
    stop("Not all cluster.genotype.levels are present in @cellColData cluster.genotype.column")
  }

  if(sum(condition %in% archrProject@cellColData[,condition.column]) == 0){
    stop("Specified condition not found in @cellColData condition.column")
  }

  if(!is.null(path)){
    if(!dir.exists(path)){
      stop("Specified folder in the path variable does not exist")
    }
  }

  peakSubset = subsetByOverlaps(getPeakSet(archrProject),coords)
  cor.diff = list()
  for(i in names(DCA)){

    loops.wt = subsetByOverlaps(DCA[[i]]$coAccessGroup1$CoAccessibility, peakSubset)
    loops.mut = subsetByOverlaps(DCA[[i]]$coAccessGroup2$CoAccessibility, peakSubset)

    idx1 = unlist(lapply(strsplit(rownames(DCA[[i]]$diffCoAccess),"_"), function(x) paste0(x[1],"_",x[2])))
    idx2 = unlist(lapply(strsplit(rownames(DCA[[i]]$diffCoAccess),"_"), function(x) paste0(x[1],"_",x[3])))
    unique.idx = unique(c(idx1,idx2))

    cor.diff[[i]] = list()

    for(j in 1:length(peakSubset)){

      message("------- Processing co-accessibility for ",peakSubset[j,]," -------")

      chr.ind = as.character(peakSubset@seqnames[j])
      unique.ind = peakSubset$idx[j]

      # For group1

      loops.index.right = GRanges(loops.wt@seqnames, ranges = IRanges(start = loops.wt@ranges@start, end = loops.wt@ranges@start+1))
      loops.index.left = GRanges(loops.wt@seqnames, ranges = IRanges(start = loops.wt@ranges@start + loops.wt@ranges@width, end = loops.wt@ranges@start + loops.wt@ranges@width + 1))

      indx1 = findOverlaps(loops.index.right, peakSubset[peakSubset@seqnames == chr.ind & peakSubset$idx == unique.ind,])
      indx2 = findOverlaps(loops.index.left, peakSubset[peakSubset@seqnames == chr.ind & peakSubset$idx == unique.ind,])

      loops.wt.temp = loops.wt[c(indx1@from,indx2@from),]

      # For group2

      loops.index.right = GRanges(loops.mut@seqnames, ranges = IRanges(start = loops.mut@ranges@start, end = loops.mut@ranges@start+1))
      loops.index.left = GRanges(loops.mut@seqnames, ranges = IRanges(start = loops.mut@ranges@start + loops.mut@ranges@width, end = loops.mut@ranges@start + loops.mut@ranges@width + 1))

      indx1 = findOverlaps(loops.index.right, peakSubset[peakSubset@seqnames == chr.ind & peakSubset$idx == unique.ind,])
      indx2 = findOverlaps(loops.index.left, peakSubset[peakSubset@seqnames == chr.ind & peakSubset$idx == unique.ind,])

      loops.mut.temp = loops.mut[c(indx1@from,indx2@from),]

      #Plotting
      if(plot){
        p = plotBrowserTrack(
          ArchRProj = archrProject[as.vector(archrProject@cellColData[,condition.column] == condition & archrProject@cellColData[,cluster.genotype.column] %in% c(cluster.genotype.levels)),],
          groupBy = cluster.genotype.column,
          region = coords,
          ylim = c(0,1),
          features = features,
          loops = GRangesList(WT = loops.wt.temp, MUT = loops.mut.temp),
          plotSummary = c("bulkTrack","featureTrack", "geneTrack","loopTrack"),
          title = paste0(i," - MUT")
        )

        if(is.null(path)){
          grid::grid.newpage()
          grid::grid.draw(p)
        }else{
          pdf(file = paste0(path,peakSubset[j,],"_coAccess.pdf"), width = 15, height = 10)
          grid::grid.draw(p)
          dev.off()
        }
      }

      # To store number of co-accessible pairs for this peak
      cor.diff[[i]][[j]] = data.frame(Peak = paste0(peakSubset[j,]),
                                      DiffGroup1 = sum(loops.wt.temp$correlation > min.cor),
                                      DiffGroup2 = sum(loops.mut.temp$correlation > min.cor),
                                      Delta = sum(loops.mut.temp$correlation > min.cor) - sum(loops.wt.temp$correlation > min.cor))
    }
  }
  for(i in names(cor.diff)){
    if(length(cor.diff[[i]])>1){
      cor.diff[[i]] = as.data.frame(do.call(rbind, cor.diff[[i]]))
      cor.diff[[i]] = cor.diff[[i]][order(cor.diff[[i]]$Delta, decreasing = T),]
      cor.diff[[i]]$Peak = factor(cor.diff[[i]]$Peak, levels =  cor.diff[[i]]$Peak)
      p2 = ggplot(cor.diff[[i]], aes(x = Peak, y = Delta))+
        geom_bar(stat = "identity") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(y = paste0("Delta in co-accessible pairs above ",min.cor," - (", cluster.genotype.levels[2]," - ",cluster.genotype.levels[1],")"), title = paste0(coords))
      print(p2)
    }
    if(length(cor.diff[[i]]) == 1){
      p2 = ggplot(cor.diff[[i]][[j]], aes(x = Peak, y = Delta))+
        geom_bar(stat = "identity") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(y = paste0("Difference in number of co-accessible pairs above ",min.cor," - (", cluster.genotype.levels[2]," - ",cluster.genotype.levels[1],")"), title = paste0(coords))
      print(p2)
    }
  }
  message("------- Done! -------")
}
