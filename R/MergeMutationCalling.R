MergeMutationCalling = function(out){

  out = paste0(out,"Split/Filtered/")

  matched_barcodes <- WTcount <- MUTcount <- WT <- MUT <- NULL # To prevent non-declared global variables

  message("------- LOADING MUTATION CALLING OUTPUTS -------")
  files = dir(path = out, pattern = "mutation_call")

  outs = list()
  for(i in files){
    outs[[i]] = readRDS(file = paste0(out,i))
    outs[[i]] <- as.data.frame(outs[[i]])
  }

  message("------- MERGING MUTATION CALLING OUTPUTS -------")
  outs = do.call(rbind,outs)
  colnames(outs) <- gsub("genotyped_reads.x.WT", "WTcount", colnames(outs))
  colnames(outs) <- gsub("genotyped_reads.x.MUT", "MUTcount", colnames(outs))
  outs.collapse = outs %>% group_by(matched_barcodes) %>% summarise(WTcount = sum(as.numeric(WTcount)),
                                                                  MUTcount = sum(as.numeric(MUTcount)))
  outs.collapse$MUTfraction = (outs.collapse$MUTcount / (outs.collapse$MUTcount + outs.collapse$WTcount))

  # remove "Too many match" and "No match" barcodes
  outs.collapse <- outs.collapse[-which(outs.collapse$matched_barcodes == "Too many matches"),]
  outs.collapse <- outs.collapse[-which(outs.collapse$matched_barcodes == "No match"),]

  message("------- COLLAPSE BARCODE METRICS -------")
  message("------- Number of matched barcodes = ", nrow(outs.collapse))

  if (!file.exists(paste0(out,'MergedOuts'))){
    system(paste0('mkdir ',out,'MergedOuts'))
  } else{
    message("WARNING: 'MergedOuts' FOLDER EXISTS AND WILL BE REWRITTEN")
  }
  saveRDS(outs.collapse, file = paste0(out,'MergedOuts/outs.collapsed.Rdata'))
  message("------- OUTPUT SAVED -------")
}
