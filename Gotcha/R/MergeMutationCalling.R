#' To merge the outputs from BatchMutationCalling function into a single data frame
#'
#' @param out Path to the fastq or filtered fastq files
#' @return Saves a merged data frame containing all the outputs of the BatchMutationCalling function in the specified directory
#' @examples
#'

MergeMutationCalling = function(out){

  out = paste0(out,"Split/Filtered/")

  WhiteListMatch <- WTcount <- MUTcount <- WT <- MUT <- NULL # To prevent non-declared global variables

  message("------- LOADING MUTATION CALLING OUTPUTS -------")
  files = dir(path = out, pattern = "_rslurm")

  outs = list()
  for(i in files){
    outs[[i]] = readRDS(file = paste0(out,i,"/results_0.RDS"))
    outs[[i]] = outs[[i]][[1]][[1]][[1]]
  }

  message("------- MERGING MUTATION CALLING OUTPUTS -------")
  outs = do.call(rbind,outs)
  outs.collapse = outs %>% group_by(WhiteListMatch) %>% summarise(WTcount = sum(as.numeric(WTcount)),
                                                                  MUTcount = sum(as.numeric(MUTcount)))
  outs.collapse$MUTfraction = (outs.collapse$MUTcount / (outs.collapse$MUTcount + outs.collapse$WTcount))

  message("------- COLLAPSE BARCODE METRICS -------")
  message("------- Number of matched barcodes = ", nrow(outs.collapse))

  system(paste0('mkdir ',out,'MergedOuts'))
  save(outs.collapse, file = paste0(out,'MergedOuts/outs.collapsed.Rdata'))
  message("------- OUTPUT SAVED -------")
}
