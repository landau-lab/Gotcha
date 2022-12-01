FilterGenotyping = function(archrProject, total.counts.column, min.reads = 100, plot = F){

  if(class(archrProject) != "ArchRProject"){
    stop("The specified object is not ArchrProject class")
  }

  if(!(total.counts.column %in% colnames(archrProject@cellColData))){
    stop("The specified total counts column is not found in the ArchRProject @cellColData slot")
  }

  message("------- GETTING TOTAL COUNTS FROM SPECIFIED COLUMN -------")

  dat = archrProject@cellColData[,total.counts.column]

  dat.subset = dat[dat >= min.reads]
  dat.subset = dat.subset[!is.na(dat.subset)]
  message(paste0("Initial minimum reads = ",min.reads))

  if(plot){
    hist(log10(dat+1), breaks=100, main = total.counts.column, xlab = "Total reads (log10)")
    abline(v = log10(min.reads), lty=2, col = "red")
  }

  archrProject@cellColData[,paste0(total.counts.column,"_above_read_th")] = archrProject@cellColData[,total.counts.column] >= min.reads

  target = NULL
  target = unlist(lapply(strsplit(total.counts.column, "_"), function(x) x[1])) # Retrieve target name
  archrProject@cellColData[,paste0(target,"_Genotype_above_read_th")] = archrProject@cellColData[,paste0(target,"_Genotype_noise_corrected")]
  archrProject@cellColData[,paste0(target,"_Genotype_above_read_th")][archrProject@cellColData[,paste0(total.counts.column,"_above_read_th")] ==FALSE] = NA

  message(paste0(total.counts.column,"_above_read_th and ", target,"_Genotype_above_read_th columns added to ArchR object @cellColData slot."))
  message("------- DONE! -------")

  return(archrProject)
}
