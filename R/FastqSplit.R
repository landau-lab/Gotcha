FastqSplit = function(path, out, reads = 2000000, ncores = 1){


  message("------- CREATING FILE INDEX -------")

  file.index = dir(path, pattern = ".fastq.gz")

  if(length(file.index)==0){
    stop("No fastq files detected in path folder")
  }

  message("------- CREATING OUTPUT FOLDER -------")

  if(!file.exists(paste0(out,"Split"))){
    out.create = paste0('mkdir ',out,"Split")
    system(command = out.create)
  }

  message("------- SPLITTING FASTQS -------")

  mclapply(file.index, function(x){
    system(paste0("gunzip -c ",path,"/",x," | split - -l ",format(reads*4,scientific=F)," --filter='gzip -f - > $FILE.gz' ",out,"Split/Split_ --additional-suffix _", gsub(x, pattern = ".fastq.gz", replacement = ""), ".fastq"))
  },mc.cores = ncores)

  message("------- DONE! -------")
}
