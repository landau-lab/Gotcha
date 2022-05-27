#' Function to split the fastq file into multiple chunks
#'
#' @param path Character string of lenght one corresponding to the full path to the folder containing the fastq files from the GoTChA library
#' @param out Character string of lenght one corresponding to the full path to the folder where to store the resulting fastq files after splitting.
#' @param reads Integer indicating the minimum number of bases below the quality threshold required to remove the read from downstream analysis
#' @param reads Which read to look into for quality filtering. It can be a pattern searc (i.e., R1|R2 or R1|R2|R3)
#' @param ncores Integer indicating the number of cores to use for parallel processing
#' @return Archr Project with added genotyping columns into the metadata
#' @examples
#'

FastqSplit = function(path,
                      out,
                      reads = 2000000,
                      ncores = 1){


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
    system(paste0("zcat ",path,"/",x," | split - -l ",format(reads*4,scientific=F)," --filter='gzip -f - > $FILE.gz' ",out,"Split/Split_ --additional-suffix _", gsub(x, pattern = ".fastq.gz", replacement = ""), ".fastq"))
  },mc.cores = ncores)

  message("------- DONE! -------")
}
