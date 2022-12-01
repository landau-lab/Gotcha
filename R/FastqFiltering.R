FastqFiltering = function(path, out="/path_to_fastqs/", 
						reads = 2000000, 
                        min.quality = 15,
                        min.bases = 1,
                        which.read = "R1",
                        read.region = NULL,
                        ncores = 1){

  message("------- STARTING FASTQSPLIT -------")
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
  
  message("------- BEGIN FASTQ FILTERING FUNCTION -------")

  out = paste0(out,"Split/")

  # get fastq file names
  fastq.files = as.list(dir(path = out, pattern = ".fastq.gz"))

  # error message if path to files is incorrect
  if(length(fastq.files) == 0){
    stop("---> No fastq files detected in folder <---")
  }

  # Define chunk groups based on split ID
  chunk.index = unlist(mclapply(fastq.files, function(x){
    unlist(lapply(strsplit(x, "_"), function(y) y[2]))
  }, mc.cores = ncores))

  # choose which fastq files to use for filtering based on which.read
  selected.files = fastq.files[grep(fastq.files, pattern = which.read)]
  selected.chunk.index = chunk.index[grep(fastq.files, pattern = which.read)]

  message("------- FASTQ FILES IDENTIFIED -------")

  # Generate index for each read in each fastq file - boolean for whether a read passes minimum quality / minimum bases
  index = mclapply(selected.files, function(x){
    temp = readFastq(dirPath = out, pattern = x)
    q = as(quality(temp),"matrix")
    if(!is.null(read.region) & min.bases <= length(read.region)){
      q = q[,read.region]
    }
    rowSums(q <= min.quality) < min.bases
  }, mc.cores = ncores)
  names(index) = selected.chunk.index

  message("------- PER READ INDEX CREATED -------")

  # read in unfiltered fastq files
  orig.files = mclapply(fastq.files, function(x){
    temp = readFastq(dirPath = out, pattern = x)
  }, mc.cores = ncores)
  names(orig.files) = unlist(fastq.files)

  message("------- RAW FASTQ FILES LOADED -------")

  # apply index to generate filtered fastq files
  filter.files = mclapply(names(orig.files), function(x){
    c.ind = unlist(lapply(strsplit(x, "_"), function(y) y[2]))
    orig.files[[x]][index[[c.ind]]]
  }, mc.cores = ncores)
  names(filter.files) = names(orig.files)
  message("------- FILTERING COMPLETE -------")

  system(paste0('mkdir ',out,"Filtered/"))

  mclapply(unique(chunk.index), function(x){
    system(paste0('mkdir ', out,"Filtered/",x))
  }, mc.cores = ncores)

  # write out filtered fastq files
  mclapply(unique(chunk.index), function(x) {
    selected = grep(names(filter.files), pattern = paste0("_", x,"_"), value = T)
    lapply(selected, function(y){
      writeFastq(object = filter.files[[y]], file = paste0(out,"Filtered/",x,"/filtered_bq",min.quality,"_bn",min.bases,"_",y), mode = "a")
    })
  }, mc.cores = ncores)

  message("------- FASTQ FILTERING METRICS -------")
  message("------- mean number of filtered reads = ", mean(unlist(lapply(filter.files, length))))
  message("------- % of remaining reads after filtering= ", round(mean(unlist(lapply(filter.files, length))/length(orig.files[[1]])*100),3)," +/- ", round(sd(unlist(lapply(filter.files, length))/length(orig.files[[1]])*100),3)," (mean +/- sd)")
  message("------- FASTQ FILTERING FUNCTION COMPLETE -------")
  
}
