MutationCalling = function(out = "/path_to_filtered_fastqs/",
                           whitelist.file.path = "/path_to_whitelist/whitelist.txt",
                           wt.max.mismatch = 0,
                           mut.max.mismatch = 0,
                           keep.raw.reads = F,
                           ncores = 1,
                           reverse.complement = T,
                           testing = F,
                           which.read = "R1",
                           atac.barcodes = T,
                           primer.sequence = "CCTCATCATCCTCCTTGTC",
                           primed.max.mismatch = 3,
                           atac.barcodes.file.path = "/path_to_singlecell.csv",
                           wt.sequence =  "CGG",
                           mut.sequence = "CAG",
                           mutation.start = 31,
                           mutation.end = 34
){

  WhiteListMatch <- WTcount <- MUTcount <- WT <- MUT <- NULL # To prevent non-declared global variables
  message("------- BEGIN MUTATION CALLING FUNCTION -------")

  if(!file.exists(whitelist.file.path)){
    stop(paste0("---> Whitelist is not available in ",whitelist.file.path," <---"))
  }
  fastq.files = as.list(dir(path = out, pattern = ".fastq.gz"))
  if(length(fastq.files) == 0){
    stop("---> No fastq files detected in folder <---")
  }

  fastq.files.sequences = mclapply(fastq.files, function(x){
    temp = readFastq(dirPath = out, pattern = x)
    as.data.frame(sread(temp))
  }, mc.cores = ncores)
  names(fastq.files.sequences) = unlist(fastq.files)

  message("------- FASTQ FILES LOADED -------")

  # Subset fastq files for testing
  if(testing == T){
    fastq.files.sequences = mclapply(fastq.files.sequences, function(x){
      if(dim(x)[1] > 1000){
        x = as.data.frame(x[1:1000,])
      }else{
        x = as.data.frame(x)
      }
    }, mc.cores = ncores)
    names(fastq.files.sequences) = unlist(fastq.files)
    message("------- FASTQ FILES SUBSAMPLED TO 1K READS FOR TESTING -------")
  }

  # Extract sequences from fastq files
  barcodes= fastq.files.sequences[grep(names(fastq.files.sequences), pattern = "_R2_")]
  R1.sequence = fastq.files.sequences[grep(names(fastq.files.sequences), pattern = "_R1_")]
  R2.sequence = fastq.files.sequences[grep(names(fastq.files.sequences), pattern = "_R3_")]
  message("------- SEQUENCES OBTAINED FROM FASTQ FILES -------")

  # Read in data and generate sample data frame
  sample.S.index = as.list(unique(mclapply(strsplit(names(fastq.files.sequences), "_"), function(x) paste0(x[1:(length(x)-3)], collapse = "_"), mc.cores = ncores)))
  Output = mclapply(sample.S.index, function(x){
    data.frame(CB = barcodes[grep(names(barcodes), pattern = x)][[1]]$x,
               R1 = R1.sequence[grep(names(barcodes), pattern = x)][[1]]$x,
               R2 = R2.sequence[grep(names(barcodes), pattern = x)][[1]]$x)
  }, mc.cores = ncores)
  message("Data merged for each lane...")

  names(Output) = lapply(unique(lapply(strsplit(names(fastq.files.sequences), "_"), function(x) paste0(x[1:(length(x)-3)]))), function(y){
    y = paste(y, collapse = "_")
  })

  sample.Name.index = names(Output)

  OutputBind = lapply(sample.Name.index, function(x){
    if(length(Output[grep(names(Output), pattern = x, value = T)])>1){
      do.call(rbind,Output[grep(names(Output), pattern = x, value = T)])}else{Output[[grep(names(Output), pattern = x, value = T)]]}
  })
  names(OutputBind) = sample.Name.index

  message("------- DATAFRAME FOR EACH SAMPLE GENERATED -------")

  # Identify cell barcodes with primed reads and subset
  primed.index = mclapply(OutputBind, function(x){
    temp = vcountPattern(primer.sequence, substr(as.character(x[,which.read]), 0,nchar(primer.sequence)), max.mismatch = primed.max.mismatch, with.indels = FALSE)
    temp = temp == 1
    return(temp)
  }, mc.cores = ncores)
  names(primed.index) = names(OutputBind)

  Primed.output = lapply(as.list(names(OutputBind)), function(x) OutputBind[[x]][primed.index[[x]],])
  names(Primed.output) = names(primed.index)
  message("------- PRIMED READS IDENTIFIED -------")
  message("number of starting reads = ", nrow(OutputBind[[1]]))
  message("number of primed reads = ", nrow(Primed.output[[1]]))
  message("% of primed reads = ", round((nrow(Primed.output[[1]])/nrow(OutputBind[[1]])*100),2))

  # Read whitelist (either 10x whitelist or barcode list from CellRanger ATAC)
  if(atac.barcodes == T){
    whitelist = read.csv(atac.barcodes.file.path)
    whitelist$barcode = substr(whitelist$barcode, 1,16)
    whitelist = as.character(whitelist[-1,1])
    message("------- 10X ATAC BARCODES IDENTIFIED -------")
  }else{
    whitelist = as.character(readLines(whitelist.file.path))
    message("------- 10X WHITELIST IDENTIFIED -------")
  }

  # Match primed barcodes to whitelist
  if(reverse.complement){
    Primed.output = mclapply(Primed.output, function(x){
      x$CB = as.character(reverseComplement(DNAStringSet(x$CB)))
      return(x)
    },mc.cores = ncores)
    names(Primed.output) = names(primed.index)
    message("------- WARNING: Cell barcodes (CB) have been reverse complemented due to sequencer used -------")
  }

  message("------- BARCODE MATCHING BEGINNING -------")
  Matched.output = lapply(Primed.output, function(y){
    y$WhiteListMatch = unlist(mclapply(y$CB, function(x){
      temp = stringdist(x, whitelist, method = "hamming")
      if(min(temp) > 2){
        temp = "No Match"
      }else{
        if(sum(temp == 0) == 1) {
          temp = whitelist[temp==0]
        }else{
          if(sum(temp == 1) == 1){
            temp = whitelist[temp==1]
          }else{
            if(sum(temp == 2) == 1){
              temp = whitelist[temp==2]
            }else{
              temp = "Too many matches"
            }
          }
        }
      }
      return(temp)
    }, mc.cores = ncores))
    return(y)
  })

  pre_match_barcodes = nrow(Primed.output[[1]])
  no_match_barcodes = nrow(Matched.output[[1]][Matched.output[[1]]$WhiteListMatch == "No Match",])
  too_many_match_barcodes = nrow(Matched.output[[1]][Matched.output[[1]]$WhiteListMatch == "Too many matches",])
  Matched.output = lapply(Matched.output, function(x) x[!(x$WhiteListMatch %in% c("No Match","Too many matches")),])
  names(Matched.output) = names(primed.index)
  post_match_barcodes = nrow(Matched.output[[1]])


  message("------- BARCODE MATCHING DONE -------")
  message("total number of primed barcodes = ", pre_match_barcodes)
  message("total number of primed matched barcodes = ", post_match_barcodes)
  message("% barcode matching = ", round((post_match_barcodes/pre_match_barcodes)*100,2))

  # Search for genotyping information

  Genotyped.output = mclapply(Matched.output, function(x){
    x$WT = vcountPattern(wt.sequence, substr(as.character(x[,which.read]),start = mutation.start, stop = mutation.end), max.mismatch = wt.max.mismatch, with.indels = FALSE)
    x$MUT = vcountPattern(mut.sequence, substr(as.character(x[,which.read]),start = mutation.start, stop = mutation.end), max.mismatch = wt.max.mismatch, with.indels = FALSE)
    x$ReadGenotype = NA
    x$ReadGenotype[rowSums(cbind(x$WT == 1,x$MUT == 1)) == 2] = "Ambiguous"
    x$ReadGenotype[rowSums(cbind(x$WT == 0,x$MUT == 0)) == 2] = "No information"
    x$ReadGenotype[rowSums(cbind(x$WT == 1,x$MUT == 0)) == 2] = "WT"
    x$ReadGenotype[rowSums(cbind(x$WT == 0,x$MUT == 1)) == 2] = "MUTANT"
    return(x)
  }, mc.cores = ncores)
  names(Genotyped.output) = names(primed.index)

  # Remove unwanted reads
  Genotyped.output = lapply(Genotyped.output, function(x){
    x = x[!is.na(x$ReadGenotype),]
    return(x)
  })
  message("Reads genotyping done...")

  # Cell genotyping
  GoTChA_out_cell = lapply(Genotyped.output, function(x){
    x = x[,c("WhiteListMatch","WT","MUT","ReadGenotype")]
    x = x %>% group_by(WhiteListMatch) %>% summarise(WTcount = sum(as.numeric(WT)),
                                                     MUTcount = sum(as.numeric(MUT)))
    x = as.data.frame(x)
    return(x)
  })

  if(keep.raw.reads == T){
    output = list(genotyped.barcodes = GoTChA_out_cell,
                  raw.reads = Primed.output)
  }else{
    output = list(genotyped.barcodes = GoTChA_out_cell)
  }
  message("DONE!")
  return(output)
}