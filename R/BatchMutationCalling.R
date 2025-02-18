#' Define read genotype and read counts per genotype for each cell barcode using parallel submission of mutation_calling function
#'
#' @param out Path to the fastq or filtered fastq files
#' @param wt.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the wild-type sequence
#' @param mut.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the mutant sequence
#' @param keep.raw.reads Logical. Whether to return the raw reads in the output file. Defaults to false
#' @param ncores Integer indicating the number of cores to use for parallel processing
#' @param reverse.complement Whether to take the reverse complement of the cell barcodes
#' @param testing Logical indicating whether to sample the first 1,000 reads for testing the function
#' @param which.read Which read to select to look for the mutation site
#' @param barcodes.file.path Path to the file containing the cell barcodes detected in the experiment
#' @param wt.sequence Character vector of length one specifying the expected wild-type sequence
#' @param mut.sequence Character vector of length one specifying the expected mutant sequence
#' @param mutation.start Position in which the expected wild-type or mutant sequence starts in the read
#' @param mutation.end Position in which the expected wild-type or mutant sequence ends in the read
#' @param max.distance Integer indicating the number of mismatches allowed between barcodes and barcode whitelist
#' @param rewrite Logical indicating whether the chunk should be rewritten if it already exists
#' @return NULL


BatchMutationCalling = function(out = "/path_to_filtered_fastqs/",
                                wt.max.mismatch = 0,
                                mut.max.mismatch = 0,
                                keep.raw.reads = F,
                                ncores = 1,
                                reverse.complement = T,
                                testing = F,
                                which.read = "R1",
                                barcodes.file.path = "/path_to_singlecell.csv",
                                wt.sequence =  "CGG",
                                mut.sequence = "CAG",
                                mutation.start = 31,
                                mutation.end = 34,
                                max.distance = 2,
                                rewrite = TRUE
){
  if (max.distance > 2){ message("WARNING: MAX DISTANCE ", max.distance, " ALLOWS MANY MISMATCHES!")}
  if (testing == TRUE){ message("------- READS WILL BE SUBSETTED FOR TESTING -------")}
  if(!file.exists(barcodes.file.path)) stop("BARCODE FILE DOES NOT EXIST")

  options(expressions = 2.5e5) # Increase the number of nested expressions to be evaluated. Limit is 5e5.

  WhiteListMatch <- WTcount <- MUTcount <- WT <- MUT <- NULL # To prevent non-declared global variables

  out = paste0(out,"Split/Filtered/")

  message("------- GENERATING CHUNK INDEX -------")
  chunk.index = unique(unlist(lapply(dir(out), function(x){
    lapply(strsplit(x,"_"), function(y) y[length(y)])
  })))

  # only select the chunk.indexes from the directory
  chunk.index = chunk.index[grepl("^[a-z]{2}$", chunk.index)]

  message("------- GENERATING PARAMETERS -------")
  pars = list(out = out,
              wt.max.mismatch = wt.max.mismatch,
              mut.max.mismatch = mut.max.mismatch,
              ncores = ncores,
              reverse.complement = reverse.complement,
              testing = testing,
              which.read = which.read,
              barcodes.file.path = barcodes.file.path,
              wt.sequence =  wt.sequence,
              mut.sequence = mut.sequence,
              mutation.start = mutation.start,
              mutation.end = mutation.end,
              max.distance = max.distance)

  pars = as.data.frame(pars,stringsAsFactors = F)

  if (!dir.exists(out)) {message("ERROR: DIRECTORY ",out ," DOES NOT EXIST")}
  setwd(out)

  message("------- CALLING MUTATIONS -------")

  mclapply(chunk.index, function(x){
    pars$out = paste0(out,x,"/")
    # run each chunk with only one core
    pars$ncores = 1
    parameters <- as.list(pars[1,])
    tryCatch({
      if (!file.exists(paste0(out, "mutation_call_", x, ".rds")) & rewrite == FALSE) {
        mutation_calling <- do.call(MutationCalling, parameters)
        saveRDS(object = mutation_calling, file = paste0(out, "mutation_call_", x, ".rds"))
      } else if (rewrite == TRUE) {
        mutation_calling <- do.call(MutationCalling, parameters)
        saveRDS(object = mutation_calling, file = paste0(out, "mutation_call_", x, ".rds"))
      } else {
        message("CHUNK ", x, " WILL NOT BE REWRITTEN. SKIPPING...")
      }
    }, error = function(e) {
      message("ERROR IN CHUNK ", x, "!!: ", conditionMessage(e), "\n")
    })
  }, mc.cores = ncores)

  message("------- FINISHED MUTATION CALLING! -------")
  return()
}


