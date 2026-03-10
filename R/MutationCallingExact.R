
#' helper function to process fastqs
#' @param out Path to the fastq files
#' @param pattern suffix to search for files in the directory
#' @param ncores number of cores to use with parallel processing
#' @return Datatable with information from fastq files
#'
read_and_process_fastq <- function(path, pattern, ncores) {
  fastq_files <- dir(path, pattern = pattern, full.names = TRUE)
  if (length(fastq_files) == 0) stop("No fastq files detected.")
  res = mclapply(fastq_files, function(file) {
    temp <- readFastq(file)
    as.data.table(sread(temp))
  }, mc.cores = ncores)
  names(res) = fastq_files
  return(res)
}


#' helper function to subset fastqs
#' @param fastq_data Datatable with fastq data
#' @param max_reads number of reads to subset
#' @param ncores number of cores to use with parallel processing
#' @return Datatable with subsetted fastq data
#'
subset_for_testing <- function(fastq_data, max_reads = 10000, ncores) {
  mclapply(fastq_data, function(dt) {
    if (nrow(dt) > max_reads) dt[1:max_reads, ] else dt
  }, mc.cores = ncores)
}


#' helper function to convert barcodes to numerical matrix for RANN clustering
#' @param input Datatable with barcodes for each read
#' @return Datatable with barcodes converted to integers
#'
convert_to_numeric_matrix <- function(input) {
  # Convert each character in the barcode to a numeric value
  char_to_num <- c(A = 1, C = 2, G = 3, T = 4)
  numeric_matrix <- t(sapply(input, function(x) {
         char_to_num[unlist(strsplit(x, ""))]
     }))
  dimnames(numeric_matrix) <- NULL
  # convert na to 5
  numeric_matrix[is.na(numeric_matrix)] <- 5
  numeric_matrix <- as.data.table(numeric_matrix)
  return(numeric_matrix)
}


#' same as previous function but barcode letters converted to different integers
#' @param input Datatable with barcodes for each read
#' @return Datatable with barcodes converted to integers
#'
convert_to_numeric_matrix_reorder <- function(input) {
  char_to_num <- c(A = 1, T = 2, C = 3, G = 4)
  numeric_matrix <- t(sapply(input, function(x) {
    char_to_num[unlist(strsplit(x, ""))]
  }))
  dimnames(numeric_matrix) <- NULL
  numeric_matrix[is.na(numeric_matrix)] <- 5
  numeric_matrix <- as.data.table(numeric_matrix)
  return(numeric_matrix)
}


#' read in file with barcode whitelist
#' @param filepath path of barcode whitlist
#' @return list of whitelist barcodes
#'
load_whitelist <- function(file_path) {
  if (!file.exists(file_path)) stop("Whitelist file not found: ", file_path)
  whitelist_data <- fread(file_path)
  if ("barcode" %in% colnames(whitelist_data)) {
    whitelist <- whitelist_data$barcode
  } else {
    # Assume a single-column CSV without a header
    whitelist <- whitelist_data[[1]]
  }
  # Trim to consistent length if needed (e.g., 16 characters for 10x barcodes)
  whitelist <- substr(whitelist, 1, 16)
  return(unique(whitelist))
}


#' perform RANN clustering to find which whitelist barcodes are closest to gotcha barcodes
#' @param barcodes Datatable with gotcha barcodes
#' @param whitelist Datatable with whitelist barcodes
#' @param order Integer specifying which function to convert barcodes to numeric matrix
#' @param radius radius with nearest neighbor RANN clustering
#' @param nearest_neighbors number of nearest neighbors to use in RANN clustering
#' @return Matrix with closest barcodes from whitelist
#'
rann_matching <- function(barcodes, whitelist, order, radius, nearest_neighbors) {
  # Convert barcodes and whitelist to numeric matrices
  if (order ==1){
    barcode_matrix <- convert_to_numeric_matrix(barcodes)
    whitelist_matrix <- convert_to_numeric_matrix(whitelist)
  } else{
    barcode_matrix <- convert_to_numeric_matrix_reorder(barcodes)
    whitelist_matrix <- convert_to_numeric_matrix_reorder(whitelist)
  }
  # Use RANN for nearest neighbor search
  nn_results <- nn2(
    data = whitelist_matrix,    # Whitelist data
    query = barcode_matrix,     # Query barcodes
    k = nearest_neighbors,      # Find num of nearest neighbor
    searchtype = "radius",      # Use radius-based search
    radius = radius       # Maximum allowed distance, will only allow for 2 mismatches
  )
  whitelist_ind <- nn_results$nn.idx
  # convert 0 to NA
  whitelist_ind <- replace(whitelist_ind, whitelist_ind==0, NA)
  # return the nearest whitelist barcodes
  mapped_matrix <- matrix(whitelist[whitelist_ind], nrow = nrow(nn_results$nn.idx), ncol = ncol(nn_results$nn.idx))
  return(mapped_matrix)
}

#' perform hamming match. check if barcode has more than one whitelist entry with one mismatch
#' @param barcodes gotcha barcode
#' @param whitelist list of whitelist barcodes
#' @return the original barcode character or "Too many matches" if there are mult match with hamming dis=1
#'
hamming_too_many_match <- function(barcodes, whitelist) {
  # get number of mismatches between strings
  dist_mat <- stringdistmatrix(barcodes, unique(whitelist), method = "hamming")
  if (length(which(dist_mat==1)) > 1){
    return("Too many matches")
  }
  return(barcodes)
}

#' perform genotyping of reads
#' @param reads gotcha barcode
#' @param wt_seq Character vector of length one specifying the expected wild-type sequence
#' @param mut_seq Character vector of length one specifying the expected mutant sequence
#' @param mutation_start Position in which the expected wild-type or mutant sequence starts in the read
#' @param mutation_end Position in which the expected wild-type or mutant sequence ends in the read
#' @param wt_max_mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the wild-type sequence
#' @param mut_max_mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the mutant sequence
#' @param ncores number of cores to use for parallel processing
#' @return Datatable with read and genotype call
#'
genotype_reads <- function(reads, wt_seq, mut_seq, mutation_start, mutation_end, wt_max_mismatch, mut_max_mismatch, ncores) {
  mclapply(reads, function(read) {
    wt_count <- vcountPattern(wt_seq, substr(read, mutation_start, mutation_end), max.mismatch = wt_max_mismatch)
    mut_count <- vcountPattern(mut_seq, substr(read, mutation_start, mutation_end), max.mismatch = mut_max_mismatch)
    genotype <- ifelse(wt_count == 1 & mut_count == 1, "Ambiguous",
                       ifelse(wt_count == 0 & mut_count == 0, "No information",
                              ifelse(wt_count == 1 & mut_count == 0, "WT", "MUTANT")))
    data.table(WT = wt_count, MUT = mut_count, Genotype = genotype)
  }, mc.cores = ncores)
}


#' Define read genotype and read counts per genotype for each cell barcode
#' @param out Path to the fastq or filtered fastq files
#' @param barcodes.file.path Path to the file containing the cell barcodes detected in the experiment
#' @param wt.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the wild-type sequence
#' @param mut.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the mutant sequence
#' @param ncores Integer indicating the number of cores to use for parallel processing
#' @param reverse.complement Whether to take the reverse complement of the cell barcodes
#' @param testing Logical indicating whether to sample the first 10,000 reads for testing the function
#' @param which.read Which read to select to look for the mutation site
#' @param wt.sequence Character vector of length one specifying the expected wild-type sequence
#' @param mut.sequence Character vector of length one specifying the expected mutant sequence
#' @param mutation.start Position in which the expected wild-type or mutant sequence starts in the read
#' @param mutation.end Position in which the expected wild-type or mutant sequence ends in the read
#' @param max.distance Maximum number of mismatches allowed between barcodes and whitelist
#' @return output Datatable with barcode and genotype calls
#'
MutationCallingExact <- function(out = "/path_to_filtered_fastqs/", barcodes.file.path = "/path_to_whitelist/whitelist.txt",
                            wt.max.mismatch = 0, mut.max.mismatch = 0, ncores = 1, reverse.complement = TRUE,
                            testing = FALSE, which.read = "R1", wt.sequence = "CGG", mut.sequence = "CAG",
                            mutation.start = 31, mutation.end = 34, max.distance = 2) {

  # make output file
  out_file <- paste0(out, "out.log")
  # get chunk name
  chunk_name <- basename(out)
  cat(paste0("------- BEGIN MUTATION CALLING ", chunk_name ," -------"), file=out_file,sep="\n")
  # Load whitelist from the specified file path
  whitelist <- load_whitelist(barcodes.file.path)
  # remove "NO_BARCODE" from whitelist
  whitelist <- whitelist[whitelist != "NO_BARCODE"]
  # Load FASTQ files
  fastq_data <- read_and_process_fastq(out, pattern = ".fastq.gz", ncores = ncores)
  if (testing) fastq_data <- subset_for_testing(fastq_data, max_reads = 10000, ncores = ncores)
  cat(paste0("------- FASTQ FILES LOADED ", chunk_name ," -------"), file=out_file, sep = "\n", append=TRUE)
  # Process sequences
  barcodes <- fastq_data[[grep(names(fastq_data), pattern = "_R2_")]]
  reads <- fastq_data[[grep(names(fastq_data), pattern = paste0("_", which.read, "_"))]]
  # reverse complement and convert to dnastringset
  if (reverse.complement) {
    barcodes <- lapply(barcodes, function(x) reverseComplement(DNAStringSet(x)))
    cat(paste0("------- ", chunk_name ," CELL BARCODES HAVE BEEN REVERSE COMPLEMENTED -------"), file=out_file, sep = "\n", append=TRUE)
  }
  cat(paste0("------- STARTING BARCODE MATCHING ", chunk_name ," -------"), file=out_file, sep = "\n", append=TRUE)
  # check which barcodes are perfect match
  matched_barcodes_ind <- (as.character(barcodes$x) %in% whitelist)
  matched_barcodes <- rep("No match", length(matched_barcodes_ind))
  # add the perfect match barcodes
  matched_barcodes[matched_barcodes_ind] <- as.character(barcodes$x)[matched_barcodes_ind]

  sink(out_file)
  cat(paste0("------- BARCODE MATCHING COMPLETED ", chunk_name ," -------"), sep = "\n", append=TRUE)
  cat(paste0("total number of ", chunk_name ," starting barcodes = ", length(matched_barcodes)), sep = "\n", append=TRUE)
  end_bc <- length(which(matched_barcodes != "No match" & matched_barcodes != "Too many matches"))
  cat(paste0("total number of ", chunk_name, " matched barcodes = ", end_bc), sep = "\n", append=TRUE)
  cat(paste0("% ", chunk_name, " barcode matching = ", round((end_bc/length(matched_barcodes))*100,2)), sep = "\n", append=TRUE)
  cat(paste0("------- STARTING PER READ GENOTYPING ", chunk_name ," -------"), sep = "\n", append=TRUE)
  sink()
  # Genotype reads
  genotyped_reads <- genotype_reads(reads, wt.sequence, mut.sequence, mutation.start, mutation.end, wt.max.mismatch, mut.max.mismatch, ncores)
  cat(paste0("------- GENOTYPING COMPLETED ", chunk_name ," -------"), file = out_file, sep = "\n", append=TRUE)
  cat(paste0("------- SAVING OUTPUT ", chunk_name ," ... -------"), file = out_file, sep = "\n", append=TRUE)
  # Output processing
  output <- list(matched_barcodes = matched_barcodes, genotyped_reads = genotyped_reads)

  cat(paste0("------- ", chunk_name ," CHUNK DONE! -------"), file = out_file, sep = "\n", append=TRUE)
  return(output)
}
