#' Define read genotype and read counts per genotype for each cell barcode using parallel job submission to slurm cluster
#'
#' @param out Path to the fastq or filtered fastq files
#' @param whitelist.file.path Path to the file containing the whitelist with possible cell barcodes
#' @param wt.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the wild-type sequence
#' @param mut.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the mutant sequence
#' @param keep.raw.reads Logical. Whether to return the raw reads in the output file. Defaults to false
#' @param ncores Integer indicating the number of cores to use for parallel processing
#' @param reverse.complement Whether to take the reverse complement of the cell barcodes
#' @param testing Logical indicating whether to sample the first 1,000 reads for testing the function
#' @param which.read Which read to select to look for the mutation site
#' @param atac.barcodes Logical indicating whether to only use the detected atac barcodes in the experiment as whitelist
#' @param primer.sequence Character vector of length one indicating the primer sequence
#' @param primed.max.mismatch  Integer indicating the maximum number of mismatches accepted when searching for the primer sequence
#' @param atac.barcodes.file.path Path to the file containing the cell barcodes detected in the experiment
#' @param wt.sequence Character vector of length one specifying the expected wild-type sequence
#' @param mut.sequence Character vector of length one specifying the expected mutant sequence
#' @param mutation.start Position in which the expected wild-type or mutant sequence starts in the read
#' @param mutation.end Position in which the expected wild-type or mutant sequence ends in the read
#' @param soptions List class object specifying the parameters for the jobs submitted to the slurm cluster
#' @return Archr Project with added genotyping columns into the metadata
#' @examples
#'


BatchMutationCalling = function(out = "/path_to_filtered_fastqs/",
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
                                mutation.end = 34,
                                soptions = list(mem = '20g', 'cpus-per-task' = 10)
){

  WhiteListMatch <- WTcount <- MUTcount <- WT <- MUT <- NULL # To prevent non-declared global variables

  out = paste0(out,"Split/Filtered/")

  message("------- GENERATING CHUNK INDEX -------")
  chunk.index = unique(unlist(lapply(dir(out), function(x){
    lapply(strsplit(x,"_"), function(y) y[length(y)])
  })))

  message("------- GENERATING PARAMETERS -------")
  pars = list(out = out,
              whitelist.file.path = whitelist.file.path,
              wt.max.mismatch = wt.max.mismatch,
              mut.max.mismatch = mut.max.mismatch,
              keep.raw.reads = keep.raw.reads,
              ncores = ncores,
              reverse.complement = reverse.complement,
              testing = testing,
              which.read = which.read,
              atac.barcodes = atac.barcodes,
              primer.sequence = primer.sequence,
              primed.max.mismatch = primed.max.mismatch,
              atac.barcodes.file.path = atac.barcodes.file.path,
              wt.sequence =  wt.sequence,
              mut.sequence = mut.sequence,
              mutation.start = mutation.start,
              mutation.end = mutation.end)

  pars = as.data.frame(pars,stringsAsFactors = F)

  setwd(out)

  message("------- SUBMITTING SLURM JOBS TO CLUSTER -------")
  sjobs <- lapply(chunk.index, function(x){
    pars$out = paste0(out,x,"/")
    slurm_apply(f = MutationCalling, params = pars, jobname = x, nodes = 1, cpus_per_node = 1, slurm_options = soptions, submit = T)
  })

  message("------- DONE! -------")
}
