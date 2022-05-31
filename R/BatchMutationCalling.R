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
