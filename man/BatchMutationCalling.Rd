\name{BatchMutationCalling}
\alias{BatchMutationCalling}
\title{
  Function to parallelize runs of the MutationCalling function in slurm servers.
}
\description{
  It generates a parameters data frame to run the MutationCalling function in each chunk
  of the output of FastqFiltering in parallel using a slurm server, and merges the
  output into a single file.
}
\usage{
  BatchMutationCalling(out = "/path_to_filtered_fastqs/",
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
                  soptions = list(mem = '20g', 'cpus-per-task' = 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
    Character vector of length one. Path to the folder containing the fastq files from the GoTChA library
  }
  \item{whitelist.file.path}{
    Character vector of length one. Path to the text file containing the barcodes from the reference.
  }
  \item{wt.max.mismatch}{
    Integer. Specifies the maximum number of mismatches allowed when searching for the wild type sequence.
    Defaults to cero.
  }
  \item{mut.max.mismatch}{
    Integer. Specifies the maximum number of mismatches allowed when searching for the mutant sequence.
    Defaults to cero.
  }
  \item{keep.raw.reads}{
    Logical. Weather to save the raw read sequences in the output file.
  }
  \item{ncores}{
    Integer. Number of cores to use in parallel processing.
  }
  \item{reverse.complement}{
    Logical. Whether to reverse-complement the cell barcode.
  }
  \item{testing}{
    Logical. Weather to subset the fastq files to 10,000 reads to test the script. Defaults to FALSE.
  }
  \item{which.read}{
    Character vector of length one. Specifies which read to use to look for the mutation.
  }
  \item{atac.barcodes}{
    Logical. Whether to use the barcodes in the provided atac.barcodes.file.path.
  }
  \item{primer.sequence}{
    Character vector of length one. Sequence of the expected primer sequence based on the amplified region.
  }
  \item{primed.max.mismatch}{
    Integer. Number of mismatches allowed when looking for the primer sequence.
  }
  \item{atac.barcodes.file.path}{
    Character vector of lenght one. Specifies the path to the singlecell.csv file from cellRanger.
  }
  \item{wt.sequence}{
    Character vector of length one. Specifies the expected sequence in wild type cells.
  }
  \item{mut.sequence}{
    Character vector of length one. Specifies the expected sequence in mutant cells.
  }
  \item{mutation.start}{
    Integer. Specifies the base position in which to start looking for exact matching of
    wild type or mutant sequence.
  }
  \item{mutation.end}{
    Integer. Specifies the base position in which to finish looking for exact matching of
    wild type or mutant sequence.
  }
  \item{soptions}{
    List containing the parameters for the submission of parallel MutationCalling jobs in a
    slurm cluster.
  }
}
\details{
  The BatchMutationCalling function facilitates parallelization of multiple MutationCalling functions in a
  slurm cluster using slurm_apply. It first creates the data frame containing the parameters for each of
  the jobs to be submitted, and then submits the jobs using the slurm_apply function. The outputs from
  each of the submitted jobs corresponding to a chunk of the fastqs are recovered via the get_slurm_out
  function and merged into a single object containing the wild type read counts, mutant read counts and
  mutant fraction for each detected barcode.
}
\value{
  Returns a list containing the mutation calling for each barcode, and the raw sequences if keep.raw.sequences is TRUE.
}
\author{
  Franco Izzo and Robert Myers
}
