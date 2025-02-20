\name{BatchMutationCalling}
\alias{BatchMutationCalling}
\title{
  Function to parallelize runs of the MutationCalling.
}
\description{
  It generates a parameters data frame to run the MutationCalling function in each chunk
  of the output of FastqFiltering in parallel, and merges the output into a single file.
}
\usage{
  BatchMutationCalling(out = "/path_to_filtered_fastqs/",
                  barcodes.file.path = "/path_to_singlecell.csv"
                  wt.max.mismatch = 0,
                  mut.max.mismatch = 0,
                  ncores = 1,
                  reverse.complement = T,
                  testing = F,
                  which.read = "R1",
                  wt.sequence =  "CGG",
                  mut.sequence = "CAG",
                  mutation.start = 31,
                  mutation.end = 34,
                  max.distance = 2,
                  rewrite = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
    Character vector of length one. Path to the folder containing the fastq files from the GoTChA library.
  }
  \item{barcodes.file.path}{
    Character vector of length one. Path to the text file containing the barcodes from the reference.
  }
  \item{wt.max.mismatch}{
    Integer. Specifies the maximum number of mismatches allowed when searching for the wild type sequence.
    Defaults to zero.
  }
  \item{mut.max.mismatch}{
    Integer. Specifies the maximum number of mismatches allowed when searching for the mutant sequence.
    Defaults to cero.
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
  \item{max.distance}{
    Integer. Specifies the number of mismatches allowed when comparing barcodes with the barcode
    whitelist.
  }
  \item{rewrite}{
    Logical. Specifies whether the output in a chunk should be rewritten if the file already exists.
    Defaults to TRUE.
  }
}
\details{
  The BatchMutationCalling function facilitates parallelization of multiple MutationCalling functions.
  It first creates the data frame containing the parameters for each of the chunks that will be run, and
  then submits each of them in parallel with one core per chunk. The outputs from each of the submitted
  jobs corresponding to a chunk of the fastqs are merged into a single object containing the wild type
  read counts, mutant read counts and mutant fraction for each detected barcode.
}
\value{
  Returns a list containing the mutation calling for each barcode.
}
\author{
  Franco Izzo, Noelle Wheeler, and Robert Myers
}
