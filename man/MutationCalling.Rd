\name{MutationCalling}
\alias{MutationCalling}
\title{
  Function to call the single cell genotypes from GoTChA experiments.
}
\description{
  Based on the provided files, matches the cell barcodes from the GoT library with the
  specified list (i.e. white list from 10X), and determines the mutational status based
  on the wild type or mutant sequences.
}
\usage{
  MutationCalling(out = "/path_to_filtered_fastqs/",
                  barcodes.file.path = "/path_to_whitelist/whitelist.txt",
                  wt.max.mismatch = 0,
                  mut.max.mismatch = 0,
                  ncores = 1,
                  reverse.complement = TRUE,
                  testing = FALSE,
                  which.read = "R1",
                  wt.sequence = "CGG",
                  mut.sequence = "CAG",
                  mutation.start = 31,
                  mutation.end = 34,
                  max.distance = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
    Character vector of length one. Path to the folder containing the fastq files from the GoTChA library
  }
  \item{barcodes.file.path}{
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
}
\details{
  The MutationCalling function is the core of the Gotcha package. It takes the fastq files of the specified
  folder, extracts the barcode, read 1 and read 2 sequences merging samples from different lanes if required,
  and proceeds to detecting the presence of the provided sequence. It then verifies that the assigned barcodes
  are present in the provided reference list (i.e. whitelist from 10X, singlecell.csv output from cellRanger).
  Once the reads with the correct barcodes are defined, it calls the genotype for each read based on the
  privided wild type and mutant sequences in the specified region. In the final step, the counts for
  wild type and mutant reads for each cell and the mutant fraction are saved as a list.
}
\value{
  Returns a list containing the mutation calling for each barcode.
}
\author{
  Noelle Wheeler, Franco Izzo, and Robert Myers.
}
