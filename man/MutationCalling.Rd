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
                  mutation.end = 34)
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
}
\details{
  The MutationCalling function is the core of the Gotcha package. It takes the fastq files of the specified
  folder, extracts the barcode, read 1 and read 2 sequences merging samples from different lanes if
  required, and proceeds to detecting the presence of the provided sequence associated with correct priming
  ("priming.sequence"). Once the correctly primed reads are found, it verifies that the assigned barcodes
  are present in the provided reference list (i.e. whitelist from 10X, singlecell.csv output from cellRanger).
  Once the reads with the correct barcodes are defined, it calls the genotype for each read based on the
  privided wild type and mutant sequences in the specified region. In the final step, the counts for
  wild type and mutant reads for each cell and the mutant fraction are saved as a list.
}
\value{
  Returns a list containing the mutation calling for each barcode, and the raw sequences if keep.raw.sequences is TRUE.
}
\author{
  Franco Izzo and Robert Myers
}
