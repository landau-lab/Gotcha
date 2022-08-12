\name{FastqFiltering}
\alias{FastqFiltering}

\title{
  Function to filter fastq files based on base quality scores
}
\description{
  Given a path to the folder containing the output of the FastqSplit function,
  it filters each fastq files for each chunk and stores them in a new folder
  named after the split suffix.
}
\usage{
  FastqFiltering(out, min.quality, min.bases, which.read,
                 read.region, ncores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
  Character string of lenght one corresponding to the full path to the folder containing
  the split fastq files from the GoTChA library.
  }
  \item{min.quality}{
  Integer. Defines the minimum quality score used to filter our reads from the split
  fastq files. Defaults to 15.
  }
  \item{min.bases}{
  Integer. Defines the minimum number of bases that can have a lower quality score than
  specified by min.quality. Defaults to 1.
  }
  \item{which.read}{
  Character vector of lenght one. Specified which read to screen for quality filtering.
  Should match the read specified in the fastq split file names and be written in
  pattern form when specifying more than one read to check for base quality
  (i.e. "R1" or "R1|R2" or "R1|R2|R3").
  }
  \item{read.region}{
  Integer vector of length two. Which part of the read to scan for low quality bases to
  filter on. If NULL, the function scans the entire read for base quality filtering.
  Defaults to NULL.
  }
  \item{ncores}{
  Integer. Specifies the number of cores to use in the mclapply function. Defaults to
  1 (no parallelization).
  }
}
\details{
  The FastqFiltering function reads the fastq.gz files contained in the specified folder,
  and uses the second element of the file names (the split suffix after running FastqSplit)
  to generate an index of unique sets of fastq files that belong to the same chunk. Then
  it reads the fastq files and generates the quality score matrix for the assigned region
  and scans for bases with quality lower than "min.quality". If the number of bases with
  lower than min.quality is higher than the "min.bases" parameter, it will flag the set
  of reads for removal. After filtering, it will generate a new folder "/Filtered/"
  within the specified "out" containing a folder with the filtered outputs for
  each chunk, as identified by the split suffix.
}
\value{
  No value returned. Stores the filtered fastq files in the /Filtered/ folder, creating
  a new folder for each chunk processed.
}
\author{
  Franco Izzo and Robert Myers
}
