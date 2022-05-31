\name{FastqSplit}
\alias{FastqSplit}

\title{
  Function to split .fastq files into chunks for downstream genotyping.
}
\description{
  Given a path to the folder containing the .fastq files from a GoTChA sequencing run, it
  splits each set of fastqs into he specified number of reads for downstream parallel
  computation. The output of the split fastqs is stored in the folder specified output
  folder.
}
\usage{
  FastqSplit(path, out, reads, ncores)
}
\arguments{
  \item{path}{
  Character string of lenght one corresponding to the full path to the folder containing
  the fastq files from the GoTChA library.
  }
  \item{out}{
  Character string of lenght one corresponding to the full path to the folder where to
  store the resulting fastq files after splitting.
  }
  \item{reads}{
  Integer. The number of reads for each chunk to split the original fastqs into. Defaults
  to 2 million reads per chunk.
  }
  \item{ncores}{
  Integer. The number of cores to use for paralellization using mclapply. Defaults to 1
  (no parallelization).
  }
}
\details{
  The FastqSplit function searches for the .fastq.gz pattern in the specified "path"" using
  the "dir" command. If the "out" directory doesn't exists, it creates it using the mkdir
  command through system(). After creating the out folder, it loops through the detected
  fastq files with mclapply and applies zcat, split and gzip commands for each of the files.
  Output folder will contain the fastq files with the added suffix from the split function.
}
\value{
  No value returned. Stores the split fastq files in the specified outs folder.
}
\author{
  Franco Izzo and Robert Myers
}
