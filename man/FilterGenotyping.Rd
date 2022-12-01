\name{FilterGenotyping}
\alias{FilterGenotyping}

\title{
  Function to filter out cells with low genotyping reads counts
}
\description{
  The function takes the specified ArchR project and retrieves the
  total counts from the specified column. It applies the cutoff
  indicated in min.reads and on the remaining distribution, it
  calculates the threshold based on the 10th percentile.
}
\usage{
FilterGenotyping(archrProject,
                 total.counts.column,
                 min.reads = 100,
                 plot = F)
}

\arguments{
  \item{archrProject}{
  ArchRproject class object containing the total genotyping counts
  in the @cellColData slot.
  }
  \item{total.counts.column}{
  Character vector of length one specifying the column present in
  the @cellColData slot of the ArchR object to take the genotyping
  counts from.
  }
  \item{min.reads}{
  Integer. Initial minimum number of reads to apply as threshold
  before defining the median and standard error of the distribution.
  Defaults to 100.
  }
  \item{plot}{
  Logical. Whether to plot the histogram of the initial read
  count distribution. Black dotted line indicates the initial
  min.reads value specified, the red dotted line indicates
  the value for the threshold (10th percentile) used.
  }
}
\value{
  Returns an ArchR project with the specified column plus the
  "_above_read_th_" suffix containing a logical value whether
  to consider the cell genotype call.
}
\author{
  Franco Izzo and Robert Myers
}
