\name{MergeMutationCalling}
\alias{MergeMutationCalling}
\title{
  Merge mutation calling outputs
}
\description{
  Takes the output of the BatchMutationCalling function and
  merges the genotyping data from each of the split fastq
  files.
}
\usage{
MergeMutationCalling(out)
}

\arguments{
  \item{out}{
  Character vector of length one specifying the path to the
  folder containing the results from the BatchMutationCalling
  function.
}
}
\details{
  The function screens for folder names containing the
  "_rslurm" patter and loads the results from each
  of the corresponding folders. It then merges the output
  data frames and collapses the counts by barcode to
  generate a merge output containing the genotyping
  information to be used for downstream applications.
}
\value{
  It creates a new folder in the specified path containing
  the data frame with the merged genotyping outputs.
}
\author{
  Franco Izzo and Robert Myers
}
