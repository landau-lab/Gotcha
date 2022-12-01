\name{GotchaLabeling}
\alias{GotchaLabeling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to run Gotcha genotype labeling from read counts
}
\description{
  Uses the reticulate packge and the r-reticulate-gotcha virtual
  environment to run the noise correction and genotype assignment
  steps in the Gotcha pipeline.
}
\usage{
GotchaLabeling(path = "", infile = "", gene_id = "", sample_id = "")
}
\arguments{
  \item{path}{
  Character vector of length one specifying the path where the metadata
  file is stored
}
  \item{infile}{
  Character vector of length one specifying the file name of the
  provided metadata
}
  \item{gene_id}{
  Character vector of length one specifying the gene or region
  targeted
}
  \item{sample_id}{
  Character vector of length one specifying the sample ID
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a data frame containing the assigned genotype labels
  for each cell and the probability for each of the possible
  genotypes
}
\author{
  Sanjay Kottapalli and Theo Botella
}
