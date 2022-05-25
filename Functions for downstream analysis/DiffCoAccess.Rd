\name{DiffCoAccess}
\alias{DiffCoAccess}
\title{
  Function to calculate co-accessibility within cell subgroups
}
\description{

}
\usage{
DiffCoAccess(archrProject,
             clusters,
             select.clusters = NULL,
             genotype,
             genotype.levels,
             peakSubset = NULL,
             maxDist,
             sample.column = NULL,
             reducedDims = "IterativeLSI",
             corCutOff = 0,
             ncores = 1)
}
\arguments{
   \item{archrProject}{
    ArchR project class object.
  }
  \item{clusters}{
    Characted vector of length one. Specifies the column name where the cluster labels are
    stored in the ArchR object metadata.
  }
  \item{select.clusters}{
    Character vector containing the names of the clusters to be included in the analysis.
  }
  \item{genotype}{
    Characted vector of length one. Specifies the column name where the genotype labels are
    stored in the ArchR object metadata.
  }
  \item{genotype.levels}{
    Character vector containing the genotypes to be included in the analysis.
  }
  \item{peakSubset}{
    SummarizedExperiment or GRanges class object specifying the genomic region to calculate the
    differential co-accessibility loops.
  }
  \item{maxDist}{
    Integer. Maximum distance from the peak of interest to consider pairs of peaks.
  }
  \item{sample.column}{
    Characted vector of length one. Specifies the column name where the sample labels are
    stored in the ArchR object metadata.
  }
  \item{reducedDims}{
    Character vector of length one specifying the dimensionality reduction present in the ArchR
    object to use. Defaults to "IterativeLSI"
  }
  \item{corCutOff}{
    Integer. Defines the minimum correlation value above which to consider co-accessible peaks
  }
  \item{ncores}{
    Integer. Defines the number of cores to use for parallel processing via mclapply.
  }
}
\details{
  The function leverages co-accessibility calculations performed via the ArchR package to re-calculate
  co-accessibility for the indicated genotypes within the specified cell clusters.
}
\value{
  List class object containing the co-accessibility values and loops for each of the cell groups.
}
\references{
 ArchR package (https://www.nature.com/articles/s41588-021-00790-6)
}
\author{
  Franco Izzo
}
