\name{DiffLMM}
\alias{DiffLMM}
\title{
  Function to run differential analysis between genotypes within specified cell clusters
}
\description{
  It performs differential analysis of the specified or provided feature matrix (i.e., motif accessibility
  z-scores or gene accessibility scores) within the specified cell cluster, comparing the conditions defined
  in treatmet.column for the labels defined in treatment.levels. Applies a linear mixture model (LMM) which
  models sample identity as random effect to account for potential batch effects.
}
\usage{
DiffLMM(metadata,
        provided.matrix,
        subset = NULL,
        sample.column,
        cluster.column,
        selected.clusters,
        treatment.column,
        treatment.levels,
        ncores = 1)
}

\arguments{
  \item{metadata}{
    A data frame containing the information per cell required for the analysis,
    rownames should be the cell barcodes.
  }
  \item{provided.matrix}{
    Sparse or regular matrix to use as input for the differential analysis.
  }
  \item{subset}{
    Character vector containing the names of a subset of features of interest to be included
    in the differential analysis.
  }
  \item{sample.column}{
    Character vector of length one specifying the name of the column present in the metadata
    of the ArchR project specifying sample identity.
  }
  \item{cluster.column}{
    Character vector of length one specifying the name of the column present in the metadata
    of the ArchR project specifying clusters to be included in the analysis.
  }
  \item{selected.clusters}{
    Character vector containing the names of the clusters to be included in the analysis.
  }
  \item{treatment.column}{
    Character vector of length one specifying the name of the column in the ArchR metadata
    containing the treatment to be compared (i.e., "Genotype").
  }
  \item{treatment.levels}{
    Character vector of length two specifying the treatment levels present in the column of
    the ArchR metadata containing the treatment (i.e., "WildType","Mutant").
  }
  \item{ncores}{
    Integer. Number of cores to use in parallel processing.
  }

}
\details{
  The DiffLMM function provides the means to perform intra-cluster differential feature analysis,
  for example within the Stem cell cluster, between wild-type and mutated cells. By including the
  sample identity as a random effect within the linear mixture model (LMM), it allows to compensate
  for sample-driven effects, and uncover the intra-cluster differences between treatments. The
  output includes P values, fold change calculated as treatment.level 2 / treatment.level 1, delta
  as calculated by treatment.level 2 - treatment.level 1 and adjusted P values via false discovery
  rate (FDR) through Benjamini-Hochberg correction.
}
\value{
  Returns a data.frame class object.
}
\author{
  Franco Izzo
}
