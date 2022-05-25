\name{PlotDCA}
\alias{PlotDCA}
\title{
  Function to plot the output of Differential Co-Accessibility (DCA) analysis.
}
\description{
  Takes the output of the DiffCoAccess function and generates track plots
  for visualization of the co-accessibility for each peak within the specified
  genomic region.
}
\usage{
  PlotDCA(archrProject,
          DCA,
          coords,
          features,
          condition.column = NULL,
          condition = NULL,
          cluster.genotype.column,
          cluster.genotype.levels,
          min.cor = 0.1,
          path = NULL,
          plot = TRUE)
}
\arguments{
  \item{archrProject}{
    ArchR project class object.
  }
  \item{DCA}{
    List obtained as output of the DifferentialCoAccess function.
  }
  \item{coords}{
    GRanges class object containing the genomic coordinates of the region of interest.
  }
  \item{features}{
    GRanges list class object containing features to be included in the track plot (i.e.,
    motif sites, gene annotations, locus control regions)
  }
  \item{condition.column}{
    Character vector of length one specifying the name of the column in the ArchR metadata
    containing a condition to subset the cells included (i.e., "Tumor type"). Defaults to NULL.
  }
  \item{condition}{
    Character vector specifying the levels of the condition column to be included in the
    analysis (i.e., "Polycythemia Vera","Myelofibrosis"). Defaults to NULL.
  }
  \item{cluster.genotype.column}{
    Character vector of length one specifying the name of the column in the ArchR metadata
    containing a cluster and treatment labels.
  }
  \item{cluster.genotype.levels}{
    Character vector of length one specifying the levels of the column in the ArchR metadata
    containing a cluster and treatment labels to be compared (i.e., "EP1_WT","EP1_MUT").
  }
  \item{min.cor}{
    Integer specifying the minimum correlation value to be considered towards the count
    of co-acceessible pairs for a given peak. Defaults to 0.1.
  }
  \item{path}{
    Character vector specifying the path to the folder where to store the track plots.
    Defaults to NULL.
  }
  \item{plot}{
    Logical. Whether to plot the tracks of the defined region. Defaults to TRUE.
  }
}
\details{
  Takes the output of the DiffCoAccess function and plots the corresponding genomic tracks.
}
\value{
  Generates track plots for each of the peaks and its corresponding pairs within the specified
  genomic region between the specified cell groups. It also outputs a plot showing the difference
  in number of co-accessible events per peak between the specified conditions.
}
\author{
  Franco Izzo
}
