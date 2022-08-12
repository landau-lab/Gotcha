\name{AddGenotypingArchr}
\alias{AddGenotypingArchr}
\title{
  Function to add the mutation calling from of GoTChA into an ArchR object.
}
\description{
  Takes a list of paths to GoTChA outputs and sample names and adds the genotyping
  information to the specified ArchR object metadata.
}
\usage{
AddGenotypingArchr(archrProject,
                    gotchaList,
                    samples,
                    target = "",
                    percentile.use = 95,
                    genotyping.thresholds = c(0.25,0.75),
                    plot = F)
}
\arguments{
  \item{archrProject}{
  An ArchR project class object.
  }
  \item{gotchaList}{
  List of paths to the GoTChA outputs generated with the MutationCalling or
  BatchMutationCalling functions.
  }
  \item{samples}{
  List of sample names matching the order of the gotchaList paths and the
  prefix used in the ArchR project.
  }
   \item{target}{
  Character vector of length one. Specifies the prefix to include in the
  column names when adding the genotyping information to the ArchR object
  metadata. Defaults to NULL.
  }
  \item{percentile.use}{
  Integer between 0 (no correction) and 99. Which percentile to use as
  limit for defining the value of the noise in empty droplets. Defaults
  to 90 (90th percentile).
  }
   \item{genotyping.thresholds}{
  Numeric integer vector of length two. Specifies the minimum and maximum
  range to define wildtype and mutant cells (i.e. c(0.25,0.75) will define
  cells with mutant fraction lower than 0.25 as wildtype and cells with
  mutant fraction higher than 0.75 as mutant.
  }
   \item{plot}{
  Logical. Whether to plot the histogram to look at distributions
  of wildtype and mutant reads.
  }
}
\details{
  The function loads the specified GoTChA outputs, adds the sample name
  prefix and "-1" suffix to the barcodes based on the provided lists. Then
  it adds the genotyping information to the ArchR object meta data at the
  @cellColData slot of the specified ArchR object. It generates the columns
  WTcount, MUTcount and MUTfraction containing the genotyping information.
  If the target is specified (i.e. "JAK2") the name will be appended to the
  column names in the metadata (JAK2_WTcounts, JAK2_MUTcounts, and
  JAK2_MUTfraction, respectively).
  In addition, noise correction is performed, the noise from wild type and mutant
  read counts will be estimated independently and for each sample, using the non-cell
  barcodes ("empty droplets"). The estimated noise counts are  substracted from the
  final genotyping counts. The resulting counts will be added as additional columns
  in the metadata with the "_noise_corrected" suffix. Estimation of the noise is made
  by taking the value of the specified percentile (defaults to 95th percentile) to avoid
  inclusion of outlier low quality cells taken as empty droplets with high genotyping
  counts.
}
\value{
  Returns an ArchR project object with the genotyping information incluede in the
  @cellColData slot.
}

\author{
  Franco Izzo and Robert Myers
}
