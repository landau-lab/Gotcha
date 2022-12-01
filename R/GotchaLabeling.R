#' Executes functions for the genotyping of cells, starting from WT and MUT read counts.
#'
#' @param path Path to the metadata file containing the genotype read counts
#' @param infile Metadata file name
#' @param gene_id Gene name
#' @param sample_id Sample name
#'
#' @return Dataframe with genotype predictions for each cell in the metadata file
#'
#' @export
#'


GotchaLabeling <-
  function(path = "",
           infile = "",
           gene_id = "",
           sample_id = "") {
    return(
      reticulate::import_from_path(module = "gotcha_labeling", path = find.package('Gotcha'))$GotchaLabeling(
        path = path,
        infile = infile,
        gene_id = gene_id,
        sample_id = sample_id
      )
    )
  }
