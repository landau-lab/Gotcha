#' Optimize the number of variants to use in the classifier given their percentage of missing data
#'
#' @param nsites Number of sites to consider. Between 1 and the total number of sites identified as in phase.
#' @param mitomutations Table containing mitomutations and annotations about their informativeness.
#' @return Dataframe with mutations which can be used to enhance the Got-cha genotyping efficiency.
#' @examples
#'
MissingDataCum <- function(nsites,
                           mitomutations, pattern="-1"){
  hetsum_percell <- mitomutations %>%
    dplyr::filter(perc_filter=="perc_pass") %>%
    dplyr::filter(wilcox=="inform.") %>%
    dplyr::arrange(-na_perc) %>%
    dplyr::select(-MutationCall) %>%
    tidyr::pivot_wider(names_from = OriginalWhiteListMatch, values_from=Heteroplasmy) %>%
    ungroup() %>%
    dplyr::select(contains(pattern))
  hetsum_percell <- hetsum_percell[1:nsites,] %>%
    colSums()
  missing=sum(is.na(hetsum_percell))/length(hetsum_percell)*100
  return(missing)
}

