#' Select mitochondrial mutations in phase with the got-cha targeted mutation
#'
#' @param mitomut_table Dataframe with all mitomutations identified

#' @return Dataframe with mutations which can be used to enhance the Got-cha genotyping efficiency.

FilterGotchaPhasedMitoCalls <- function(mitomut_table){


  MitoVarsAnnotatedFiltered <- mitomut_table %>%
    dplyr::filter(wilcox == "inform." & perc_filter=="perc_pass") #%>%
    #dplyr::filter(MutationCall %in% c("WT","MUT")) #%>%
    #dplyr::mutate(MissingVAF=ifelse(is.na(Heteroplasmy),1,0)) %>%
    #dplyr::group_by(OriginalWhiteListMatch) %>%
    #dplyr::mutate(AtleastOneMissingVAF=sum(MissingVAF)) %>%
    #dplyr::filter(AtleastOneMissingVAF<1) %>%
    #dplyr::select(-MissingVAF,-AtleastOneMissingVAF)

  return(MitoVarsAnnotatedFiltered)

}
