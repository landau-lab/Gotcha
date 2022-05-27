#devtools::document()
library(tidyr)
library(dplyr)
library(tidymodels)

source("R/variant_calling.R")


#' Identify mitochondrial mutations in phase with the got-cha targeted mutation
#'
#' @param mgatk RangedSummarizedExperiment object which can be generated following the pipeline in \link[mgatk]{https://github.com/caleblareau/mgatk}.
#' @param genotyping_table Got-cha genotyping table. MutationCall allowed values are only WT,HET,MUT or NA. Rownames should contain the cell barcodes.
#' @param gene Name of the targeted gene of interest. It has to match the column name within the got-cha genotyping table.
#' @param min_cellcov_calculate_af Minimum cell coverage at a given site to generate an heteroplasmic value.
#' @param pattern Suffix added to all the cell/sample names.
#' @param frac_thres Minimum percentage of cells required in at least two out of the three categories of cells containing 0\% (1), 100\% (2) or intermediate heteroplasmic values (3).

#' @return Dataframe with mutations which can be used to enhance the Got-cha genotyping efficiency.
#' @examples
#' load("data/Input.RData", verbose = T)
#' IdentifyGotchaPhasedMitoCalls(mgatk=mgatk, genotyping_table = genotyping_table, gene = "JAK2")
#print(paste("Num cells carrying a variant to trust it:",numcellsupportingvar))
#print(paste("Minimum strand AF correlation to trust a variant:",strand_correlation))
#print(paste("Variance-mean ratio (log-scaled) higher than:",log10_vmr))
#print(paste("Min cell coverage to stabilize VAF:",low_coverage_threshold))
#print(paste("Min coverage to change VAF into NA: ",min_cellcov_calculate_af))
IdentifyGotchaPhasedMitoCalls <- function(mgatk,
                           genotyping_table,
                           gene,
                           min_cellcov_calculate_af = 10,
                           pattern = "-1",
                           frac_thres = 1,
                           numcellsupportingvar=5,
                           strand_correlation=0.65,
                           log10_vmr=-2,
                           low_coverage_threshold = 10,
                           stabilize_variance = TRUE){


  ###########################################################
  # Create a dataframe containing ONLY the genotyped cells  #
  print("Intersecting the mgatk object with the genotyping table...")
  genotyped_cells_table <- tibble::rownames_to_column(genotyping_table, var = "OriginalWhiteListMatch") %>%
    dplyr::mutate(OriginalWhiteListMatch=gsub(".*#","",OriginalWhiteListMatch))
  genotyped_cells <- genotyped_cells_table %>%
    dplyr::select(OriginalWhiteListMatch) %>%
    unlist()
  mgatk_se_subset <- mgatk[,!is.na(match(colnames(mgatk),genotyped_cells))]

  ###########################################################
  # Perform the variant calling                             #
  print("Calling the mitochondrial mutations...")
  mut_se <- call_mutations_mgatk(mgatk_se_subset,
                                 stabilize_variance = stabilize_variance,
                                 low_coverage_threshold = low_coverage_threshold)

  misc_df <- data.frame(rowData(mut_se))
  filter_df <- misc_df %>%
      dplyr::filter(n_cells_conf_detected >= numcellsupportingvar & strand_correlation >= strand_correlation & log10(vmr) > -5) %>%
      # remove known false positives that occur from mapping error, just removing all sites
      # A302C, C309T, C311T, C312T, C313T, G316C, C514A, A515G, A523C, C524G,C3106A, T3109C, C3110A, 3107
      dplyr::filter(!grepl("^302A|^309C|^311C|^312C|^313C|^316G|^514C|^515A|^523A|^524C|^3106C|^3107N|^3109T|^3110C", variant))

  ###########################################################
  # Obtain the AF table                                     #
  # Rows are sites (Site names in last column)              #
  # Columns are cells                                       #
  print("Obtaining heteroplasmies...")

  GetCountsNucleotide <- function(nucleotide) {
    print(nucleotide)
    counts_fw <- as.data.frame(data.matrix(assays(mgatk_se_subset)[[paste(nucleotide,"_counts_fw",sep="")]]))
    counts_rev <- as.data.frame(data.matrix(assays(mgatk_se_subset)[[paste(nucleotide,"_counts_rev",sep="")]]))
    counts <- counts_fw +  counts_rev
    SiteNames <- paste(paste(paste0(data.frame(rowRanges(mgatk_se_subset))[,c(2)]),paste0(rowRanges(mgatk_se_subset)$refAllele),sep=""),nucleotide,sep=">")
    counts$Site <- SiteNames
    counts %>%
      dplyr::filter(!grepl(paste("*",nucleotide,">",nucleotide,sep=""),Site))
  }

  CountsNucleotides <- lapply(c("A","C","G","T"),GetCountsNucleotide)
  AllelesAllSamples <- do.call(rbind,CountsNucleotides)

  ############################################################
  # Get the MT position numbers                              #
  PositionIndexes <- AllelesAllSamples %>%
    dplyr::mutate(Pos=gsub("[A|T|C|G|N].*","",Site)) %>%
    dplyr::select(Pos) %>%
    unlist() %>% paste0() %>% as.numeric()

  ############################################################
  # Get the coverage at all the positions                    #
  cov_sites <- as.data.frame(data.matrix(assays(mgatk_se_subset)[["coverage"]]))
  # downsample the coverage table
  cov <- cov_sites[PositionIndexes,]
  # Avoid 0/0 division errors
  cov[cov == 0] <- NA

  # Last column is "Site" in AllelesAllSamples
  af <- as.matrix(AllelesAllSamples[,-ncol(AllelesAllSamples)]) / as.matrix(cov)
  rownames(af) <- AllelesAllSamples$Site

  MINIMUMCOV_SITE <- min_cellcov_calculate_af
  af_matrix <- data.matrix(af)
  # Change AF at cell sites with less than MINIMUMCOV_SITE reads to NA
  af_matrix[data.matrix(cov) < MINIMUMCOV_SITE] <- NA
  af_table <- as.data.frame(af_matrix)
  af_table$Site <- AllelesAllSamples$Site

  af_filtered <- af_table[filter_df$variant,] %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = contains(pattern),
                        names_to="OriginalWhiteListMatch",
                        values_to=c("Heteroplasmy"))

  table_heteroplasmy <-  genotyped_cells_table %>%
      dplyr::rename(MutationCall=paste(gene,"_Genotype_above_read_th",sep="")) %>%
      dplyr::select(OriginalWhiteListMatch,MutationCall) %>% # CellLine
      dplyr::inner_join(af_filtered, by="OriginalWhiteListMatch") %>%
      dplyr::mutate(MutationCall=ifelse(is.na(MutationCall),"Non-genotyped",ifelse(MutationCall=="NA","Non-genotyped",MutationCall)))

  ############################################################
  # Find informative variants using wilcoxon                 #
  print("Calculating metrics to discern informative variants...")
  table_heteroplasmy2 <- table_heteroplasmy %>%
      dplyr::filter(MutationCall %in% c("MUT","WT")) %>%
      # perform a t-test to check if each site is informative or not
      group_by(Site) %>% do(tidy(wilcox.test(Heteroplasmy ~ MutationCall, data = .))) %>%
      dplyr::mutate(wilcox=ifelse(p.value < 0.05,"inform.","non-inform.")) %>%
      dplyr::select(Site,wilcox) %>%
      dplyr::right_join(table_heteroplasmy)

  #############################################################
  # Find informative variants using prop.test                 #
  #out <- table_heteroplasmy %>%
  #  dplyr::filter(MutationCall %in% c("MUT","WT")) %>%
  #  dplyr::mutate(heteroplasmypres=ifelse(Heteroplasmy > 0,1,Heteroplasmy)) %>%
  #  dplyr::select(-Heteroplasmy) %>%
  #  group_by(Site,MutationCall) %>% dplyr::summarize(counts=sum(heteroplasmypres, na.rm = T)) %>%
  #  tidyr::pivot_wider(names_from = MutationCall, values_from=counts) %>%
  #  dplyr::mutate(WT=sum(MUT,WT)) %>%
  #  # perform a proportion test to check if each site is informative or not (absence presence)
  #  mutate(rate = map2(MUT, WT, ~ prop.test(.x, .y) %>% broom::tidy())) %>%
  #  unnest(rate) %>%
  #  dplyr::select(Site,p.value) %>%
  #  dplyr::mutate(prop.test=ifelse(p.value<0.05,"informative","non-informative")) %>%
  #  dplyr::right_join(table_heteroplasmy2)

  ##############################################################
  # Remove variants not present in a given proportion of cells #
  # Until now present in at least numcellsupportingvar cells   #
  #AF_table1 <- out
  AF_table1 <- table_heteroplasmy2
  AF_table2 <- AF_table1 %>%
    dplyr::mutate(a=ifelse(Heteroplasmy==0,1,0)) %>%
    dplyr::mutate(b=ifelse(Heteroplasmy==1,1,0)) %>%
    dplyr::mutate(c=ifelse(Heteroplasmy>0 & Heteroplasmy<1,1,0)) %>%
    dplyr::mutate(d=ifelse(is.na(Heteroplasmy),1,0)) %>%
    dplyr::group_by(Site) %>%
    dplyr::mutate(a_sum=sum(a, na.rm = T)) %>%
    dplyr::mutate(b_sum=sum(b, na.rm = T)) %>%
    dplyr::mutate(c_sum=sum(c, na.rm = T)) %>%
    dplyr::mutate(na_sum=sum(d, na.rm = T)) %>%
    dplyr::mutate(a_perc=a_sum/(a_sum+b_sum+c_sum)*100) %>%
    dplyr::mutate(b_perc=b_sum/(a_sum+b_sum+c_sum)*100) %>%
    dplyr::mutate(c_perc=c_sum/(a_sum+b_sum+c_sum)*100) %>% # intermediate: 0<x>1
    dplyr::mutate(na_perc=na_sum/(a_sum+b_sum+c_sum+na_sum)*100) %>%
    dplyr::select(-a,-b,-c,-a_sum,-b_sum,-c_sum) %>%
    dplyr::mutate(gen="") %>%
    dplyr::mutate(gen=ifelse(a_perc>=frac_thres,paste0(as.character(gen), "1"),gen)) %>%
    dplyr::mutate(gen=ifelse(b_perc>=frac_thres,paste0(as.character(gen), "1"),gen)) %>%
    dplyr::mutate(gen=ifelse(c_perc>=frac_thres,paste0(as.character(gen), "1"),gen)) %>%
    ungroup() %>%
    dplyr::mutate(genlength=stringr::str_length(gen)) %>%
    dplyr::mutate(perc_filter=ifelse(gen>1,"perc_pass","perc_fail")) %>%
    dplyr::select(-a_perc,-b_perc,-c_perc,-gen,-genlength) %>%
    dplyr::group_by(Site,MutationCall) %>%
    dplyr::mutate(hetmutcall=mean(Heteroplasmy,na.rm=T)) %>%
    dplyr::select(-OriginalWhiteListMatch, -Heteroplasmy) %>%
    ungroup() %>%
    unique() %>%
    tidyr::pivot_wider(names_from = MutationCall, values_from=hetmutcall) %>%
    dplyr::mutate(diff=abs(MUT-WT)) %>%
    dplyr::select(-MUT,-`Non-genotyped`,-WT,-HET,-na_sum,-d) %>% unique()

  AF_table <- AF_table1 %>%
    dplyr::left_join(AF_table2) %>%
    arrange(-diff,na_perc) %>%  # sort dataframe by differences and % of missing data
    mutate(Site = factor(Site, unique(Site)))

  return(AF_table)

}
