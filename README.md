# [GoT-ChA](https://www.biorxiv.org/content/10.1101/2022.05.11.491515v1): Genotyping of Targeted loci with single-cell Chromatin Accessibility
### Why use GoT-ChA?
Somatic mutations are crucial for cancer initiation and evolution, and have been identified across a number of healthy tissues in the human body. These mutations can disrupt normal cellular functions, leading to aberrant clonal expansions via acquired fitness advantages or skewed differentiation topologies. GoT-ChA and similar methods (e.g. GoT, TARGET-seq) aim to pair targeted genotyping with single-cell sequencing approaches in order to understand the impact of somatic mutations directly in human patient samples, in both malignant and non-malignant contexts.

GoT-ChA is unique in that it is a high-throughput single-cell method that pairs targeted genotyping with chromatin accessibility measurements based on the broadly utilized scATAC-seq platform from 10x Genomics. Previous single-cell genotyping approaches were largely based on scRNA-seq protocols, utilizing expressed and captured mRNA transcripts as sources for genotyping information. This results in a limiting dependence on gene expression and mutation location (due to 3' end bias of most scRNA-seq methods), precluding the usage of such technologies on lowly expressed or inconveniently located loci of interest. GoT-ChA surmounts these limitations via direct utilization of genomic DNA for genotyping, whilst simultaneously assaying chromatin accessibility.  

### How does GoT-ChA work?
![image](https://user-images.githubusercontent.com/38476687/170100937-117d5c2e-78cf-4f68-a710-4cb079e7a471.png)
In order to capture genotypes within droplet-based scATAC-seq, two GoT-ChA primers are added into the cell barcoding PCR reaction that are designed to flank the locus of interest. One primer contains the partial Nextera Read 1N sequence in its handle, which allows for GoT-ChA amplicons to interact with the 10x Genomics gel bead oligonucleotide and obtain a unique cell barcode, just as tagmented chromatin fragments do. Further, the second GoT-ChA primer allows for exponential amplification of GoT-ChA amplicons while tagmented chromatin fragments are only linearly amplified. 

![image](https://user-images.githubusercontent.com/38476687/170101600-a33bab72-7b26-436a-b042-4ea40fa0fa4d.png)
After the single-cell emulsion is broken, a small portion of the sample is taken for GoT-ChA library construction, comprised of a hemi-nested PCR, biotin-streptavidin pull-down, and an on-bead sample indexing PCR. The final GoT-ChA library can be pooled with scATAC-seq libraries and sequenced together using standard ATAC sequencing parameters.

### Designing a GoT-ChA experiment
![image](https://user-images.githubusercontent.com/38476687/170109264-8010c7cc-8ee7-4149-98f8-e1273a69d7d5.png)
To utilize GoT-ChA, three primers need to be designed that flank the genomic region of interest. GoT-ChA_R1 (containing the partial Nextera Read1 sequence in its handle) and GoT-ChA_Rev flank the loci of interest, ideally forming an amplicon between 200-500 bps in size. GoT-ChA_nested is utilized in the hemi-nested PCR during GoT-ChA library construction, and crucially needs to bind within 50bp (inclusive of binding site) of the mutation to be covered with standard scATAC-seq sequencing parameters.

## Overview of the Gotcha R package
The Gotcha R package is designed to provide a pipeline for end-to-end processing of the sequecing libraries generated through the GoT-ChA method. It containes functions for both data pre-processing (i.e., filtering of fastq files by base quality, split of fastqs into chunks for parallel processing) as well as mutation calling and noise correction (estimation of noise from empty droplets) and genotype calling. In addition, we have included functions for downstream analysis that work seemingly with the ArchR (https://www.archrproject.com) single cell ATAC-seq processing package, and allow to calculate differential motif accessibility and gene accessibility scores, as well as calculation of co-accessibility for subgroups of cells and visualization of genome tracks. 

### Data pre-processing
In order to pre-process the fastq files generated through GoT-ChA genotyping, we included two functions within the Gotcha R package as follows:

1 - FastqFiltering: reads in the fastqs generated by gotcha, and filters them based on pre-determined base quality scores.

2 - FastqSplit: if parallelization is available through a slurm cluster, this function allows to split the fastq files into chunks for parallel computing of by the BatchMutationCalling function. 


These two functions allow to set up the initial files for running the mutation calling step.

### Mutation calling
To perform mutation calling, we provide two main functions. The core of the Gotcha R package is the ability to genotype each read and map it to the corresponding cell barcode. We apply pattern matching including the primer sequence, to guarantee the detection of wild-type or mutant reads containing the primer sequence ("primed reads"). Only those reads that were correctly primed are kept in downstream steps. Once the genotype for each read has been defined, we calculate the counts of reads per cell for each of the specified genotypes. As output, we provide a data frame containing the cell barcode, wild-type read counts and mutant reads counts and fraction of mutant reads for each cell. 

3 - MutationCalling: this is the core function of the Gotcha R package. It allows local parallelization via multiple cores.

4 - BatchMutationCalling: this function allows to speed up computation times in a slurm cluster. It allows to apply the MutationCalling function in the fastq files split by the FastqSplit function, and submit a job to the cluster for each file. Using this function, processing a GoT-ChA genotyping library of 40 million reads takes ~ 4h. As output, it generates a data frame within each folder of the fastq split files, that can be merged using the followin function.

5 - MergeMutationCalling: this function allows to merge the outputs of the BatchMutationCalling present in each split folder into one single file. This is the final output of the mutation calling step, and can be merged to the metadata of scATAC-seq objects as generated by Signac or ArchR.

### Defining genotype calls
Once we have measured the abundance of wild-type and mutant reads for each cell barcode present in the experiment, we need to define the genotype for each cell. This can be done via two different noise correction methods: Cluster-based genotyping or Empty droplet-based genotyping.

#### Cluster-based background noise correction and genotype assignment


#### Empty droplet-based noise correction and genotype assignment

### Downstream analysis
## Testing GoT-ChA
### Sample data for testing the Gotcha pipeline
Initial fastq files for the genotyping GoT-ChA library as well as metadata containing the final genotype calls are publicly available at GEO (GSE203251).
