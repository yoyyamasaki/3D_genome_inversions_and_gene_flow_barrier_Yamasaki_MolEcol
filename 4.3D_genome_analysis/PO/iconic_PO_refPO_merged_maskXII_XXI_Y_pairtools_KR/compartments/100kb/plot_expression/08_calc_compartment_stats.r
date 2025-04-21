#Calculate basic statistics of AB compartment.
#Dependency: tidyverse v2.0.0, , bedtoolsr v.2.30.0-6, bedtools v2.31.1
library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/software/bedtools-2.31.1/bedtools2/bin/")

#File path
gff_file <- "~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.mod.maskXII_XXI_Y.gff"
AB_file <- "../iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR_100kb_domains.bed"
AB_settings <- "KR_50kb"
genome <- "~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.maskXII_XXI_Y.fa.order.genome"

#Import data
AB_boundary <- read_tsv(AB_file, col_names = c("chr", "start", "end", "compartment", "eigen")) %>%
  dplyr::mutate(size = end - start + 1)%>%
  dplyr::filter(!(chr %in% c("chrY"))) %>%
  dplyr::filter(eigen != 0)

chr_size <- read_tsv(genome, col_names = c("chr", "size")) %>%
  dplyr::filter(!(chr %in% c("chrY")))  

#A compartment size, percentage, and mean size
A_size_total <- AB_boundary %>%
  dplyr::filter(compartment == "A") %>%
  dplyr::select(size) %>%
  pull(size) %>%
  sum()

A_size_total/sum(chr_size$size)

A_size_mean<- AB_boundary %>%
  dplyr::filter(compartment == "A") %>%
  dplyr::select(size) %>%
  pull(size) %>%
  mean()
  
#B compartment size, percentage, and mean size
B_size_total <- AB_boundary %>%
  dplyr::filter(compartment == "B") %>%
  dplyr::select(size) %>%
  pull(size) %>%
  sum()

(B_size_total)/sum(chr_size$size)

B_size_mean<- AB_boundary %>%
  dplyr::filter(compartment == "B") %>%
  dplyr::select(size) %>%
  pull(size) %>%
  mean()


#x squared test for gene concentration in A compartment
obs_gene_A <- nrow(PO_gff_gene_table_A)
obs_gene_B <- nrow(PO_gff_gene_table_B)
exp_gene_A <- (nrow(PO_gff_gene_table_A) + nrow(PO_gff_gene_table_B)) * (A_size_total/sum(chr_size$size))
exp_gene_B <- (nrow(PO_gff_gene_table_A) + nrow(PO_gff_gene_table_B)) * (B_size_total/sum(chr_size$size))

AB_gene_matrix <- matrix(c(obs_gene_A, obs_gene_B, 
                                exp_gene_A, exp_gene_B), nrow = 2, byrow = T)

chisq.test(AB_gene_matrix)


