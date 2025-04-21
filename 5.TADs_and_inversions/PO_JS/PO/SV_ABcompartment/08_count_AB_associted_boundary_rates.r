#Test association betwen inversion breakpoints and AB boundaries.
#Dependency: tidyverse v2.0.0
library(tidyverse)

#Read file
genome <- read_tsv("~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.maskXII_XXI_Y.fa.order.genome", col_names = c("chr", "length")) %>%
  dplyr::mutate(No_bins = trunc(length/25000)) %>%
  dplyr::filter(!(chr %in% c("chrY")))

mask_file <- "~/db/PO_merged_maskXII_XXI_Y_20240517/mask_regions_PAR_chrXIIrep_chrXXIrep.bed"
mask <- read_tsv(mask_file, col_names = c("chr", "start", "end")) %>%
  dplyr::mutate(size = end - start) %>%
  dplyr::filter(!(chr %in% c("chrY"))) %>%
  dplyr::mutate(No_bins = trunc(size/25000))

AB <- read_tsv("/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR/compartments/100kb/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR_100kb_domains.bed",
                               col_names = c("chr", "start", "end", "compartment", "eigenvalue", "note")) 
AB$start <- AB$start -1
AB_boundary_unique <- AB %>%
  pivot_longer(cols = c("start", "end"), names_to = "start_end", values_to = "bp") %>%
  dplyr::select(chr, bp) %>%
  unique()

No_total_bins = sum(genome$No_bins) - sum(mask$No_bins)

No_INV = 16
No_INV_breaks = No_INV *2
No_bins_associated_INV_breaks = No_INV_breaks *2

nrow(AB_boundary_unique)*2 / No_total_bins
 
#x squared test based on the number of bins
# Inversions >25kb
# Number of bins
Unique_ABBoundary <- nrow(AB_boundary_unique)
InvBreak_ABBoundary <- 7*2
InvBreak_NonBoundary <- 32*2 - InvBreak_ABBoundary
InvNonBreak_ABBoundary <- Unique_ABBoundary * 2 - InvBreak_ABBoundary
InvNonBreak_ABNonBoundary <- No_total_bins - InvNonBreak_ABBoundary - InvBreak_NonBoundary

INV_Boundary_matrix <- matrix(c(InvBreak_ABBoundary, InvBreak_NonBoundary, 
                                InvNonBreak_ABBoundary, InvNonBreak_ABNonBoundary), nrow = 2, byrow = T)

chisq.test(INV_Boundary_matrix)

# Inversions All
InvBreakAll_ABBoundary <- 14*2
InvBreakAll_NonBoundary <- 118*2 - InvBreakAll_ABBoundary
InvNonBreak_ABBoundary <- Unique_ABBoundary * 2 - InvBreakAll_ABBoundary
InvNonBreak_ABNonBoundary <- No_total_bins - InvNonBreak_ABBoundary - InvBreak_NonBoundary

INV_Boundary_matrixAll <- matrix(c(InvBreakAll_ABBoundary, InvBreakAll_NonBoundary, 
                                   InvNonBreak_ABBoundary, InvNonBreak_ABNonBoundary), nrow = 2, byrow = T)

chisq.test(INV_Boundary_matrixAll)


