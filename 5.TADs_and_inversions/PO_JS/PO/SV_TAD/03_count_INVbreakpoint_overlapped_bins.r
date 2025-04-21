#Test association betwen inversion breakpoints and TAD boundaries.
#Dependency: tidyverse v2.0.0, bedtoolsr v 2.30.0-6
library(tidyverse)
library(bedtoolsr)
options(scipen = 100)

################################################################################
#Read files
inversion <- read_tsv("../syri20240917_PO_JS_asm5syri.INV.bed", col_names = c("chr", "start", "end"))

TAD1 <- read_tsv("~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR/spectraltad/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level1.bed",
                 col_names = c("chr", "start", "end", "level"))
TAD2 <- read_tsv("~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/PO_hifi/iconic_PO_ref_PO_merge_maskXII_XXI_Y_paitrools/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR/spectraltad/iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level2.bed",
                 col_names = c("chr", "start", "end", "level"))

genome <- read_tsv("~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.maskXII_XXI_Y.fa.order.genome", col_names = c("chr", "length")) %>%
  dplyr::mutate(No_bins = trunc(length/25000) + 1) %>%
  dplyr::filter(!(chr %in% c("chrY")))

mask_file <- "~/db/PO_merged_maskXII_XXI_Y_20240517/mask_regions_PAR_chrXIIrep_chrXXIrep.bed"
mask <- read_tsv(mask_file, col_names = c("chr", "start", "end")) %>%
  dplyr::mutate(size = end - start) %>%
  dplyr::filter(!(chr %in% c("chrY"))) %>%
  dplyr::mutate(No_bins = trunc(size/25000))

No_total_bins <- sum(genome$No_bins) - sum(mask$No_bins)

CHR <- paste0("chr", as.roman(seq(1, 21)))

################################################################################
#Mask inversions on the masked sites
inversion <- bt.intersect(inversion, mask, "-wa", "-v")
colnames(inversion) <- c("chr", "start", "end")
inversion_ordered <- tibble()

#Order inversion by chromome
for(i in CHR){
  inversion_ordered <- bind_rows(inversion_ordered, dplyr::filter(inversion, chr == i))
}

#Calculate size and attach ID
inversion_ordered <- inversion_ordered %>%
  dplyr::mutate(size = end - start) %>%
  dplyr::mutate(inversionID = paste0("INVPO_", seq(1, nrow(inversion_ordered))))

#Make bed for all bin
chr_bins <- tibble()
for(i in CHR){
  genome_chr <- genome %>%
    dplyr::filter(chr == i)
  start_bin <- seq(0, genome_chr$length[1], 25000)
  end_bin <- c(seq(25000, genome_chr$length[1], 25000), genome_chr$length[1])
  chr_bins <- bind_rows(chr_bins, tibble(chr = rep(i, length(start_bin)), start = start_bin, end = end_bin))
}

##Make unique TAD boundary list
TAD_boundaries <- bind_rows(TAD1, TAD2) %>%
  dplyr::select(chr, start, end) %>%
  pivot_longer(cols = c(start, end), names_to = "position", values_to = "bp") %>%
  dplyr::mutate(start = bp - 25000, end = bp + 25000) %>%
  dplyr::select(chr, start, end) %>%
  unique()
TAD_boundaries$start[grep("-25000", TAD_boundaries$start)] <- 0

################################################################################
# For all inversions
##Make bed of inversion breakpoints position
inversion_breakpoints_onebase <- inversion_ordered %>%
  pivot_longer(cols = c("start", "end"), names_to = "pos", values_to = "bp") %>%
  dplyr::mutate(start = bp - 1, end = bp ) %>%
  dplyr::select(chr, start, end, inversionID) 

##Count inversion breakpoint overlapped bin number
breakpoint_overlapped_bins <- bt.intersect(chr_bins, inversion_breakpoints_onebase, "-wa")
nrow(breakpoint_overlapped_bins)
breakpoint_overlapped_bins_unique <- unique(breakpoint_overlapped_bins)
Total_bin_No_overlap_INV <- nrow(breakpoint_overlapped_bins_unique)

##Count TAD boundary bin number overlapped with inversion breakpoints
TADboundary_InversionBreakpoint_overlapped_bins <- bt.intersect(TAD_boundaries, inversion_breakpoints_onebase, "-wa")
nrow(TADboundary_InversionBreakpoint_overlapped_bins)
TADboundary_InversionBreakpoint_overlapped_bins_unique <- unique(TADboundary_InversionBreakpoint_overlapped_bins)
bin_TADboundary_overlapped_INV <- nrow(TADboundary_InversionBreakpoint_overlapped_bins_unique)


## x-squired test
InvBreakAll_TadBoundary <- bin_TADboundary_overlapped_INV * 2
InvBreakAll_NonBoundary <- Total_bin_No_overlap_INV * 2 - InvBreakAll_TadBoundary
InvNonBreak_TadBoundary <- nrow(TAD_boundaries) * 2 - InvBreakAll_TadBoundary
InvNonBreak_TadNonBoundary <- No_total_bins - InvNonBreak_TadBoundary - InvBreakAll_NonBoundary - InvBreakAll_TadBoundary

INV_Boundary_matrixAll <- matrix(c(InvBreakAll_TadBoundary, InvBreakAll_NonBoundary, 
                                   InvNonBreak_TadBoundary, InvNonBreak_TadNonBoundary), nrow = 2, byrow = T)

chisq.test(INV_Boundary_matrixAll)

################################################################################
# For 25kb< inversions
##Make bed of inversion breakpoints position
inversion_breakpoints_onebase_25kb <- inversion_ordered %>%
  dplyr::filter(size >= 25000) %>%
  pivot_longer(cols = c("start", "end"), names_to = "pos", values_to = "bp") %>%
  dplyr::mutate(start = bp - 1, end = bp ) %>%
  dplyr::select(chr, start, end, inversionID) 

##Count inversion breakpoint overlapped bin number
breakpoint_overlapped_bins_25kb <- bt.intersect(chr_bins, inversion_breakpoints_onebase_25kb, "-wa")
nrow(breakpoint_overlapped_bins_25kb)
breakpoint_overlapped_bins_25kb_unique <- unique(breakpoint_overlapped_bins_25kb)
Total_bin_No_overlap_INV_25kb <- nrow(breakpoint_overlapped_bins_25kb_unique)

##Count TAD boundary bin number overlapped with inversion breakpoints
TADboundary_InversionBreakpoint_overlapped_bins_25kb <- bt.intersect(TAD_boundaries, inversion_breakpoints_onebase_25kb, "-wa")
nrow(TADboundary_InversionBreakpoint_overlapped_bins_25kb)
TADboundary_InversionBreakpoint_overlapped_bins_25kb_unique <- unique(TADboundary_InversionBreakpoint_overlapped_bins_25kb)
bin_TADboundary_overlapped_INV_25kb <- nrow(TADboundary_InversionBreakpoint_overlapped_bins_25kb_unique)

## x-squired test
InvBreak_25kb_TadBoundary <- (bin_TADboundary_overlapped_INV_25kb * 2) - 1
InvBreak_25kb_NonBoundary <- Total_bin_No_overlap_INV_25kb * 2 - InvBreak_25kb_TadBoundary
InvNonBreak_25kb_TadBoundary <- nrow(TAD_boundaries) * 2 - InvBreak_25kb_TadBoundary
InvNonBreak_25kb_TadNonBoundary <- No_total_bins - InvNonBreak_25kb_TadBoundary - InvBreak_25kb_NonBoundary - InvBreak_25kb_TadBoundary

INV_25kb_Boundary_matrix <- matrix(c(InvBreak_25kb_TadBoundary, InvBreak_25kb_NonBoundary, 
                                   InvNonBreak_25kb_TadBoundary, InvNonBreak_25kb_TadNonBoundary), nrow = 2, byrow = T)

chisq.test(INV_25kb_Boundary_matrix)

