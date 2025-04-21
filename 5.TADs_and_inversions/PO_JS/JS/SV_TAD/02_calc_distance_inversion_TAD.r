#Make table and plot which describe the distance between inversion breakpoints and the nearest TAD boundaries.
#Dependency: tidyverse v2.0.0, bedtoolsr v 2.30.0-6
library(tidyverse)
library(bedtoolsr)
options(scipen = 100)

################################################################################
#Read files
inversion <- read_tsv("../syri20240917_PO_JS_asm5syri.inversion.bedpe", col_names = c("chr_PO", "start_PO", "end_PO", "chr_JS", "start_JS", "end_JS")) %>%
  dplyr::mutate(length_PO = end_PO - start_PO, length_JS = end_JS - start_JS)

TAD1 <- read_tsv("~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/JS/iconhic_JS_refJS_merge_maskXII_XXI_Y_paitrools/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_KR/spectraltad/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level1.bed",
                 col_names = c("chr", "start", "end", "level"))
TAD2 <- read_tsv("~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/JS/iconhic_JS_refJS_merge_maskXII_XXI_Y_paitrools/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_KR/spectraltad/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level2.bed",
                 col_names = c("chr", "start", "end", "level"))

TAD <- bind_rows(TAD1, TAD2)

mask_file <- "~/db/JS_male_merged_mask_chrXII_chrXXI_chrYXI_20240517/mask_chrXII_chrXXI_chrY_IX.bed"
mask <- read_tsv(mask_file, col_names = c("chr", "start", "end")) %>%
  dplyr::mutate(size = end - start) %>%
  dplyr::filter(!(chr %in% c("chrY_IX"))) 

CHR <- paste0("chr", as.roman(seq(1, 21)))

#Define function to calculate the distance between inversion breakpoint and the nearest TAD boundary
find_closest_to_zero <- function(vec){
  vec[which.min(abs(vec))]
}

#Order inversions in chromosome number
inversion_ordered <- tibble()
for(i in CHR){
  inversion_ordered <- bind_rows(inversion_ordered, dplyr::filter(inversion, chr_PO == i))
}

#Add inversion ID
inversion_ordered <- inversion_ordered %>%
  dplyr::mutate(inversionID = paste0("INV", seq(1, nrow(inversion_ordered))))

#Extract 25kb < inversions
inversion_25kb <- inversion_ordered %>%
  dplyr::filter(length_PO >= 25000)
write_tsv(inversion_25kb, "./syri20240917_PO_JS_asm5syri.25kb.INV.POref.bed", col_names = F)

#Extract G. nipponics based inversion corrdinate
inversion_25kb_JS <- inversion_25kb %>%
  dplyr::select(chr_JS, start_JS, end_JS, length_JS, inversionID)
colnames(inversion_25kb_JS) <- c("chr", "start", "end", "length", "inversionID")

#Make bed file which record inversion breakpoint +- 25kb 
inversion_JS_25kb_breakpoints <- inversion_25kb %>%
  dplyr::select(chr_JS, start_JS, end_JS, length_JS, inversionID) %>%
  pivot_longer(cols = c("start_JS", "end_JS"), names_to = "pos", values_to = "bp") %>%
  dplyr::mutate(start = bp - 25000, end = bp + 25000) %>%
  dplyr::select(chr_JS, start, end, inversionID)  %>%
  dplyr::arrange(chr_JS, start)

colnames(inversion_JS_25kb_breakpoints) <- c("chr", "start", "end", "inversionID")
  
inversion_JS_25kb_breakpoints_merged <- bt.merge(inversion_JS_25kb_breakpoints)

write_tsv(inversion_JS_25kb_breakpoints_merged, "./syri20240917_PO_JS_asm5syri.25kb.INV_JSbreakpoint.POref.bed", col_names = F)

################################################################################
#Inversion over 25000bp
#Make a table which describe the distance between inbersion breakpoints and the nearest TAD boundaries.
TAD_overlap <- tibble()
inversion_TAD_distance <- tibble()
for(i in 1:nrow(inversion_25kb_JS)){
  
  #Extract inversion overlapped TADs
  TAD_tmp <- bt.intersect(TAD, inversion_25kb_JS[i, ], "-wa")
  
  #If overlapped TAD existed
  if(nrow(TAD_tmp) > 0){
    colnames(TAD_tmp) <- c("chr", "start", "end", "level")
    
    TAD_tmp <- TAD_tmp %>%
      dplyr::mutate(inversionID = pull(inversion_25kb_JS[i, "inversionID"]))
    
    TAD_overlap <- bind_rows(TAD_overlap, TAD_tmp)
    
    #Calculate distance between inversion breakpoint and the nearest TAD
    front_distance <- find_closest_to_zero(pull(inversion_25kb_JS[i, "start"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_25kb_JS[i, "end"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))
    
    #Convert distance from bp to bin size
    inversion_TAD_dist_tmp <- inversion_25kb_JS[i, ] %>%
      dplyr::mutate(front_distance_JS = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance_JS = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_TAD_distance <- bind_rows(inversion_TAD_distance, inversion_TAD_dist_tmp)
  
  #If overlapped TAD did'nt existed
  }else if(nrow(TAD_tmp) == 0){
    
    #Extract TAD which is on the same chromosome with the focal inversion
    TAD_tmp <- TAD %>%
      dplyr::filter(chr == pull(inversion_25kb_JS[i, "chr"]))
    
    #Calculate distance between inversion breakpoint and the nearest TAD
    front_distance <- find_closest_to_zero(pull(inversion_25kb_JS[i, "start"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_25kb_JS[i, "end"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))

    #Convert distance from bp to bin size
    inversion_TAD_dist_tmp <- inversion_25kb_JS[i, ] %>%
      dplyr::mutate(front_distance = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_TAD_distance <- bind_rows(inversion_TAD_distance, inversion_TAD_dist_tmp)
  }
}

write_tsv(inversion_TAD_distance, "./inversion_JS_asm5_TAD_PO_dist.POref.txt")
write_tsv(TAD_overlap, "./inversion_JS_asm5_overlapped_TAD_PO.POref.txt")

#Plot histogram of distance between inversion breakpoints and the nearest TAD boundaries
inversion_TAD_distance_longer <- inversion_TAD_distance %>%
  pivot_longer(cols = c("front_bin_JS", "rear_bin_JS"), names_to = "front_rear", values_to = "N_bin") 
inversion_TAD_distance_longer$N_bin <- abs(inversion_TAD_distance_longer$N_bin)
inversion_TAD_distance_longer_orthologous <- inversion_TAD_distance_longer

inversion_TAD_distance_hist <- hist(inversion_TAD_distance_longer_orthologous$N_bin, plot = F, breaks = seq(0, 73))

INV_TAD_dist_hist <- ggplot(tibble(breaks = inversion_TAD_distance_hist$breaks[-1], counts = inversion_TAD_distance_hist$counts)) +
  geom_col(aes(x = breaks, y = counts)) + 
  ggtitle("Distance between breakpoint and TAD boundary") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20), limits = c(0.5, 20)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(0, 13)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1, 
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )

ggsave("INV_TAD_distans_bin_JS.POref.eps", INV_TAD_dist_hist, width =5, height = 5)

inversion_TAD_distance_hist$counts[1]/sum(inversion_TAD_distance_hist$counts)

################################################################################
#All inversion
#Extract only G. nipponics corrdinate information
inversion_ordered_JS <- inversion_ordered %>%
  dplyr::select(chr_JS, start_JS, end_JS, length_JS, inversionID)
colnames(inversion_ordered_JS) <- c("chr", "start", "end", "length", "inversionID")

#Remove inversions on masked regions
inversion_ordered_JS <- bt.intersect(inversion_ordered_JS, mask, "-wa", "-v")
colnames(inversion_ordered_JS) <- c("chr", "start", "end", "length", "inversionID")
inversion_ordered_JS <- as_tibble(inversion_ordered_JS)

#Plot inversion size distribution
Inversion_dist_JS <- ggplot(inversion_ordered_JS) +
  geom_histogram(aes(x = length), binwidth = 0.25) +
  scale_x_log10(breaks = c(1e2, 1e3, 1e4, 25000, 1e5, 1e6, 1e7)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12), limits = c(0, 11)) +
  ggtitle("Size and number of inversions: JS_PO") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1, 
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )

ggsave("INV_dist_JS.POref.eps", Inversion_dist_JS, width =5, height = 5)

#Make a table which describe the distance between inbersion breakpoints and the nearest TAD boundaries.
TAD_overlap_all <- tibble()
inversion_TAD_distance_all <- tibble()
for(i in 1:nrow(inversion_ordered_JS)){

  #Extract inversion overlapped TADs
  TAD_tmp <- bt.intersect(TAD, inversion_ordered_JS[i, ], "-wa")

  #If overlapped TAD existed
  if(nrow(TAD_tmp) > 0){
    colnames(TAD_tmp) <- c("chr", "start", "end", "level")
    
    TAD_tmp <- TAD_tmp %>%
      dplyr::mutate(inversionID = pull(inversion_ordered_JS[i, "inversionID"]))
    
    TAD_overlap_all <- bind_rows(TAD_overlap_all, TAD_tmp)

    #Calculate distance between inversion breakpoint and the nearest TAD
    front_distance <- find_closest_to_zero(pull(inversion_ordered_JS[i, "start"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_ordered_JS[i, "end"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))

    #Convert distance from bp to bin size
    inversion_TAD_dist_tmp <- inversion_ordered_JS[i, ] %>%
      dplyr::mutate(front_distance = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_TAD_distance_all <- bind_rows(inversion_TAD_distance_all, inversion_TAD_dist_tmp)

  #If overlapped TAD did'nt existed
    }else if(nrow(TAD_tmp) == 0){

    #Extract TAD which is on the same chromosome with the focal inversion
      TAD_tmp <- TAD %>%
      dplyr::filter(chr == pull(inversion_ordered_JS[i, "chr"]))

    #Calculate distance between inversion breakpoint and the nearest TAD
    front_distance <- find_closest_to_zero(pull(inversion_ordered_JS[i, "start"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_ordered_JS[i, "end"]) - unique(c(TAD_tmp$start, TAD_tmp$end)))

    #Convert distance from bp to bin size
    inversion_TAD_dist_tmp <- inversion_ordered_JS[i, ] %>%
      dplyr::mutate(front_distance = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_TAD_distance_all <- bind_rows(inversion_TAD_distance_all, inversion_TAD_dist_tmp)
  }
}

write_tsv(inversion_TAD_distance_all, "./inversion_JS_asm5_TAD_JS_dist_all.POref.txt")
write_tsv(TAD_overlap_all, "./inversion_JS_asm5_overlapped_TAD_JS_all.POref.txt")

#Plot histogram of distance between inversion breakpoints and the nearest TAD boundaries
inversion_TAD_distance_longer_all <- inversion_TAD_distance_all %>%
  pivot_longer(cols = c("front_bin_JS", "rear_bin_JS"), names_to = "front_rear", values_to = "N_bin") 
inversion_TAD_distance_longer_all$N_bin <- abs(inversion_TAD_distance_longer_all$N_bin)

inversion_TAD_distance_all_hist <- hist(inversion_TAD_distance_longer_all$N_bin, plot = F, breaks = seq(0, 73))

INV_TAD_dist_all_hist <- ggplot(tibble(breaks = inversion_TAD_distance_all_hist$breaks[-1], counts = inversion_TAD_distance_all_hist$counts)) +
  geom_col(aes(x = breaks, y = counts)) + 
  ggtitle("Distance between breakpoint and TAD boundary") +
  scale_x_continuous(breaks = c(1, seq(5, 75, 5)), limits = c(0.5, 20)) +
  scale_y_continuous(breaks = seq(0, 35, 5)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1, 
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )

ggsave("INV_TAD_distans_bin_JS_all.POref.eps", INV_TAD_dist_all_hist, width =5, height = 5)

