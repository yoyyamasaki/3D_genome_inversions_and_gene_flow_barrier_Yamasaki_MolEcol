#Make table and plot which describe the distance between inversion breakpoints and the nearest AB compartment boundaries.
#Dependency: tidyverse v2.0.0, bedtoolsr v 2.30.0-6
library(tidyverse)
library(bedtoolsr)

################################################################################
#Read files
inversion <- read_tsv("./inversion_JS_PObased.bed", col_names = c("chr", "start", "end"))

AB <- read_tsv("/home/yo/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/JS/iconhic_JS_refJS_merge_maskXII_XXI_Y_paitrools/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_KR/compartments/100kb/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_KR_100kb_domains.bed",
                 col_names = c("chr", "start", "end", "compartment", "eigenvalue", "note")) %>%
  dplyr::mutate(size = end - start + 1)%>%
  dplyr::filter(!(chr %in% c("chrIX_Y"))) %>%
  dplyr::filter(eigenvalue != 0)

mask_file <- "~/db/JS_male_merged_mask_chrXII_chrXXI_chrYXI_20240517/mask_chrXII_chrXXI_chrY_IX.bed"
mask <- read_tsv(mask_file, col_names = c("chr", "start", "end")) %>%
  dplyr::mutate(size = end - start) %>%
  dplyr::filter(!(chr %in% c("chrY"))) 

CHR <- paste0("chr", as.roman(seq(1, 21)))

#Define function to calculate the distance between inversion breakpoint and the nearest TAD boundary
find_closest_to_zero <- function(vec){
  vec[which.min(abs(vec))]
}
  
#Extact inversions which are not overlapped with masked regions
inversion <- bt.intersect(inversion, mask, "-wa", "-v")
colnames(inversion) <- c("chr", "start", "end")

#Order inversions in chromosome number
inversion_ordered <- tibble()
for(i in CHR){
  inversion_ordered <- bind_rows(inversion_ordered, dplyr::filter(inversion, chr == i))
}

#Add inversion ID
inversion_ordered <- inversion_ordered %>%
  dplyr::mutate(size = end - start) %>%
  dplyr::mutate(inversionID = paste0("INVPO_", seq(1, nrow(inversion_ordered))))

#Extract 25kb < inversions
inversion_25kb <- inversion_ordered %>%
  dplyr::filter(size >= 25000)
write_tsv(inversion_25kb, "./syri20240917_PO_JS_asm5syri.25kb.INV.bed", col_names = F)

#Make bed file which record inversion breakpoint +- 25kb 
inversion_25kb_breakpoints <- inversion_25kb %>%
  pivot_longer(cols = c("start", "end"), names_to = "pos", values_to = "bp") %>%
  dplyr::mutate(start = bp - 25000, end = bp + 25000) %>%
  dplyr::select(chr, start, end, inversionID) 

inversion_25kb_breakpoints_merged <- bt.merge(inversion_25kb_breakpoints)

write_tsv(inversion_25kb_breakpoints_merged, "./syri20240917_JS_asm5syri.25kb.INVJS_breakpoint.POref.bed")

################################################################################
#Inversion over 25000bp
#Make a table which describe the distance between inbersion breakpoints and the nearest AB boundaries.
AB_overlap <- tibble()
inversion_AB_distance <- tibble()
for(i in 1:nrow(inversion_25kb)){
  #Extract inversion overlapped AB    
  AB_tmp <- bt.intersect(AB, inversion_25kb[i, ], "-wa")

  #If overlapped AB existed    
  if(nrow(AB_tmp) > 0){
  colnames(AB_tmp) <- c("chr", "start", "end", "compartment", "eigenvalue", "note", "size")
  
  AB_tmp <- AB_tmp %>%
    dplyr::mutate(inversionID = pull(inversion_25kb[i, "inversionID"]))
  
  AB_overlap <- bind_rows(AB_overlap, AB_tmp)

  #Calculate distance between inversion breakpoint and the nearest AB    
  front_distance <- find_closest_to_zero(pull(inversion_25kb[i, "start"]) - unique(c(AB_tmp$start, AB_tmp$end)))
  rear_distance <- find_closest_to_zero(pull(inversion_25kb[i, "end"]) - unique(c(AB_tmp$start, AB_tmp$end)))

  #Convert distance from bp to bin size    
  inversion_AB_dist_tmp <- inversion_25kb[i, ] %>%
    dplyr::mutate(front_distance_JS = front_distance,
                  front_bin_JS = front_distance/25000,
                  rear_distance_JS = rear_distance,
                  rear_bin_JS = rear_distance/25000)

  inversion_AB_distance <- bind_rows(inversion_AB_distance, inversion_AB_dist_tmp)
  
  #If overlapped AB did'nt existed
  }else if(nrow(AB_tmp) == 0){

    #Extract AB which is on the same chromosome with the focal inversion    
    AB_tmp <- TAD %>%
      dplyr::filter(chr == pull(inversion_25kb[i, "chr"]))

    #Calculate distance between inversion breakpoint and the nearest AB        
    front_distance <- find_closest_to_zero(pull(inversion_25kb[i, "start"]) - unique(c(AB_tmp$start, AB_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_25kb[i, "end"]) - unique(c(AB_tmp$start, AB_tmp$end)))

    #Convert distance from bp to bin size        
    inversion_AB_dist_tmp <- inversion_25kb[i, ] %>%
      dplyr::mutate(front_distance_JS = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance_JS = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_AB_distance <- bind_rows(inversion_AB_distance, inversion_AB_dist_tmp)
  }
}

write_tsv(inversion_AB_distance, "./inversion_JS_asm5_AB_PO_dist.POref.txt")
write_tsv(AB_overlap, "./inversion_asm5_overlapped_AB_JS.txt")

#Plot histogram of distance between inversion breakpoints and the nearest AB boundaries
inversion_AB_distance_longer <- inversion_AB_distance %>%
  pivot_longer(cols = c("front_bin_JS", "rear_bin_JS"), names_to = "front_rear", values_to = "N_bin") 
inversion_AB_distance_longer$N_bin <- abs(inversion_AB_distance_longer$N_bin)
inversion_AB_distance_longer_orthologous <- inversion_AB_distance_longer

inversion_AB_distance_hist <- hist(inversion_AB_distance_longer_orthologous$N_bin, plot = F, breaks = seq(0, 35))

INV_AB_dist_hist <- ggplot(tibble(breaks = inversion_AB_distance_hist$breaks[-1], counts = inversion_AB_distance_hist$counts)) +
  geom_col(aes(x = breaks, y = counts)) + 
  ggtitle("Distance between breakpoint and AB boundary") +
  scale_x_continuous(breaks = c(1, seq(5, 85, 5)), limits = c(0.5, 83)) +
  scale_y_continuous(breaks = seq(0, 20, 2)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1, 
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )

ggsave("INV_AB_distans_bin_JS.eps", INV_AB_dist_hist, width =5, height = 5)

################################################################################
#All inversion
#Make a table which describe the distance between inbersion breakpoints and the nearest AB boundaries.
AB_overlap_all <- tibble()
inversion_AB_distance_all <- tibble()
for(i in 1:nrow(inversion_ordered)){
  #Extract inversion overlapped AB      
  AB_tmp <- bt.intersect(AB, inversion_ordered[i, ], "-wa")

  #If overlapped AB existed      
  if(nrow(AB_tmp) > 0){
    colnames(AB_tmp) <- c("chr", "start", "end", "compartment", "eigenvalue", "note", "size")
    
    AB_tmp <- AB_tmp %>%
      dplyr::mutate(inversionID = pull(inversion_ordered[i, "inversionID"]))
    
    AB_overlap_all <- bind_rows(AB_overlap_all, AB_tmp)

    #Calculate distance between inversion breakpoint and the nearest AB            
    front_distance <- find_closest_to_zero(pull(inversion_ordered[i, "start"]) - unique(c(AB_tmp$start, AB_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_ordered[i, "end"]) - unique(c(AB_tmp$start, AB_tmp$end)))

    #Convert distance from bp to bin size            
    inversion_AB_dist_tmp <- inversion_ordered[i, ] %>%
      dplyr::mutate(front_distance_JS = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance_JS = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_AB_distance_all <- bind_rows(inversion_AB_distance_all, inversion_AB_dist_tmp)
    
  #If overlapped AB did'nt existed
  }else if(nrow(AB_tmp) == 0){
     
  #Extract AB which is on the same chromosome with the focal inversion    
    AB_tmp <- TAD %>%
      dplyr::filter(chr == pull(inversion_ordered[i, "chr"]))

  #Calculate distance between inversion breakpoint and the nearest AB        
    front_distance <- find_closest_to_zero(pull(inversion_ordered[i, "start"]) - unique(c(AB_tmp$start, AB_tmp$end)))
    rear_distance <- find_closest_to_zero(pull(inversion_ordered[i, "end"]) - unique(c(AB_tmp$start, AB_tmp$end)))

  #Convert distance from bp to bin size          
    inversion_AB_dist_tmp <- inversion_ordered[i, ] %>%
      dplyr::mutate(front_distance_JS = front_distance,
                    front_bin_JS = front_distance/25000,
                    rear_distance_JS = rear_distance,
                    rear_bin_JS = rear_distance/25000)
    
    inversion_AB_distance_all <- bind_rows(inversion_AB_distance_all, inversion_AB_dist_tmp)
  }
}

write_tsv(inversion_AB_distance_all, "./inversion_JS_asm5_AB_JS_dist_all.POref.txt")
write_tsv(AB_overlap_all, "./inversion_asm5_overlapped_AB_JS_all.txt")

#Plot histogram of distance between inversion breakpoints and the nearest AB boundaries
inversion_AB_distance_longer_all <- inversion_AB_distance_all %>%
  pivot_longer(cols = c("front_bin_JS", "rear_bin_JS"), names_to = "front_rear", values_to = "N_bin") 
inversion_AB_distance_longer_all$N_bin <- abs(inversion_AB_distance_longer_all$N_bin)

max(inversion_AB_distance_longer_all$N_bin)

inversion_AB_distance_all_hist <- hist(inversion_AB_distance_longer_all$N_bin, plot = F, breaks = seq(0, 73))

INV_AB_dist_all_hist <- ggplot(tibble(breaks = inversion_AB_distance_all_hist$breaks[-1], counts = inversion_AB_distance_all_hist$counts)) +
  geom_col(aes(x = breaks, y = counts)) + 
  ggtitle("Distance between breakpoint and AB boundary") +
  scale_x_continuous(breaks = c(1, seq(5, 160, 5)), limits = c(0.5, 73)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1, 
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )

ggsave("INV_AB_distans_bin_JS_all.POref.eps", INV_AB_dist_all_hist, width =5, height = 5)

