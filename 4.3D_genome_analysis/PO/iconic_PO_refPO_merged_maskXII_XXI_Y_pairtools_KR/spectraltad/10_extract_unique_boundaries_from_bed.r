#Caluculate unique TAD boundary number and TAD size.
#Dependency: tidyverse v2.0.0
library(tidyverse)
options(scipen = 100)

#Import TAD data
level1 <- read_tsv("./iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level1.bed", col_names = c("chr", "start", "end", "level")) %>%
  pivot_longer(cols = c("start", "end"), names_to = "position", values_to = "bp")  %>%
  unique()

level2 <- read_tsv("./iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level2_unique.bed", col_names = c("chr", "start", "end", "level"))%>%
  pivot_longer(cols = c("start", "end"), names_to = "position", values_to = "bp") %>%
  unique()

CHR <- paste0("chr", as.roman(1:21))

#Extract unique boundaries in each TAD level
level1_boundary <- tibble()
level2_boundary <- tibble()
all_uniq_boundary <- tibble()

for(i in CHR){
  
  level1_tmp <- level1 %>%
    dplyr::filter(chr == i)
  
  level2_tmp <- level2 %>%
    dplyr::filter(chr == i)
  
  level1_tmp_end <- level1_tmp %>%
    dplyr::filter(position == "end")
  
  level2_tmp_end <- level2_tmp %>%
    dplyr::filter(position == "end")
  
  #level1 boundary
  level1_uniqe_tmp <- level1_tmp_end[!(level1_tmp_end$bp %in% level2_tmp_end$bp), ]
  level2_uniqe_tmp <- level2_tmp_end[!(level2_tmp_end$bp %in% level1_tmp_end$bp), ]
  unique_tmp <- bind_rows(level1_uniqe_tmp, level2_uniqe_tmp)
  
  level1_boundary <- bind_rows(level1_boundary, unique_tmp)
  
  #level2 boundary
  shared_tmp <- level1_tmp_end[level1_tmp_end$bp %in% level2_tmp_end$bp, ]
  
  level2_boundary <- bind_rows(level2_boundary, shared_tmp)

}

#all unique boundary
uniq_boundary <- bind_rows(level1, level2) %>%
  dplyr::select(chr, bp) %>%
  unique()

#Caluculate mean and median TAD size
mean(read_tsv("./iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level1.bed", col_names = c("chr", "start", "end", "level")) %>%
       dplyr::mutate(length = end - start) %>%
       pull(length))
median(read_tsv("./iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level1.bed", col_names = c("chr", "start", "end", "level")) %>%
       dplyr::mutate(length = end - start) %>%
       pull(length))

mean(read_tsv("./iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level2_unique.bed", col_names = c("chr", "start", "end", "level")) %>%
       dplyr::mutate(length = end - start) %>%
       pull(length))
median(read_tsv("./iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_spectraltad_25kb_Level2_unique.bed", col_names = c("chr", "start", "end", "level")) %>%
         dplyr::mutate(length = end - start) %>%
         pull(length))

