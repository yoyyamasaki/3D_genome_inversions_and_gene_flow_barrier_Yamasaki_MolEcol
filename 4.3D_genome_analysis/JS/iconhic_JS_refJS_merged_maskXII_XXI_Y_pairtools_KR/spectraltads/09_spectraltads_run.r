#Run SpectralTAD.
#Dependency: tidyverse v2.0.0, SpectralTAD v1.22.0
library(tidyverse)
library(SpectralTAD)

#Define variables
BIN <- c("25kb", "50kb", "100kb")
RESOLUTION <- c(25000, 50000, 100000)
PREFIX <- "iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_spectraltad"

CHR <- c(paste0("chr", as.roman(seq(1:21))))

#Run spectral TAD in defferent bin size
for(j in 1:length(BIN)){
  
  BIN_tmp <- BIN[j]
  RESOLUTION_tmp <- RESOLUTION[j]
  
  contact <- read_tsv(paste0("../hic/binned/iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_", BIN_tmp, ".txt"), 
                      col_names = c("chr1", "chr1_start", "chr1_end", "chr2", "chr2_start", "chr2_end", "contact")) 
  
  assign(paste0("contact_tad_", BIN_tmp, "_level1"), tibble())
  assign(paste0("contact_tad_", BIN_tmp, "_level2"), tibble())
  
  for(i in CHR){
    contact_chr <- contact %>% dplyr::filter(chr1 == i) %>%
      dplyr::select("chr1_start", "chr2_start", "contact")
    
    assign(paste0("contact_tad_", BIN_tmp, "_", i), SpectralTAD(contact_chr, 
                                                                     chr = i, 
                                                                     resolution = RESOLUTION_tmp, 
                                                                     qual_filter = FALSE, 
                                                                     z_clust = FALSE, 
                                                                     levels = 2
                                                                     )
           )
    
    assign(paste0("contact_tad_", BIN_tmp, "_level1"), dplyr::bind_rows(get(paste0("contact_tad_", BIN_tmp, "_level1")),
                                                             get(paste0("contact_tad_", BIN_tmp, "_", i))[["Level_1"]]
                                                             )
    )
    assign(paste0("contact_tad_", BIN_tmp, "_level2"), dplyr::bind_rows(get(paste0("contact_tad_", BIN_tmp, "_level2")),
                                                                        get(paste0("contact_tad_", BIN_tmp, "_", i))[["Level_2"]]
                                                                        )
    )
  }
  
  write_tsv(get(paste0("contact_tad_", BIN_tmp, "_level1")), paste0(PREFIX, "_", BIN_tmp, "_Level1", ".bed"), col_names = F)
  write_tsv(get(paste0("contact_tad_", BIN_tmp, "_level2")), paste0(PREFIX, "_", BIN_tmp, "_Level2", ".bed"), col_names = F)
}

#Count TAD number
contact_tad_25kb_level1 %>% 
  pivot_longer(cols = c(start, end), names_to = "position", values_to = "bp") %>%
  dplyr::select(chr, bp) %>%
  unique() %>%
  nrow()

contact_tad_25kb_level2 %>% 
  pivot_longer(cols = c(start, end), names_to = "position", values_to = "bp") %>%
  dplyr::select(chr, bp) %>%
  unique() %>%
  nrow()

