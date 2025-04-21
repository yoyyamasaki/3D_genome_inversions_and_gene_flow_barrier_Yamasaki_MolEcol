#Make inversion bed files. They were based on G. aculeatus and G. nipponicus genome coordinate. Also make 25kb < inversion bed files.
#Dependency: tidyverse v2.0.0
library(tidyverse)

#Read file
inversion_PO_JS <- read_tsv("../syri20240917_PO_JS_asm5syri.inversion.bedpe", col_names = c("chr_PO", "start_PO", "end_PO", "chr_JS", "start_JS", "end_JS")) %>%
  dplyr::mutate(length_PO = end_PO - start_PO, length_JS = end_JS - start_JS)

CHR <- paste0("chr", as.roman(seq(1, 21)))

#Make inversion bed file for G. aculeatus genome corrdinate
inversion_PO <- inversion_PO_JS %>%
  dplyr::select(chr_PO, start_PO, end_PO)

#Make 25kb < inversion bed file for G. aculeatus genome corrdinate
inversion_PO_25kb <- inversion_PO_JS %>%
  dplyr::filter(length_PO >= 25000) %>%
  dplyr::select(chr_PO, start_PO, end_PO)

#Make inversion bed file for G. nipponicus genome corrdinate
inversion_JS <- inversion_PO_JS %>%
  dplyr::select(chr_JS, start_JS, end_JS)

#Make 25kb < inversion bed file for G. nipponicus genome corrdinate
inversion_JS_25kb <- inversion_PO_JS %>%
  dplyr::filter(length_PO >= 25000) %>%
  dplyr::select(chr_JS, start_JS, end_JS)

#Write diles
write_tsv(inversion_PO_25kb, "./inversion_PO_25kb_PObased.bed", col_names = F)
write_tsv(inversion_JS_25kb, "./inversion_JS_25kb_PObased.bed", col_names = F)

write_tsv(inversion_PO, "./inversion_PO_PObased.bed", col_names = F)
write_tsv(inversion_JS, "./inversion_JS_PObased.bed", col_names = F)
