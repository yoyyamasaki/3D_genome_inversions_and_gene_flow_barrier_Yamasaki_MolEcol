#Assign compartment status to each gene.
#Dependency: tidyverse v2.0.0, bedtoolsr v.2.30.0-6
library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/software/bedtools-2.31.1/bedtools2/bin/")

#File path
gff_file <- "~/db/JS_male_merged_mask_chrXII_chrXXI_chrYXI_20240517/JS_male_merged.v20240517.mod.maskXII_XXI_Y.gff"
AB_file <- "../iconhic_JS_refJS_merged_maskXII_XXI_Y_pairtools_KR_100kb_domains.bed"
mask_file <- "~/db/JS_male_merged_mask_chrXII_chrXXI_chrYXI_20240517/mask_chrXII_chrXXI_chrY_IX.bed"
AB_settings <- "KR_100kb"

#Import gff
gff_header <-c("chr", "source", "feature", "start", "end", "score", "direction", "frame", "notes")
JS_gff <- read_tsv(gff_file, skip = 1 , col_names = gff_header) %>%
  dplyr::filter(!(chr %in% c("chrY_IX")))

#Extract genes
JS_gff_gene <- JS_gff %>% 
  dplyr::filter(feature == "gene") 

#Convert notes column to table
gff_notes <- JS_gff_gene$notes %>% str_split(pattern = ";", simplify = T) 
No_col <- ncol(gff_notes)

for(i in 1:No_col){

  assign(paste0("col", i), gff_notes[, i] %>% str_split(pattern = "=", simplify = T))
  assign(paste0("col", i, "val"), get(paste0("col", i))[, 2])
  assign(paste0("colname", i), unique(get(paste0("col", i))[, 1]))
  
}

gff_notes_table <- bind_cols(mget(paste0("col", seq(1, No_col), "val")))

colnames(gff_notes_table) <- unlist(mget(paste0(paste0("colname", seq(1, No_col)))))

#Merge tables
JS_gff_gene_table <- bind_cols(JS_gff_gene[, 1:8], gff_notes_table)

JS_bed_gene_table <- JS_gff_gene_table[, c(1, 4, 5, 9, 13)]
JS_bed_gene_table$start <- JS_bed_gene_table$start - 1 

AB_boundary <- read_tsv(AB_file, col_names = c("chr", "start", "end", "compartment", "eigen")) %>%
  dplyr::mutate(size = end - start + 1)%>%
  dplyr::filter(!(chr %in% c("chrIX_Y"))) %>%
  dplyr::filter(eigen != 0)

A_compartment <- AB_boundary %>% dplyr::filter(compartment == "A")
A_compartment_boundary <- A_compartment %>% 
  pivot_longer(cols= c("start", "end"), names_to = "position", values_to = "start") %>%
  dplyr::mutate(end = start + 1) %>%
  dplyr::select(chr, start, end)

B_compartment <- AB_boundary %>% dplyr::filter(compartment == "B")
B_compartment_boundary <- B_compartment %>% 
  pivot_longer(cols= c("start", "end"), names_to = "position", values_to = "start") %>%
  dplyr::mutate(end = start + 1) %>%
  dplyr::select(chr, start, end)

#Extract genes in each compartment
JS_gff_gene_table_A <- bedtoolsr::bt.intersect(JS_bed_gene_table, A_compartment, "-wa") 
JS_gff_gene_table_A <- bedtoolsr::bt.intersect(JS_gff_gene_table_A, A_compartment_boundary, "-v")
JS_gff_gene_table_A <- bedtoolsr::bt.intersect(JS_gff_gene_table_A, mask_file, "-v")

write_tsv(JS_gff_gene_table_A, paste0("./A_genes_", AB_settings, ".txt"))

JS_gff_gene_table_B <- bedtoolsr::bt.intersect(JS_bed_gene_table, B_compartment)
JS_gff_gene_table_B <- bedtoolsr::bt.intersect(JS_gff_gene_table_B, B_compartment_boundary, "-v")
JS_gff_gene_table_B <- bedtoolsr::bt.intersect(JS_gff_gene_table_B, mask_file, "-v")

write_tsv(JS_gff_gene_table_B, paste0("./B_genes_", AB_settings, ".txt"))

