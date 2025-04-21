#Assign compartment status to each gene.
#Dependency: tidyverse v2.0.0, bedtoolsr v.2.30.0-6, bedtools v2.31.1
library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "~/software/bedtools-2.31.1/bedtools2/bin/")

#File path
gff_file <- "~/db/PO_merged_maskXII_XXI_Y_20240517/PO_male_merged.v20240517.mod.maskXII_XXI_Y.gff"
AB_file <- "../iconic_PO_refPO_merged_maskXII_XXI_Y_pairtools_KR_100kb_domains.bed"
mask_file <- "~/db/PO_merged_maskXII_XXI_Y_20240517/mask_regions_PAR_chrXIIrep_chrXXIrep.bed"
AB_settings <- "KR_100kb"

#Import gff
gff_header <-c("chr", "source", "feature", "start", "end", "score", "direction", "frame", "notes")
PO_gff <- read_tsv(gff_file, skip = 1 , col_names = gff_header)  %>%
  dplyr::filter(!(chr %in% c("chrXIX", "chrY")))

#Extract genes
PO_gff_gene <- PO_gff %>% dplyr::filter(feature == "gene") 

#Convert notes column to table
gff_notes <- PO_gff_gene$notes %>% str_split(pattern = ";", simplify = T) 
No_col <- ncol(gff_notes)

for(i in 1:No_col){
  
  assign(paste0("col", i), gff_notes[, i] %>% str_split(pattern = "=", simplify = T))
  assign(paste0("col", i, "val"), get(paste0("col", i))[, 2])
  assign(paste0("colname", i), unique(get(paste0("col", i))[, 1]))
  
}

gff_notes_table <- bind_cols(mget(paste0("col", seq(1, No_col), "val")))

colnames(gff_notes_table) <- unlist(mget(paste0(paste0("colname", seq(1, No_col)))))

#Merge tables
PO_gff_gene_table <- bind_cols(PO_gff_gene[, 1:8], gff_notes_table)

PO_bed_gene_table <- PO_gff_gene_table[, c(1, 4, 5, 9, 13)]
PO_bed_gene_table$start <- PO_bed_gene_table$start - 1 

AB_boundary <- read_tsv(AB_file, col_names = c("chr", "start", "end", "compartment", "eigen"))%>%
  dplyr::mutate(size = end - start + 1)%>%
  dplyr::filter(!(chr %in% c("chrY"))) %>%
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
PO_gff_gene_table_A <- bedtoolsr::bt.intersect(PO_bed_gene_table, A_compartment, "-wa") 
PO_gff_gene_table_A <- bedtoolsr::bt.intersect(PO_gff_gene_table_A, A_compartment_boundary, "-v")
PO_gff_gene_table_A <- bedtoolsr::bt.intersect(PO_gff_gene_table_A, mask_file, "-v")

write_tsv(PO_gff_gene_table_A, paste0("./A_genes_", AB_settings, ".txt"))

PO_gff_gene_table_B <- bedtoolsr::bt.intersect(PO_bed_gene_table, B_compartment)
PO_gff_gene_table_B <- bedtoolsr::bt.intersect(PO_gff_gene_table_B, B_compartment_boundary, "-v")
PO_gff_gene_table_B <- bedtoolsr::bt.intersect(PO_gff_gene_table_B, mask_file, "-v")

write_tsv(PO_gff_gene_table_B, paste0("./B_genes_", AB_settings, ".txt"))

