#Plot expression level for each compartment and conduct statistical test.
#Dependency: tidyverse v2.0.0, ggsignif v.0.6.4, brunnermunzel v2.0
library(tidyverse)
library(ggsignif)
library(brunnermunzel)

#File path
tpm_data <- "~/HDD/yo_analysis/HDD4_20230223/stickleback/HiC/chromatin_analysis/JS/RNAseq/JS_male_maskXII_XXI_Y/RNAseq_analysis/counts_tidy_tpm_filtered.txt"
A_compartment_table <- "./A_genes_KR_100kb.txt"
B_compartment_table <- "./B_genes_KR_100kb.txt"
bed_col <- c("chr", "start", "end", "Geneid", "symbol")
CHR <- c(paste0("chr", as.roman(1:21)))

#Read gene names for each compartment
JS_gff_gene_table_A <- read_tsv(A_compartment_table, skip = 1, col_names = bed_col)%>%
  dplyr::filter(!(chr %in% c("chrIX", "chrXIX", "chrY_IX")))
JS_gff_gene_table_B <- read_tsv(B_compartment_table, skip = 1, col_names = bed_col)%>%
  dplyr::filter(!(chr %in% c("chrIX", "chrXIX", "chrY_IX")))

overlapped_genes <- intersect(JS_gff_gene_table_A$Geneid, JS_gff_gene_table_B$Geneid)

#Read TPM data
counts_tidy_tpm_filtered <- read_tsv(tpm_data)

#Merge compartment and TPM data
A_tpm <- left_join(JS_gff_gene_table_A, counts_tidy_tpm_filtered)
A_tpm <- A_tpm %>% dplyr::mutate(compartment = rep("A", nrow(A_tpm)))

B_tpm <- left_join(JS_gff_gene_table_B, counts_tidy_tpm_filtered)
B_tpm <- B_tpm %>% dplyr::mutate(compartment = rep("B", nrow(B_tpm)))

compartment_tpm <- bind_rows(A_tpm, B_tpm) %>%
  drop_na() %>%
  dplyr::filter(! Geneid %in% overlapped_genes) 

#Brunner-Munzel test
brunnermunzel.test(dplyr::filter(compartment_tpm, tissue == "brain") %>% dplyr::filter(compartment == "A") %>% pull(tpm_mean),
                   dplyr::filter(compartment_tpm, tissue == "brain") %>% dplyr::filter(compartment == "B") %>% pull(tpm_mean))

brunnermunzel.test(dplyr::filter(compartment_tpm, tissue == "gill") %>% dplyr::filter(compartment == "A") %>% pull(tpm_mean),
                   dplyr::filter(compartment_tpm, tissue == "gill") %>% dplyr::filter(compartment == "B") %>% pull(tpm_mean))

brunnermunzel.test(dplyr::filter(compartment_tpm, tissue == "liver") %>% dplyr::filter(compartment == "A") %>% pull(tpm_mean),
                   dplyr::filter(compartment_tpm, tissue == "liver") %>% dplyr::filter(compartment == "B") %>% pull(tpm_mean))

brunnermunzel.test(dplyr::filter(compartment_tpm, tissue == "pecmusAD") %>% dplyr::filter(compartment == "A") %>% pull(tpm_mean),
                   dplyr::filter(compartment_tpm, tissue == "pecmusAD") %>% dplyr::filter(compartment == "B") %>% pull(tpm_mean))

brunnermunzel.test(dplyr::filter(compartment_tpm, tissue == "pit") %>% dplyr::filter(compartment == "A") %>% pull(tpm_mean),
                   dplyr::filter(compartment_tpm, tissue == "pit") %>% dplyr::filter(compartment == "B") %>% pull(tpm_mean))

#Plot expression of genes in each compartment in box plot
tissue <- unique(compartment_tpm$tissue)

for(j in tissue){
  dir.create(j)
  compartment_tpm_tissue <- compartment_tpm %>% dplyr::filter(tissue == j) %>%
    dplyr::select(chr, Geneid, symbol, tissue, tpm_mean, compartment) %>%
    unique()
  
  #All chromosome merged
  plot_AB_expression <- ggplot(compartment_tpm_tissue, aes(x = compartment, y = tpm_mean + 1))+
    geom_violin(aes(color = compartment), width = 1, linewidth = 0.3) +
    geom_boxplot(aes(fill = compartment), width = 0.1, color = "black", linewidth = 1) +
    geom_signif(comparisons = list(c("A", "B")),
                test = "brunner.munzel.test", na.rm = FALSE, map_signif_level = F, col = "red") +
    scale_y_log10(limits = c(1, max(compartment_tpm_tissue$tpm_mean) + 100000)) +
    scale_fill_manual(values = c("red3", "cyan")) +
    ggtitle(paste0(j, "_all")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 1,
                                      linetype = NULL, color = NULL, inherit.blank = T),
          legend.position = "none"
    )
  ggsave(paste0(j, "/expression_", j, "_compartment_all.png"), plot_AB_expression, height = 3, width = 3)

  #Per chromosome in one figure
  plot_AB_expression <- ggplot(compartment_tpm_tissue, aes(x = compartment, y = tpm_mean + 1))+
    geom_violin(aes(color = compartment), width = 1, linewidth = 0.3) +
    geom_boxplot(aes(fill = compartment), width = 0.1, color = "black", linewidth = 1) +
    geom_signif(comparisons = list(c("A", "B")),
                test = "brunner.munzel.test", na.rm = FALSE, map_signif_level = F, col = "red") +
    facet_wrap( ~ chr, ncol = 4)+
    scale_y_log10(limits = c(1, max(compartment_tpm_tissue$tpm_mean) + 100000)) +
    scale_fill_manual(values = c("red3", "cyan")) +
    ggtitle(paste0(j, "_all")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", size = 1,
                                      linetype = NULL, color = NULL, inherit.blank = T),
          legend.position = "none"
    )
  ggsave(paste0(j, "/expression_", j, "_compartment_facet.png"), plot_AB_expression, height = 18, width = 12)
  
  #Per chromosome and one plot in one file
  for(i in CHR){
    plot_AB_expression <- ggplot(dplyr::filter(compartment_tpm_tissue, chr == i), aes(x = compartment, y = tpm_mean + 1))+
      geom_violin(aes(color = compartment), width = 1, linewidth = 0.3) +
      geom_boxplot(aes(fill = compartment), width = 0.1, color = "black", linewidth = 1) +
      geom_signif(comparisons = list(c("A", "B")),
                  test = "brunner.munzel.test", na.rm = FALSE, map_signif_level = F, col = "red") +
      scale_y_log10(limits = c(1, max(compartment_tpm_tissue$tpm_mean) + 100000)) +
      scale_fill_manual(values = c("red3", "cyan")) +
      ggtitle(paste0(j, "_", i)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black", size = 1,
                                        linetype = NULL, color = NULL, inherit.blank = T),
            legend.position = "none"
            )
    ggsave(paste0(j, "/expression_", j, "_compartment_", i, ".png"), plot_AB_expression, height = 3, width = 3)
  }
}
