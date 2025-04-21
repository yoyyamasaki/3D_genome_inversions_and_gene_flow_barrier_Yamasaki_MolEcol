#Conduct TPM normalization.
#Dependency: tidyverse v2.0.0

library(tidyverse)

#Read count data and reformat
dat_raw <- read_tsv("../PO_count.PO_male_merged.v20240517.mod.maskXII_XXI_Y.txt", skip = 1)[, c(-2, -3, -4, -5)]

colnames(dat_raw) <- colnames(dat_raw) %>%
  str_replace_all(pattern = "./mapping/", replacement = "") %>%
  str_replace_all(pattern = "Aligned.sortedByCoord.out.bam", replacement = "")

counts_tidy <- tidyr::pivot_longer(dat_raw, cols = !c(1, 2), names_to = "ID", values_to = "reads")
sample_vec <- unlist(str_split(counts_tidy$ID, pattern ="_"))

counts_tidy <- bind_cols(counts_tidy, 
                         tissue = sample_vec[seq(1, length(sample_vec), 3)],
                         sample = sample_vec[seq(2, length(sample_vec), 3)],
                         sex = sample_vec[seq(3, length(sample_vec), 3)])

#Calculation of TPM
counts_tidy_tpm <- dplyr::mutate(counts_tidy, reads_1K = reads * (1000 / Length))

readsum_tmp <- c()
for(i in unique(counts_tidy_tpm$ID)){
  extracted_tmp <- counts_tidy_tpm %>% dplyr::filter(ID == i)
  readsum_tmp <- c(readsum_tmp, sum(extracted_tmp$reads))
}

total_reads <- tibble(ID = unique(counts_tidy_tpm$ID), readsum = readsum_tmp)
counts_tidy_tpm <-left_join(counts_tidy_tpm, total_reads, by = "ID") %>%
  dplyr::mutate(tpm = reads_1K * (10^6 / readsum))

#Calculate mean TPM per genes and tissues
counts_tidy_tpm <- counts_tidy_tpm %>% 
  dplyr::group_by(Geneid, tissue)
mean_tpm <- summarize(counts_tidy_tpm, tpm_mean = mean(tpm))

counts_tidy_tpm_filtered <- left_join(counts_tidy_tpm, mean_tpm, by = c("Geneid" = "Geneid", "tissue" = "tissue")) 

#Write output file
write_tsv(counts_tidy_tpm_filtered, "./counts_tidy_tpm_filtered.txt")
