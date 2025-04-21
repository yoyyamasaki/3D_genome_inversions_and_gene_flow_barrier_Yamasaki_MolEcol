#Plot Fst and dxy, and conduct statistical test.

library(tidyverse)
library(ggsignif)
library(lawstat)
library(brunnermunzel)


#Fst
#Read files
inINV_fst <- read_tsv("../PO_JS_Quebec_WLD_bcftools.inINV.filtered.POmaskXIIXXIYv20240517_fst.txt") %>%
  dplyr::filter(!(chromosome %in% c("chrIX", "chrXIX", "chrY"))) %>%
  dplyr::mutate(inversion = "in")

outINV_fst <- read_tsv("../PO_JS_Quebec_WLD_bcftools.SYNAL.filtered.POmaskXIIXXIYv20240517_fst.txt") %>%
  dplyr::filter(!(chromosome %in% c("chrIX", "chrXIX", "chrY"))) %>%
  dplyr::mutate(inversion = "out")

INV_fst <- bind_rows(inINV_fst, outINV_fst) %>%
  dplyr::filter(no_snps > 25) %>%
  dplyr::filter(pop1 == "JS" & pop2 == "PO")

#Plot box plot
plot_box_fst <- ggplot(data = INV_fst, aes(x = inversion, y = avg_wc_fst)) +
  geom_violin(aes(color = inversion), width = 1, linewidth = 0.3) +
  geom_boxplot(aes(fill = inversion), width = 0.2, color = "black", linewidth = 1) +
  geom_signif(comparisons = list(c("in", "out")),
              test = "brunner.munzel.test", na.rm = FALSE, map_signif_level = F, col = "red") +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_fill_manual(values = c("red3", "cyan")) +
  ggtitle("Fst inversion") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1, 
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )
ggsave("inversion_fst_JS_PO.eps", plot_box_fst, height = 5, width = 5)

#Calculate mean and median
drop_na(INV_fst) %>% 
  dplyr::group_by(inversion) %>%
  dplyr::summarize(mean = mean(avg_wc_fst))

drop_na(INV_fst) %>% 
  dplyr::group_by(inversion) %>%
  dplyr::summarize(median = median(avg_wc_fst))


#dxy
#Read files
inINV_dxy <- read_tsv("../PO_JS_Quebec_WLD_bcftools.inINV.filtered.POmaskXIIXXIYv20240517_dxy.txt") %>%
  dplyr::filter(!(chromosome %in% c("chrIX", "chrXIX", "chrY"))) %>%
  dplyr::mutate(inversion = "in")

outINV_dxy <- read_tsv("../PO_JS_Quebec_WLD_bcftools.SYNAL.filtered.POmaskXIIXXIYv20240517_dxy.txt") %>%
  dplyr::filter(!(chromosome %in% c("chrIX", "chrXIX", "chrY"))) %>%
  dplyr::mutate(inversion = "out")

INV_dxy <- bind_rows(inINV_dxy, outINV_dxy) %>%
  dplyr::filter(no_sites > 2500) %>%
  dplyr::filter(pop1 == "JS" & pop2 == "PO")

max(INV_dxy$avg_dxy)

#Plot box plot
plot_box_dxy <- ggplot(data = INV_dxy, aes(x = inversion, y = avg_dxy)) +
  geom_violin(aes(color = inversion), width = 1, linewidth = 0.3) +
  geom_boxplot(aes(fill = inversion), width = 0.2, color = "black", linewidth = 1) +
  geom_signif(comparisons = list(c("in", "out")),
              test = "brunner.munzel.test", na.rm = FALSE, map_signif_level = F, col = "red") +
  scale_y_continuous(limits = c(0, 0.035)) +
  scale_fill_manual(values = c("red3", "cyan")) +
  ggtitle("dxy inversion") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 1,
                                    linetype = NULL, color = NULL, inherit.blank = T),
        legend.position = "none"
  )
ggsave("inversion_dxy_JS_PO.eps", plot_box_dxy, height = 5, width = 5)

#Calculate mean and median
drop_na(INV_dxy) %>% 
  dplyr::group_by(inversion) %>%
  dplyr::summarize(mean = mean(avg_dxy))

drop_na(INV_dxy) %>% 
  dplyr::group_by(inversion) %>%
  dplyr::summarize(median = median(avg_dxy))

