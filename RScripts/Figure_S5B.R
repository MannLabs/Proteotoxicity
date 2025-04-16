## -- Figure S5B

gsea_f1 <- limma_level_f %>%
  select(Genes, F1_logFC)

write.table(gsea_f1, file = "../output/Tables/GSEA_F1.rnk", sep = "\t", quote = F, row.names = F, col.names = F)

gsea_f1_result <- read_tsv("../data/meta/GSEA_output/F1_gsea__enrichment_results_wg_result1733327946.txt")

gsea_f4 <- limma_level_f %>%
  select(Genes, F4_logFC)

write.table(gsea_f4, file = "../output/Tables/GSEA_F4.rnk", sep = "\t", quote = F, row.names = F, col.names = F)

gsea_f4_result <- read_tsv("../data/meta/GSEA_output/F4_gsea__enrichment_results_wg_result1733328043.txt")

gsea_f_result <- data.frame(gsea_f1_result %>% select(description, normalizedEnrichmentScore, FDR) %>% rename(NES_F1 = normalizedEnrichmentScore, FDR_F1 = FDR))  %>%
  left_join(gsea_f4_result %>% select(description, normalizedEnrichmentScore, FDR) %>% rename(NES_F4 = normalizedEnrichmentScore, FDR_F4 = FDR)) %>%
  mutate(significant = FDR_F1 < 0.1 | FDR_F4 < 0.1)

ggplot(gsea_f_result, aes(x = NES_F1, y = NES_F4, fill = significant))+
  geom_point(pch = 21)+
  geom_text_repel(data = gsea_f_result %>% filter(significant == TRUE), aes(x = NES_F1, y = NES_F4, label = description))+
  geom_abline(yintercept = 0, slope = 1, lty = "dotted")+
  geom_hline(yintercept = 0, lty = "dotted") +
  geom_vline(xintercept = 0, lty = "dotted") +
  theme_classic()+
  coord_fixed()+
  scale_fill_manual(values = c("white", "black"))

ggsave(file = "../output/Figures/Figure_S5B.pdf", width = 8, height = 8)
