# Figure S2A

# The cor function is commented here to speed up processing the whole pipeline

#cor(t(d_stats), method = "pearson", use = "pairwise.complete.obs") -> cor
#save(cor, file = "../output/Variables/Biopsies_correlation.RData")

load("../output/Variables/Biopsies_correlation.RData")

as.data.frame(cor) %>%
  dplyr::select("P01009") %>%
  dplyr::rename(cor = P01009) %>%
  rownames_to_column("protein") %>%
  left_join(meta_pg) %>%
  left_join(pg_n_detected) %>%
  arrange(cor) %>%
  mutate(rank = c(1:nrow(.))) %>%
  mutate(protein = str_replace_all(protein, ";.*", "")) -> cor_alpha1

ggplot(cor_alpha1, aes(x = rank, y = cor)) +
  geom_point()+
  theme_classic()+
  geom_hline(yintercept = 0, lty = "dotted")

ggsave(file = "../output/Figures/Figure_S3A.pdf", width = 5, height = 5)
