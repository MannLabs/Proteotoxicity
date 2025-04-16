## -- Figure S2E

plasma <- read_tsv("../data/meta/HPA_plasma-proteome-liver.txt")

limma_level %>%
  mutate(plasma = Genes %in% plasma$Gene) -> limma_plasma

  ggplot()+
    geom_boxplot(data = limma_plasma, aes(x = plasma, y = logFC), outlier.shape = NA) +
    geom_jitter(data = limma_plasma %>% filter(group == "not significant"), aes(x = plasma, y = logFC), color = "grey80", width = 0.1, alpha = 0.5) +
    geom_jitter(data = limma_plasma %>% filter(group != "not significant"), aes(x = plasma, y = logFC), color = "black", width = 0.1, alpha = 0.5)+
    theme_classic()+
    geom_hline(yintercept = 0, lty = "dotted") -> plot

ggsave(plot, file = "../output/Figures/Figure_S2E.pdf", width = 4, height = 5)

table((limma_level %>%
        mutate(plasma = Genes %in% plasma$Gene))$plasma)

table((limma_level %>%
        mutate(plasma = Genes %in% plasma$Gene)%>%
        filter(plasma == TRUE))$group)
