## -- Figure 2E

go_peroxisome <- read_tsv("../data/meta/GO_terms/GO0005777.tsv")

limma_level %>%
  mutate(is.peroxisomal = Protein.Group %in% go_peroxisome$Protein.Group) -> limma_level_peroxisomal

ggplot()+
  geom_point(data = limma_level_peroxisomal %>% filter(is.peroxisomal == FALSE),
             aes(x = logFC, y = -log10(adj.P.Val)), pch = 21, fill = "grey80", size = 3, alpha = 0.5)+
  geom_point(data = limma_level_peroxisomal %>% filter(is.peroxisomal == TRUE),
             aes(x = logFC, y = -log10(adj.P.Val)), pch = 21, fill = "orange", size = 3)+
  geom_text_repel(data = limma_level_peroxisomal %>% filter(is.peroxisomal == TRUE, -log10(adj.P.Val) > 5),
                  aes(x = logFC, y = -log10(adj.P.Val), label = Genes),
                  nudge_x = 0.2, nudge_y = 0.2, segment.curvature = -0.2, segment.ncp = 1, segment.angle = 90, seed = '1234')+
  theme_classic()

ggsave(file = "../output/Figures/Figure_2E.pdf", width = 10, height = 10)
