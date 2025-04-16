## -- Figure S7B

stats_biopsies <- read_tsv("../output/Tables/Table-S1-limma.tsv") %>%
  select(Protein.Group, logFC, adj.P.Val) %>%
  dplyr::rename(logFC_biopsies = logFC, adj.P.Val_biopsies = adj.P.Val)

limma_border_firstrow %>%
  left_join(stats_biopsies) %>%
  ggplot(aes(x = logFC, y = logFC_biopsies))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, lty = "dotted")+
  geom_smooth(method = "lm", se = F, lty = "dotted")+
  geom_hline(yintercept = 0, lty = "dotted")+
  geom_vline(xintercept = 0, lty = "dotted")+
  theme_classic()+
  ggrepel::geom_text_repel(data = . %>% filter(adj.P.Val < 0.05),
                           aes(x = logFC, y = logFC_biopsies, label = Genes),
                           nudge_x = 0.2, nudge_y = 0.2, segment.curvature = -0.2, segment.ncp = 1, segment.angle = 90, seed = '1234')+
  scale_y_continuous(breaks = seq(-2,4,by = 1))+
  scale_x_continuous(breaks = seq(-2,4,by = 1))


ggsave("../output/Figures/Figure_S7B.pdf", width = 5, height = 5)
