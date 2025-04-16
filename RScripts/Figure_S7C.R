## -- Figure S7C

limma_border_firstrow %>%
  left_join(stats_biopsies) %>%
  ggplot(aes(x = -log10(adj.P.Val), y = -log10(adj.P.Val_biopsies)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), lty = "dotted")+
  geom_vline(xintercept = -log10(0.05), lty = "dotted")+
  theme_classic()+
  ggrepel::geom_text_repel(data = . %>% filter(adj.P.Val < 0.05),
                           aes(x = -log10(adj.P.Val), y = -log10(adj.P.Val_biopsies), label = Genes),
                           nudge_x = 0.2, nudge_y = 0.2, segment.curvature = -0.2, segment.ncp = 1, segment.angle = 90, seed = '1234')+
  scale_y_continuous(breaks = seq(0,30,by = 5))+
  scale_x_continuous(breaks = seq(0,30,by = 5))
  
ggsave("../output/Figures/Figure_S7C.pdf", width = 5, height = 5)
