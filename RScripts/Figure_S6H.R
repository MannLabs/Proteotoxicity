## -- Figure S6H

d_norm %>%
  drop_na(int_norm) %>%
  group_by(Protein.Group) %>%
  summarise(n = n()) -> pg_n_detected

pg_n_detected %>%
  filter(n > 0.3*length(SA_incl)) %>%
  pull(Protein.Group) -> pg_30

d_norm %>%
  filter(Protein.Group %in% pg_30) %>%
  filter(biology == "Border", spatial %in% c("plusone", "minusone")) %>%
  mutate(well = paste(batch, well, sep = "_")) %>%
  dplyr::select(well, Protein.Group, int_norm) %>%
  spread(well, int_norm) %>%
  column_to_rownames("Protein.Group") -> d_borderstats

condition <- as.factor((meta %>%
                          mutate(well = paste(batch, well, sep = "_")) %>%
                          column_to_rownames("well"))[colnames(d_borderstats),]$spatial)

sample <- as.factor((meta %>%
                       mutate(well = paste(batch, well, sep = "_")) %>%
                       column_to_rownames("well"))[colnames(d_borderstats),]$sample_short)

design <- model.matrix(~ condition + sample)

fit <- limma::lmFit(d_borderstats, design)
fit <- limma::eBayes(fit)
limma_border_firstrow <- limma::topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein.Group") %>%
  left_join(meta_pg)

ggplot(data = limma_border_firstrow , aes(x = logFC, y = -log10(adj.P.Val), fill = -log10(adj.P.Val)))+
  geom_point(pch = 21, color = "black", size = 3)+
  #geom_hline(yintercept = -log10(0.05), lty = "dotted")+
  ggrepel::geom_text_repel(data = limma_border_firstrow %>% filter(adj.P.Val < 0.05),
                           aes(x = logFC, y = -log10(adj.P.Val), label = Genes),
                           nudge_x = 0.2, nudge_y = 0.2, segment.curvature = -0.2, segment.ncp = 1, segment.angle = 90, seed = '1234')+
  theme_classic()+
  viridis::scale_fill_viridis(option = "inferno") -> plot_volcano_border

ggsave(plot_volcano_border, file = "../output/Figures/Figure_S6H.pdf", width = 10, height = 8)
