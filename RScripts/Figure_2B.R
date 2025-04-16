## -- Figure 2B

SA_stats <- meta_biopsies %>%
  pull(well_id)

d_long %>%
  filter(Protein.Group == "P01009") %>%
  drop_na(int) %>%
  arrange(int) %>%
  pull(well_id) -> well_serpina1

pg_significant <- limma_level %>%
  filter(adj.P.Val < 0.01) %>%
  pull(protein)

d_long %>%
  filter(well_id %in% SA_stats) %>%
  filter(Protein.Group %in% pg_30) %>%
  mutate(int = log2(int)) %>%
  dplyr::select(- ms_id) %>%
  spread(well_id, int) %>%
  column_to_rownames("Protein.Group") -> d_stats

myBreaks <- c(seq(-2,2, by = 0.1))
myColor <- colorRampPalette(viridis(100, option = "magma"))(length(myBreaks))

pdf("../output/Figures/Figure_2B.pdf")
pheatmap::pheatmap(t(scale(t(d_stats[pg_significant,well_serpina1]))), show_rownames = F, show_colnames = F,
                   cutree_rows = 7,
                   breaks = myBreaks,
                   color = myColor, cluster_cols = F)
dev.off()
