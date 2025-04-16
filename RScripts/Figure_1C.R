## -- Figure 1C

d_long %>%
  drop_na(int) %>%
  group_by(Protein.Group) %>%
  summarise(n = n()) -> pg_n_detected

pg_n_detected %>%
  filter(n > 0.3*length(SA_incl)) %>%
  pull(Protein.Group) -> pg_30

# Volcano plot
SA_stats <- meta_biopsies %>%
  pull(well_id)

d_long %>%
  filter(Protein.Group %in% pg_30) %>%
  filter(well_id %in% SA_stats) %>%
  mutate(int = log2(int)) %>%
  select(- ms_id) %>%
  spread(well_id, int) %>%
  column_to_rownames("Protein.Group") -> d_stats

patient <- as.factor((meta_biopsies %>% column_to_rownames("well_id"))[colnames(d_stats),]$sample)
level <- as.factor((meta_biopsies %>% column_to_rownames("well_id"))[colnames(d_stats),]$level)

design <- model.matrix(~ level + patient)

fit <- lmFit(d_stats, design)
fit <- eBayes(fit)
limma_level <- topTable(fit, number = Inf, confint = TRUE, coef = "level2_high", adjust.method = "fdr") %>%
  rownames_to_column("Protein.Group") %>%
  left_join(meta_pg) %>%
  mutate(group =  ifelse(adj.P.Val > 0.05, "not significant", ifelse(logFC > 0, "upregulated", "downregulated")))

(
  ggplot(data = limma_level , aes(x = logFC, y = -log10(adj.P.Val), fill = -log10(adj.P.Val)))+
    geom_point(pch = 21, color = "black", size = 3)+
    geom_hline(yintercept = -log10(0.05), lty = "dotted")+
    geom_text_repel(data = limma_level %>% filter(adj.P.Val < 0.05),
                    aes(x = logFC, y = -log10(adj.P.Val), label = Genes),
                    nudge_x = 0.2, nudge_y = 0.2, segment.curvature = -0.2, segment.ncp = 1, segment.angle = 90, seed = '1234')+
    theme_bw()+
    scale_fill_viridis(option = "inferno")
) %>% ggsave(., file = "../output/Figures/Figure_1C.pdf", width = 10, height = 8)

write_tsv(limma_level, "../output/Tables/Table-S1-limma.tsv")

# Numbers in the paper
limma_level %>%
  filter(Genes == "SERPINA1") %>%
  mutate(elevation = 2^logFC) %>%
  pull(elevation)

limma_level %>%
  filter(Genes == "HSPA5") %>%
  mutate(elevation = 2^logFC) %>%
  pull(elevation)

limma_level %>%
  filter(Genes == "LMAN1") %>%
  mutate(elevation = 2^logFC) %>%
  pull(elevation)

limma_level %>%
  filter(adj.P.Val < 0.05) %>%
  summarise(n = n()/nrow(limma_level))
