## -- Figure 2G

d_long %>%
  mutate(int = log2(int)) %>%
  filter(Protein.Group %in% pg_30) %>%
  left_join(meta_biopsies) %>%
  drop_na(kleiner_score) %>%
  mutate(Protein.Group_kleiner = paste(Protein.Group, kleiner_score, sep = "__")) %>%
  dplyr::select(Protein.Group_kleiner, ms_id, int) %>%
  spread(ms_id, int) %>%
  column_to_rownames("Protein.Group_kleiner")  -> d_wide_kleiner

as.data.frame(t(scale(t(d_wide_kleiner)))) %>%
  rownames_to_column("Protein.Group_kleiner") %>%
  gather(ms_id, z_score, !Protein.Group_kleiner) %>%
  mutate(Protein.Group = str_replace(Protein.Group_kleiner, "__.*", ""), kleiner_score = as.numeric(as.character(str_replace(Protein.Group_kleiner, ".*__", "")))) %>%
  left_join(d_long) %>%
  left_join(meta_pg) %>%
  left_join(meta_biopsies) %>%
  left_join(level_of_alpha1) %>%
  drop_na(z_score) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                         keys=Genes, 
                                         column="ENSEMBL", 
                                         keytype="SYMBOL",
                                         multiVals="first")) -> d_group_graph_kleiner

db_kegg_specific <- data.frame(ENTREZ = str_replace(db_kegg$kg.sets[["hsa04146 Peroxisome"]], ".*_", "")) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                          keys=ENTREZ, 
                          column="ENSEMBL", 
                          keytype="ENTREZID",
                          multiVals="first"))

d_scaled %>%
  filter(ENSEMBL %in% db_kegg_specific$ENSEMBL) %>%
  drop_na(z_score) -> tmp

ggplot() +
  #facet_grid(.~ Protein.Group)+
  #coord_fixed(ratio = 1) +
  geom_smooth(data = tmp %>% filter(kleiner_score %in% 1),
              aes(y = z_score, x = log2(ref_int)), method = "loess", se = FALSE, color = "purple", linewidth = 2) +
  geom_smooth(data = tmp %>% filter(kleiner_score %in% 2),
              aes(y = z_score, x = log2(ref_int)), method = "loess", se = FALSE, color = "grey50", linewidth = 2) +
  geom_smooth(data = tmp %>% filter(kleiner_score %in% 3),
              aes(y = z_score, x = log2(ref_int)), method = "loess", se = FALSE, color = "orange", linewidth = 2) +
  geom_smooth(data = tmp %>% filter(kleiner_score %in% 4),
              aes(y = z_score, x = log2(ref_int)), method = "loess", se = FALSE, color = "darkred", linewidth = 2) +
  scale_x_continuous(breaks = seq(20,30, by = 1)) +
  theme_classic()+
  theme(legend.position="none") +
  labs(title = "hsa04146 Peroxisome") -> plot

ggsave(plot, file = "../output/Figures/Figure_2G.pdf")
