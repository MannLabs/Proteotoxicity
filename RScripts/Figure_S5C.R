## -- Figure S5

fibrosis_colors = c("#4150a2", "#5ab4e5", "#7699a8", "#000000")

gsea_f_result <- data.frame(gsea_f1_result %>% select(geneSet, description, normalizedEnrichmentScore, FDR) %>% rename(NES_F1 = normalizedEnrichmentScore, FDR_F1 = FDR))  %>%
  left_join(gsea_f4_result %>% select(description, normalizedEnrichmentScore, FDR) %>% rename(NES_F4 = normalizedEnrichmentScore, FDR_F4 = FDR)) %>%
  mutate(significant = FDR_F1 < 0.1 | FDR_F4 < 0.1) %>%
  mutate(pathways_id = paste(geneSet, description))

d_long %>%
  mutate(int = log2(int)) %>%
  left_join(meta) %>%
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
  left_join(meta) %>%
  left_join(level_of_alpha1) %>%
  drop_na(z_score) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                         keys=Genes, 
                                         column="ENSEMBL", 
                                         keytype="SYMBOL",
                                         multiVals="first")) -> d_group_graph_kleiner

pathways_to_map_extended <- gsea_f_result %>% filter(significant == TRUE) %>% pull(pathways_id)

for(i in pathways_to_map_extended){
  
  db_kegg_subset <- data.frame(ENTREZ = str_replace(db_kegg$kg.sets[[i]], ".*_", "")) %>%
    mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                           keys=ENTREZ, 
                                           column="ENSEMBL", 
                                           keytype="ENTREZID",
                                           multiVals="first"))
  
  ggplot() +
    geom_smooth(data = d_group_graph_kleiner %>% filter(ENSEMBL %in% db_kegg_subset$ENSEMBL),
                aes(y = z_score, x = log2(ref_int), color =  as.factor(kleiner_score)), method = "loess", se = FALSE, linewidth = 2) +
    scale_x_continuous(breaks = seq(20,30, by = 1)) +
    labs(x = "Intensity Alpha-1", y = "POI") +
    theme_classic() +
    scale_color_manual(values = fibrosis_colors)+
    labs(title = i)
  
  ggsave(paste("../output/Figures/Figure_S5C_", i, ".pdf", sep = ""), width = 6, height = 4)
}
