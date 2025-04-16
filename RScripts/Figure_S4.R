## -- Figure S4

as.data.frame(t(scale(t(d_wide)))) %>%
  rownames_to_column("Protein.Group") %>%
  gather(ms_id, z_score, !Protein.Group) %>%
  left_join(d_long) %>%
  left_join(meta_pg) %>%
  left_join(meta_biopsies) %>%
  filter(Protein.Group %in% pg_70) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                          keys=Genes, 
                          column="ENSEMBL", 
                          keytype="SYMBOL",
                          multiVals="first")) %>%
  left_join(level_of_alpha1) %>%
  drop_na(ref_int) %>%
  filter(ENSEMBL != "NANA") -> d_scaled

for(i in pathways_to_map){
  
  db_kegg_specific <- data.frame(FLYBASECG = str_replace(db_kegg$kg.sets[[i]], ".*_", "")) %>%
    mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                            keys=FLYBASECG, 
                            column="ENSEMBL", 
                            keytype="ENTREZID",
                            multiVals="first"))
  
  d_scaled %>%
    filter(ENSEMBL %in% db_kegg_specific$ENSEMBL) %>%
    drop_na(z_score) -> tmp
  
  ggplot() +
    #facet_grid(.~ Protein.Group)+
    #coord_fixed(ratio = 1) +
    geom_smooth(data = tmp,
                aes(y = z_score, x = log2(ref_int), fill = Protein.Group), method = "loess", se = FALSE, color = "grey80") +
    geom_smooth(data = tmp,
                aes(y = z_score, x = log2(ref_int)), method = "loess", se = FALSE, color = "#0072B2", linewidth = 2) +
    scale_x_continuous(breaks = seq(20,30, by = 1)) +
    labs(x = "Intensity Alpha-1", y = "POI") +
    theme_classic()+
    theme(legend.position="none") +
    #scale_y_continuous(limits = c(-1,1)) +
    labs(title = i) -> plot
  
  i_name = str_replace_all(i, " / ", "")
  
  ggsave(plot, file = paste("../output/Figures/Figure_S4_", i_name, ".pdf", sep = ""), width = 5, height = 5)
  
}
