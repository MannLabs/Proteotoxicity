## -- Figure 1E

limma_level %>%
  select(Genes, logFC) %>%
  write.table(., file = "../output/Tables/GSEA_Biopsies.rnk", col.names = F, row.names = F, quote = F, sep = "\t")

db_kegg <- gage::kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)

limma_level %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                          keys=Genes, 
                          column="ENSEMBL", 
                          keytype="SYMBOL",
                          multiVals="first")) %>%
  filter(ENSEMBL != "NANA") -> limma_level_ENSEMBL

dir.create("../output/Figures/Figure_1E/")

for(i in names(db_kegg$kg.sets)){
  db_kegg_specific <- data.frame(protein = str_replace(db_kegg$kg.sets[[i]], ".*_", "")) %>%
    mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                            keys=protein, 
                            column="ENSEMBL", 
                            keytype="ENTREZID",
                            multiVals="first"))
  
  limma_level_ENSEMBL %>%
    mutate(category = as.factor(ENSEMBL %in% db_kegg_specific$ENSEMBL)) -> limma_specific
  
  plot_data <- limma_specific %>%
    dplyr::select(ENSEMBL, logFC, category) %>%
    mutate(ymin = 0, ymax = 1)
  
  # Generate the plot
  ggplot(plot_data) +
    geom_vline(xintercept = (plot_data %>% filter(category == TRUE) %>% pull(logFC))) +
    theme_bw()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
    )+
    geom_vline(xintercept = 0, color = "blue", linewidth = 1)+
    scale_x_continuous(limits = c(-2,5)) -> plot
  
  top_hit <- limma_specific %>%
    filter(category == TRUE) %>%
    top_n(wt = logFC, n = 1) %>%
    pull(Genes)
  
  bottom_hit <- limma_specific %>%
    filter(category == TRUE) %>%
    top_n(wt = -logFC, n = 1) %>%
    pull(Genes)
  
  i_name = str_replace_all(i, " / ", "")
  
  ggsave(plot, file= paste("../output/Figures/Figure_1E/", i_name, "__t_", top_hit, "_b_",  bottom_hit, ".pdf", sep = ""), width = 5, height = 1)
}


