## -- Figure S3C

#poi = c("hsa04146 Peroxisome", "hsa03040 Spliceosome", hsa03010 Ribosome, )

for(i in pathways_to_map){
  limma_staging %>%
    mutate(ENTREZ = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                          keys=Genes, 
                                          column="ENTREZID", 
                                          keytype="SYMBOL",
                                          multiVals="first")) %>%
    mutate(path_of_interest = ENTREZ %in% str_replace(db_kegg$kg.sets[[i]], ".*_", "")) -> limma_subset
  
  ggplot()+
    theme(
      axis.text = element_blank(),        # Remove axis labels
      axis.title = element_blank(),       # Remove axis titles
      axis.ticks = element_blank(),       # Remove axis ticks
      axis.line = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 1), # Black box
      panel.background = element_rect(fill = "white"))+
    geom_density_2d(data = limma_subset %>% filter(path_of_interest == FALSE), aes(x = low_moderate_fc, y = moderate_high_fc), color = "grey40", alpha = 0.3)+
    geom_point(data = limma_subset %>% filter(path_of_interest == TRUE), aes(x = low_moderate_fc, y = moderate_high_fc),
               fill = "#0072B2", color = "black", alpha = 1, pch = 21, size = 3)+
    geom_hline(yintercept = 0, lty = "dotted")+
    geom_vline(xintercept = 0, lty = "dotted") +
    scale_y_continuous(limits = c(-1,1))+
    scale_x_continuous(limits = c(-1,1))+
    geom_abline(intercept = 0, slope = 1, lty = "dotted")+
    labs(title = i)
  
  ggsave(file = paste("../output/Figures/Figure_S3C", i, ".pdf", sep = ""), width = 3, height = 3.2)
}