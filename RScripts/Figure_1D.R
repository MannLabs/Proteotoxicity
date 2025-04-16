## -- Figure 1D

plot_protein_levels <- function(target){
  
  d_long %>%
    left_join(meta_pg) %>%
    filter(Genes == target) %>%
    left_join(meta_biopsies) %>%
    drop_na(kleiner_score) %>%
    ggplot(aes(x = level, y = log2(int), fill = level))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size = 1, width = 0.1)+
    scale_fill_manual(values = c("white", " grey80", "grey60"))+
    theme_bw()+
    labs(x = paste("Level of", target), y = "MS Intensity (log2)", title = target) -> plot_tmp
  
  ggsave(plot_tmp, file = paste("../output/Figures/Figure_1D_", target, ".pdf", sep = ""), width = 5, height = 5)
}

plot_protein_levels("SERPINA1")
plot_protein_levels("HSPA5")
plot_protein_levels(target = "LMAN1")
plot_protein_levels("TNFSF10")
plot_protein_levels("LGALS3BP")
plot_protein_levels(target = "TXNDC5")
