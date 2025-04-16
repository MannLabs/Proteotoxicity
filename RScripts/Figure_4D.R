## -- Figure 4D

umap_explant01 <- read_csv("../data/UMAP_data/Explant_DVP_01__UMAP_with_spectral_cluster_assignment.csv")
umap_explant01b <- read_csv("../data/UMAP_data/Explant_DVP_01b__UMAP_with_spectral_cluster_assignment.csv")
umap_explant03 <- read_csv("../data/UMAP_data/Explant_DVP_03__UMAP_with_spectral_cluster_assignment.csv")
umap_explant05 <- read_csv("../data/UMAP_data/Explant_DVP_05__UMAP_with_spectral_cluster_assignment.csv")

myBreaks <- c(seq(-2,2, by = 0.1))
myColor <- colorRampPalette(c('#FFFFFF', '#4D8365'))(length(myBreaks))

plot_umap <- function(patient, poi, umap){
  
  d_long %>%
    filter(slide == patient) %>%
    mutate(cluster = as.numeric(as.character(cluster))) %>%
    filter(protein == poi) %>%
    left_join(umap) -> umap_tstbot_poi
  
  ggplot()+
    geom_point(data = umap_tstbot_poi, aes(x = UMAP1, y = UMAP2, color = log2(int)), pch = 1, alpha = 0.1)+
    geom_point(data = umap_tstbot_poi %>% filter(selected_cluster != 0), aes(x = UMAP1, y = UMAP2, color = log2(int)), pch = 16, alpha = 1, size = 3)+
    coord_fixed(ratio = 1)+
    theme_classic()+
    scale_color_viridis()+
    scale_fill_viridis() +
    labs(title = poi) +
    theme(
      # Remove axis titles
      axis.title = element_blank(),
      # Remove axis text
      axis.text = element_blank(),
      # Remove axis ticks
      axis.ticks = element_blank(),
      # Remove panel grid
      panel.grid = element_blank(),
      # Remove the legend
      legend.position = "none",
      axis.line = element_blank()) -> plot
  
  ggsave(plot, file = paste("../output/Figures/Figure_4D_", patient, "_", poi, ".png"), width = 6, height = 5)
  return(plot)
  
}

plot_umap(patient = "Explant_DVP_01", poi = "P01009", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "P02741", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "Q96HE7", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "P50591", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "O75795", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "Q9UHF1", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "P27797", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "P30040", umap = umap_explant01)

plot_umap(patient = "Explant_DVP_01", poi = "P49257", umap = umap_explant01)
plot_umap(patient = "Explant_DVP_01", poi = "O15254", umap = umap_explant01)

plot_umap(patient = "Explant_DVP_01b", poi = "P01009", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "P02741", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "Q96HE7", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "P50591", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "O75795", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "Q9UHF1", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "P27797", umap = umap_explant01b)
plot_umap(patient = "Explant_DVP_01b", poi = "P30040", umap = umap_explant01b)

plot_umap(patient = "Explant_DVP_03", poi = "P01009", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "P02741", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "Q96HE7", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "P50591", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "O75795", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "Q9UHF1", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "P27797", umap = umap_explant03)
plot_umap(patient = "Explant_DVP_03", poi = "P30040", umap = umap_explant03)

plot_umap(patient = "Explant_DVP_05", poi = "P01009", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "P02741", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "Q96HE7", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "P50591", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "O75795", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "Q9UHF1", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "P27797", umap = umap_explant05)
plot_umap(patient = "Explant_DVP_05", poi = "P30040", umap = umap_explant05)
