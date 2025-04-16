## -- Figure S9

plot_cluster <- function(umap){
  ggplot()+
    geom_point(data = umap, aes(x = UMAP1, y = UMAP2, fill = cluster), pch = 21, alpha = 0.1)+
    geom_point(data = umap %>% filter(selected_cluster != 0) %>% distinct(selected_cluster, .keep_all = T), aes(x = UMAP1, y = UMAP2, color = selected_cluster), pch = 16, alpha = 1, size = 3)+
    coord_fixed(ratio = 1)+
    theme_classic()+
    scale_fill_viridis() +
    scale_color_viridis() +
    ggrepel::geom_text_repel(data = umap %>% filter(selected_cluster != 0) %>% distinct(selected_cluster, .keep_all = T), aes(x = UMAP1, y = UMAP2, label = selected_cluster), size = 3)
  
}

plot_cluster(umap_explant01)
plot_cluster(umap_explant01b)
plot_cluster(umap_explant03)
plot_cluster(umap_explant05)

d_long %>%
  dplyr::select(protein, well, int) %>%
  filter(protein %in% pg_30) %>%
  spread(well, int) %>%
  column_to_rownames("protein") -> d_stats
  
protein_of_interest <- "P50591"  # Replace with your protein name or index

# Extract the values of the protein of interest
protein_values <- as.numeric(d_stats[protein_of_interest, ])

# Compute correlations with all other proteins
correlations <- apply(d_stats, 1, function(row) cor(row, protein_values, use = "pairwise.complete.obs"))

# Convert correlations into a data frame for easier viewing
correlation_results <- data.frame(
  Protein.Group = rownames(d_stats),
  Correlation = correlations
) %>%
  left_join(pg)

plot_umap(patient = "Explant_DVP_01", poi = "O15304", umap = umap_explant01) #SIVA
plot_umap(patient = "Explant_DVP_01", poi = "O95831", umap = umap_explant01) #AIFM1
plot_umap(patient = "Explant_DVP_01", poi = "Q07812", umap = umap_explant01) #BAX
plot_umap(patient = "Explant_DVP_01", poi = "Q8IX12", umap = umap_explant01) #CCAR1
plot_umap(patient = "Explant_DVP_01", poi = "Q8N163", umap = umap_explant01) #CCAR2
plot_umap(patient = "Explant_DVP_01", poi = "Q96IZ0", umap = umap_explant01) #PAWR
plot_umap(patient = "Explant_DVP_01", poi = "Q9BZZ5", umap = umap_explant01) #API5
plot_umap(patient = "Explant_DVP_01", poi = "Q9ULZ3", umap = umap_explant01) #PYCARD



