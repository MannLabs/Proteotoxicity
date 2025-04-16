## -- Figure S3D - E

meta_biopsies %>%
  filter(well_id %in% colnames(d_stats)) %>%
  column_to_rownames("well_id") %>%
  select(level) %>%
  arrange(level) -> meta_heatmap

reactome <- read_tsv("../data/meta/Reactome_pathway_of_interest.txt")

terms <- unique(reactome$`term name`)

for(i in terms){
  unlist(
    reactome%>%
      filter(`term name` == i) %>%
      pull(fg_protein) %>%
      str_split(";")
  ) -> candidates
  
  d_long %>%
    left_join(meta_pg) %>%
    filter(tolower(str_replace_all(Genes, ";.*", "")) %in% candidates) %>%
    select(Genes, int, well_id) %>%
    spread(well_id, int) %>%
    column_to_rownames("Genes") -> d_stats_subset

  pheatmap(t(scale(t(d_stats_subset[,rownames(meta_heatmap)]))),
           show_rownames = T, show_colnames = F,
           breaks = myBreaks, color = myColor,
           annotation_col = meta_heatmap,
           cluster_cols = F, border_color = NA, cellwidth = 8, cellheight = 8) -> plot_tmp
  
  filename = paste("../output/Figures/Figure_S3DEF_", str_replace_all(i, "%", "_"), ".pdf")
  
  ggsave(plot_tmp, file = filename, height = 10, width = 20)
}