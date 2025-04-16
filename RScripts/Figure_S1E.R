## -- Figure S1E

d_long %>%
  mutate(int = log2(int)) %>%
  select(-well_id) %>%
  spread(ms_id, int) %>%
  column_to_rownames("Protein.Group") %>%
  filter(complete.cases(.)) -> d_complete

d_long %>%
  filter(Protein.Group == "P01009") %>%
  select(well_id, ms_id, int) %>%
  mutate(int = log2(int)) %>%
  distinct(ms_id, .keep_all = T) %>%
  left_join(meta_biopsies) %>%
  drop_na(ms_id) %>%
  column_to_rownames("ms_id") -> meta_pca

p <- PCAtools::pca(mat = d_complete, metadata = meta_pca[colnames(d_complete),], removeVar = 0.1)

(
  PCAtools::biplot(p ,
                   colby = 'kleiner_score',
                   hline = 0, vline = 0,
                   labSize = 3,
                   lab = NA,
                   encircle = F,
                   #legendPosition = NA,
                   encircleFill = F,
                   showLoadings = F)+
    scale_color_viridis(option = "E")
) %>% ggsave(., file = "../output/Figures/Figure_S1E.pdf", width = 3, height = 5)
