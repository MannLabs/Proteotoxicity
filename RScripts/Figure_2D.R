## -- Figure 2D

nature_colors <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Green
  "#F0E442", # Yellow
  "#0072B2", # Dark Blue
  "#D55E00", # Vermillion
  "#CC79A7"  # Reddish-Purple
)

d_long %>%
  mutate(int = log2(int)) %>%
  dplyr::select(-well_id) %>%
  spread(ms_id, int) %>%
  column_to_rownames("Protein.Group")  -> d_wide

as.data.frame(t(scale(t(d_wide)))) %>%
  rownames_to_column("Protein.Group") %>%
  gather(ms_id, z_score, !Protein.Group) %>%
  left_join(d_long) %>%
  left_join(meta_pg) %>%
  left_join(meta_biopsies) %>%
  filter(Protein.Group %in% pg_30) %>%
  mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                          keys=Genes, 
                          column="ENSEMBL", 
                          keytype="SYMBOL",
                          multiVals="first")) %>%
  left_join(level_of_alpha1) %>%
  drop_na(ref_int) %>%
  filter(ENSEMBL != "NANA") -> d_scaled

db_kegg <- gage::kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)

pathways_to_map = c("hsa03008 Ribosome biogenesis in eukaryotes",
                    "hsa03030 DNA replication",
                    "hsa03040 Spliceosome",
                    "hsa04141 Protein processing in endoplasmic reticulum",
                    "hsa04146 Peroxisome",
                    "hsa04610 Complement and coagulation cascades")

db_kegg_specific <- c()
for(i in pathways_to_map){
  
  tmp <- data.frame(ENTREZID = str_replace(db_kegg$kg.sets[[i]], ".*_", "")) %>%
    mutate(ENSEMBL = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                            keys=ENTREZID, 
                            column="ENSEMBL", 
                            keytype="ENTREZID",
                            multiVals="first")) %>%
    mutate(pathway = i)
  db_kegg_specific <- rbind(db_kegg_specific, tmp)
}

db_kegg_specific %>%
  left_join(d_scaled %>%
              mutate(ENSEMBL = as.character(ENSEMBL)), relationship = "many-to-many") %>%
  drop_na(z_score) -> tmp

ggplot() +
  geom_smooth(data = tmp,
              aes(y = z_score, x = log2(ref_int), color = pathway), method = "loess", se = FALSE, linewidth = 1) +
  scale_x_continuous(breaks = seq(20,30, by = 1)) +
  labs(x = "Intensity Alpha-1", y = "z score") +
  theme_classic()+
  scale_color_manual(values = nature_colors)

ggsave(file = "../output/Figures/Figure_2D.pdf", width = 8, height = 5)
