## -- Figure S2D

precursors <- read_tsv("../data/protein_quantification/Biopsies_PiZ_precursor-report.tsv")

precursors %>%
  gather(ms_id, int, contains("Evo12")) %>%
  mutate(int = as.numeric(int)) %>%
  filter(int != 0) %>%
  filter(ms_id %in% SA_incl) %>%
  mutate(well_id = str_replace(str_replace(ms_id, ".*__", ""), "_.*", "")) %>%
  left_join(meta_biopsies) %>%
  left_join(meta_pg) -> precursors_long

precursors_long %>%
  group_by(ion, Genes, level) %>%
  summarise(median = median(log2(int))) -> precursors_summarised

precursors_summarised %>%
  filter(Genes %in% c("SERPINA3", "SERPINA5", "SERPINC1", "SERPINA1")) %>%
  ggplot(aes(x = level, y = median, group = ion))+
  geom_point(size = 2, alpha = 2)+
  geom_line(lty = "dotted")+
  theme_classic()+
  facet_grid(.~ Genes)

ggsave(file = "../output/Figures/Figure_S2D.pdf", width = 7, height = 5)

# Check for redundancy
precursors_long %>%
  filter(Genes == "SERPINA1") %>%
  distinct(ion) %>%
  pull(ion) -> reporter_ions

precursors_long %>%
  filter(ion %in% reporter_ions) %>%
  distinct(Genes) %>%
  pull(Genes)
