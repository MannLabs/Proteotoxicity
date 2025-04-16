## -- Figure 3A

polygon_to_map(xml = "../data/meta/xml_files/scDVP-S1.xml", data = d_norm, protein = "P01009", metadata = meta, subset = "Island")

ggsave("../output/Figures/Figure_3A-map.pdf", width = 5, height = 5)

ggplot(d_norm %>% filter(Protein.Group == "P01009", biology == "Island", sample_short == "scDVP-S1"), aes(x = sample_short, y = int_norm, fill = int_norm))+
  #geom_boxplot()+
  geom_jitter(width = 0.1, color = "black", size = 3, pch = 21)+
  scale_y_continuous(breaks = seq(0,40, by = 1))+
  theme_classic()+
  viridis::scale_fill_viridis() -> plot_alpha1_persistor

ggsave("../output/Figures/Figure_3A-boxplot.pdf", width = 3, height = 3)