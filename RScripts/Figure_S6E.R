## -- Figure 3B

polygon_to_map(xml = "../data/meta/xml_files/scDVP-S2.xml", data = d_norm, protein = "P01009", metadata = meta, subset = "Island")
ggsave("../output/Figures/Figure_S6E-S2-map.pdf", width = 5, height = 5)

polygon_to_map(xml = "../data/meta/xml_files/scDVP-S3.xml", data = d_norm, protein = "P01009", metadata = meta, subset = "Island")
ggsave("../output/Figures/Figure_S6E-S3-map.pdf", width = 5, height = 5)

ggplot(d_norm %>% filter(Protein.Group == "P01009", biology == "Island", sample_short == "scDVP-S2"), aes(x = sample_short, y = int_norm, fill = int_norm))+
  #geom_boxplot()+
  geom_jitter(width = 0.1, color = "black", size = 3, pch = 21)+
  scale_y_continuous(breaks = seq(0,40, by = 1))+
  theme_classic()+
  viridis::scale_fill_viridis() -> plot_alpha1_persistor

ggsave("../output/Figures/Figure_S6E-S2-boxplot.pdf", width = 3, height = 3)

ggplot(d_norm %>% filter(Protein.Group == "P01009", biology == "Island", sample_short == "scDVP-S3"), aes(x = sample_short, y = int_norm, fill = int_norm))+
  #geom_boxplot()+
  geom_jitter(width = 0.1, color = "black", size = 3, pch = 21)+
  scale_y_continuous(breaks = seq(0,40, by = 1))+
  theme_classic()+
  viridis::scale_fill_viridis() -> plot_alpha1_persistor

ggsave("../output/Figures/Figure_S6E-S3-boxplot.pdf", width = 3, height = 3)
