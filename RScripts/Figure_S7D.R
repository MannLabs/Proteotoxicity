## -- Figure S7D

ggplot(data = d_norm %>%
         left_join(meta_pg) %>%
         drop_na(spatial) %>%
         filter(Genes %in% c("ASS1", "HAL", "CYP2E1", "ADH1A")),
       aes(x = spatial, y = int_norm, fill = spatial))+
  
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1) +
  theme_classic()+
  scale_y_continuous(breaks = seq(from = 10, to = 40, by = 1))+
  labs(x = "Spatial region", y = "Intensity") +
  scale_fill_manual(values = c("white", "grey85", "grey60", "grey45")) +
  facet_grid(Genes ~ biology, scales = "free")

ggsave("../output/Figures/Figure_S7D.pdf", width = 8, height = 12)

d_norm %>%
  left_join(meta_pg) %>%
  drop_na(spatial) %>%
  filter(Genes %in% c("ASS1", "HAL", "CYP2E1", "ADH1A")) %>%
  group_by(spatial, biology, Genes) %>%
  summarise(n = n()) -> test
