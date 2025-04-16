## -- Figure S6G

ggplot(data = d_norm %>% filter(Protein.Group == "P01009"), aes(x = spatial, y = int_norm, fill = spatial))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1) +
  theme_classic()+
  scale_y_continuous(breaks = seq(from = 12, to = 40, by = 1))+
  labs(x = "Spatial region", y = "Intensity") +
  scale_fill_manual(values = c("white", "grey85", "grey60", "grey45")) +
  geom_text_repel(aes(label = paste(sample_short, rank_id)))+
  facet_grid(.~ biology)

ggsave("../output/Figures/Figure_S6G.pdf", width = 6, height = 4)
