## -- Figure 2A

level_of_alpha1 <- d_long %>%
  filter(Protein.Group == "P01009") %>%
  dplyr::rename(ref_int = int) %>%
  dplyr::select(well_id, ref_int)

cor_alpha1 %>%
  arrange(cor) %>%
  top_n(11, cor) %>%
  pull(Protein.Group) -> cor_alpha1_top10

d_long %>%
  left_join(meta_pg) %>%
  filter(Protein.Group %in% cor_alpha1_top10) %>%
  dplyr::select(Genes, well_id, int) %>%
  left_join(level_of_alpha1) %>%
  group_by(Genes) %>%
  mutate(z_score = scale(int, scale = TRUE)) -> d_cor_top10

ggplot(d_cor_top10, aes(x = log2(ref_int), y = z_score, group = Genes))+
  geom_line(alpha = 0.4)+
  geom_smooth(se = F, color = viridis(3)[1])+
  theme_classic()+
  scale_x_continuous(breaks = seq(20,30, by = 1))

ggsave("../output/Figures/Figure_2A-linechart.pdf", width = 5, height = 5)

ggplot(level_of_alpha1, aes(y = 1, x = log2(ref_int)))+
  geom_boxplot()+
  geom_jitter(width = 0.01)+
  scale_y_continuous(limits = c(0,2))+
  scale_x_continuous(breaks = seq(20,30, by = 1))+
  theme_classic()

ggsave("../output/Figures/Figure_2A-boxplot.pdf", width = 5, height = 1)