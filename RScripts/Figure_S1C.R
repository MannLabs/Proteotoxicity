## -- Figure_S1C

d_long %>%
  left_join(meta_biopsies) %>%
  drop_na(kleiner_score) %>%
  group_by(Protein.Group) %>%
  mutate(cv = 100 * sd(int)/mean(int)) %>%
  ungroup() %>%
  ggplot(aes(x = as.factor(kleiner_score), y = cv)) +
  geom_violin(draw_quantiles = 0.5, fill = viridis(3)[2])+
  scale_y_continuous(limit = c(0,200)) +
  theme_classic() -> plot_cv

ggsave(plot_cv, file = "../output/Figures/Figure_S1C.pdf", width = 5, height = 5)