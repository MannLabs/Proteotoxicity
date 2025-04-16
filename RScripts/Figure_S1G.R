## Figure S1D

(
  ggplot(levels_of_alpha %>% drop_na(kleiner_score), aes(x = level, y = log10(int), fill = as.factor(kleiner_score)))+
  geom_boxplot(outlier.shape = NA) +
  #geom_line(aes(group = sample), alpha = 0.5) +
  #geom_point(pch = 21, size = 3) +
  scale_fill_manual(values = viridis(6, option = "magma")[5:2])+
  theme_bw()+
  scale_y_continuous(breaks = seq(6,9,by = 0.5)) +
  labs(x = "Levels of Alpha-1", y = "MS Intensity (log10)")
) %>% ggsave(., file = "../output/Figures/Figure_S1G.pdf", width = 7, height = 5)
