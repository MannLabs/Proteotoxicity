## -- Figure S1D

d_long %>%
  filter(Protein.Group == "P01009") %>%
  left_join(meta_biopsies) -> levels_of_alpha


(
  ggplot(levels_of_alpha, aes(x = level, y = log10(int), fill = level))+
    geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = sample), alpha = 0.5) +
    geom_point(pch = 21, size = 3) +
    scale_fill_manual(values = viridis(5)[2:4])+
    theme_classic()+
    scale_y_continuous(breaks = seq(6,9,by = 0.5)) +
    labs(x = "Levels of Alpha-1", y = "MS Intensity (log10)")
) %>% ggsave(., file = "../output/Figures/Figure_S1D.pdf", width = 5, height = 5)