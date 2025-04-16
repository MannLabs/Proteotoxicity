## -- Figure S1B

d_long_noexcl %>%
  group_by(ms_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(include = ms_id %in% too_low_n) %>%
  ggplot(aes(x = ms_id, y = n, pch = include, color = include))+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0,6000))+
  theme_classic()+
  geom_hline(yintercept = stats$mean, lty = "dotted", size = 1, color = "black")+
  geom_hline(yintercept = stats_lower_n, lty = "dotted", size = 1, color = "grey50")+
  scale_shape_manual(values = c(4, 21))+
  scale_color_manual(values = c("grey50", "blue")) +
  labs(x = "Run ID", y = "Number of proteins detected")+
  theme(axis.text.x=element_blank()) -> plot_n

ggsave(plot_n, file = "../output/Figures/Figure_S1b.pdf", width = 5, height = 5)