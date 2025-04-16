## Figure S6C

d_long_noexcl %>%
  group_by(sample, batch) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(sample) %>%
  mutate(include = sample %in% SA_incl) %>%
  mutate(rank = c(1:nrow(d_stats_n))) %>%
  ggplot(aes(x = rank, y = n, shape = include, fill = batch))+
  geom_point()+
  labs(x = "Run number", y = "Number of proteins")+
  scale_shape_manual(values = c(4, 21))+
  theme_classic()+
  geom_hline(yintercept = c(2644, 2962), lty = "dotted")+
  scale_y_continuous(breaks = seq(from = 0, to = 4500, by = 500)) +
  scale_fill_manual(values = c("darkblue", "white")) -> plot_pg_run

ggsave("../output/Figures/Figure_S6C.pdf", width = 5, height = 4)
                     
                    