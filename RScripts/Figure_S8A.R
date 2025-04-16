## -- Figure S8A

d %>%
  gather(sample, int, !protein) %>%
  mutate(well = str_replace(sample, ".*__", "")) %>%
  mutate(int = ifelse(int == 0, NA, int)) %>%
  drop_na(int) %>%
  left_join(meta) %>%
  drop_na(slide) -> d_long

stats_mean <- d_long %>%
  group_by(well) %>%
  summarise(n = n()) %>%
  summarise(mean = mean(n)) %>%
  pull(mean)

stats_sd<- d_long %>%
  group_by(well) %>%
  summarise(n = n()) %>%
  summarise(sd = sd(n)) %>%
  pull(sd)

SA_incl <- d_long %>%
  group_by(well) %>%
  summarise(n = n()) %>%
  filter(n > (stats_mean - 0.5 * stats_sd)) %>%
  pull(well)

d_long %>%
  filter(well %in% SA_incl) %>%
  group_by(well) %>%
  summarise(n = n()) %>%
  summarise(mean = mean(n)) %>%
  pull(mean) -> stats_mean_incl

d_long %>%
  mutate(include = well %in% SA_incl) %>%
  group_by(well, slide, include) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(well) %>%
  mutate(rank = c(1:length(well))) %>%
  ggplot(aes(x = rank, y = n, color = slide, shape = include)) + 
  geom_point(size = 3)+
  theme_bw()+
  scale_y_continuous(limits = c(0,7000), breaks = seq(0,7000, by = 1000))+
  scale_color_manual(values = viridis(5)) +
  scale_shape_manual(values = c(4, 16))+
  labs(y = "# Proteins detected", x = "Run order")+
  theme_classic() +
  geom_hline(yintercept = stats_mean_incl, linetype = 2)

ggsave(file = "../output/Figures/Figure_S8A.pdf", width = 5, height = 5)

d_long -> d_long_ori

d_long %>%
  drop_na(int) %>%
  group_by(protein) %>%
  summarise(n = n()) -> pg_n_detected

pg_n_detected %>%
  filter(n > 0.3*length(SA_incl)) %>%
  pull(protein) -> pg_30

d_long %>% 
  filter(well %in% SA_incl) -> d_long
