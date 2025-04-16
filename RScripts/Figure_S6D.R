## -- Figure S6D

d_stats_n %>%
  group_by(sample_short) %>%
  summarise(N = n()) -> d_stats_N

d_stats_n %>%
  left_join(d_stats_N) %>%
  filter(sample %in% SA_incl) %>%
  group_by(sample_short) %>%
  summarise(samples_included = n(), area_total = sum(area_corrected), area_mean = mean(area_corrected), protein = median(n))

# Precursors
pre %>%
  gather(sample, int, contains("OA")) %>%
  filter(int != 0) %>%
  filter(sample %in% SA_incl) %>%
  left_join(d_long %>% distinct(sample, sample_short)) %>%
  group_by(sample, sample_short) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(sample_short) %>%
  summarise(mean = median(n))
