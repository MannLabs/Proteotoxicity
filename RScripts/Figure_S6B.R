## --Figure S6B

ext_stats <- read_tsv("../data/protein_quantification/scDVP_windows.tsv") %>%
  mutate(Method = factor(Method, levels = c("8Th18ms", "20Th40ms", "30vW-40ms", "45vW-28ms"))) %>%
  filter(Included == "yes")

ext_summary <- ext_stats %>%
  group_by(Method) %>%
  summarise(mean_protein = mean(`Protein Groups`), sd_protein = sd(`Protein Groups`),
            mean_precursor = mean(Precursor), sd_precursor = sd(Precursor))

ggplot()+
  geom_bar(data = ext_summary, aes(x = Method, y = mean_protein), stat = "identity") +
  geom_jitter(data = ext_stats, aes(x = Method, y = `Protein Groups`), width = 0.3, alpha = 0.7)+
  geom_errorbar(data = ext_summary, aes(x = Method, ymin = mean_protein, ymax = mean_protein + sd_protein), 
                width = 0.2)+
  theme_classic()+
  scale_y_continuous(breaks = seq(0,3000, by = 500))

ggsave('../output/Figures/Figure_S6b.pdf', width = 5, height = 4)
