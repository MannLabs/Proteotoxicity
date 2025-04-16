## -- Figure S1F

(
  PCAtools::biplot(p,
                   x = 'PC2', y = 'PC3',
                   colby = 'int',
                   hline = 0, vline = 0,
                   labSize = 3,
                   lab = NA,
                   encircle = F,
                   encircleFill = F,
                   showLoadings = F)+
  scale_color_viridis()
) %>% ggsave(., file = "../output/Figures/Figure_S1F.pdf", width = 5, height = 5)
