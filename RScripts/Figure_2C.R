## -- Figure 2C

plot_to_alpha <- function(poi, title){
  
  d_long %>%
    left_join(meta_pg) %>%
    left_join(meta_biopsies) %>%
    filter(Genes %in% c(poi)) %>%
    left_join(level_of_alpha1) %>%
    ggplot(aes(y = log2(int), x = log2(ref_int))) +
    geom_point()+
    #facet_grid(.~ Protein.Group)+
    theme_classic()+
    #coord_fixed(ratio = 1) +
    labs(x = "Intensity Alpha-1", y = poi, title = title) +
    geom_smooth(method = "loess", se = FALSE) +
    scale_x_continuous(breaks = seq(20,30, by = 1)) -> plot
  
  ggsave(plot, file = paste("../output/Figures/Figure_2C-", poi, ".pdf", sep = ""), width = 3, height = 3)
}

## Early up
plot_to_alpha(poi = "LGALS3BP", title = "early up")

## Late up
plot_to_alpha(poi = "DNAJB11", title = "late up")
plot_to_alpha(poi = "TNFSF10", title = "late up")

## Early down
plot_to_alpha(poi = "ASGR1", title = "early down")

## Late down
plot_to_alpha(poi = "SFPQ", title = "late down")