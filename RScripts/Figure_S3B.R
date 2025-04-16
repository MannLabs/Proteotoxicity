## -- Figure S3C

pg_n_detected %>%
  filter(n > 0.7*length(SA_incl)) %>%
  pull(Protein.Group) -> pg_70

limma_comparison <- function(level1, level2, coef){
  SA_stats <- meta_biopsies %>%
    filter(level %in% c(level1, level2)) %>%
    pull(well_id)
  
  d_long %>%
    filter(Protein.Group %in% pg_70) %>%
    filter(well_id %in% SA_stats) %>%
    mutate(int = log2(int)) %>%
    dplyr::select(- ms_id) %>%
    spread(well_id, int) %>%
    column_to_rownames("Protein.Group") -> d_stats
  
  patient <- as.factor((meta_biopsies %>% column_to_rownames("well_id"))[colnames(d_stats),]$sample)
  level <- as.factor((meta_biopsies %>% column_to_rownames("well_id"))[colnames(d_stats),]$level)
  
  design <- model.matrix(~ level + patient)
  
  fit <- lmFit(d_stats, design)
  fit <- eBayes(fit)
  limma <- topTable(fit, number = Inf, confint = TRUE, coef = coef, adjust.method = "fdr") %>%
    rownames_to_column("Protein.Group") %>%
    left_join(meta_pg) %>%
    mutate(comparison = coef)
  
  return(limma)
}

limma_comparison(level1 = "0_low", level2 = "1_moderate", coef = "level1_moderate") -> limma_0_1
limma_comparison(level1 = "1_moderate", level2 = "2_high", coef = "level2_high") -> limma_1_2

pg_significant_extremes <- unique(c(limma_level %>%
                                      filter(adj.P.Val < 0.05) %>%
                                      pull(Protein.Group),
                                    limma_0_1 %>%
                                      filter(adj.P.Val < 0.05) %>%
                                      pull(Protein.Group),
                                    limma_1_2 %>%
                                      filter(adj.P.Val < 0.05) %>%
                                      pull(Protein.Group)))

# Plot both categories against each other
limma_0_1 %>%
  dplyr::select(Protein.Group, logFC, adj.P.Val, B) %>%
  dplyr::rename(low_moderate_fc = logFC, low_moderate_p = adj.P.Val, low_moderate_B = B) -> limma_low_moderate

limma_1_2 %>%
  dplyr::select(Protein.Group, logFC, adj.P.Val, B) %>%
  dplyr::rename(moderate_high_fc = logFC, moderate_high_p = adj.P.Val, moderate_high_B = B) -> limma_moderate_high

full_join(limma_low_moderate, limma_moderate_high) %>%
  mutate(diff_fc = low_moderate_fc - moderate_high_fc) %>%
  mutate(diff_B = low_moderate_B - moderate_high_B) %>%
  left_join(meta_pg) %>%
  mutate(annotate = abs(diff_B) > 5) %>%
  mutate(significant = Protein.Group %in% pg_significant_extremes) -> limma_staging

ggplot()+
  geom_point(data = limma_staging %>% filter(significant == FALSE), aes(x = low_moderate_fc, y = moderate_high_fc), color = "grey40", alpha = 0.3)+
  geom_point(data = limma_staging %>% filter(significant == TRUE), aes(x = low_moderate_fc, y = moderate_high_fc), color = "black", alpha = 0.8)+
  geom_hline(yintercept = 0, lty = "dotted")+
  geom_vline(xintercept = 0, lty = "dotted") +
  geom_text_repel(data = limma_staging %>% filter(significant == TRUE),
                  aes(x = low_moderate_fc, y = moderate_high_fc, label = Genes), size = 1.3)+
  theme_classic()+
  scale_y_continuous(limits = c(-1,3))+
  scale_x_continuous(limits = c(-1,3))+
  scale_color_manual(values = c("grey70", "black"))+
  geom_abline(intercept = 0, slope = 1, lty = "dotted")

ggsave("../output/Figures/Figure_S3B.pdf", width = 5, height = 5)
